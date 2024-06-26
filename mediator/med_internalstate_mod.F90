module med_internalstate_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use ESMF         , only : ESMF_RouteHandle, ESMF_FieldBundle, ESMF_State, ESMF_Field, ESMF_VM
  use ESMF         , only : ESMF_GridComp, ESMF_Mesh, ESMF_MAXSTR, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
  use med_kind_mod , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_utils_mod, only : chkerr => med_utils_ChkErr

  implicit none
  private

  ! public routines
  public :: med_internalstate_init
  public :: med_internalstate_coupling
  public :: med_internalstate_defaultmasks

  integer, public :: logunit            ! logunit for mediator log output
  integer, public :: diagunit           ! diagunit for budget output (med main only)
  logical, public :: maintask = .false. ! is this the maintask
  integer, public :: med_id             ! needed currently in med_io_mod and set in esm.F90

  ! Components
  integer, public :: compmed = 1
  integer, public :: compatm = 2
  integer, public :: complnd = 3
  integer, public :: compocn = 4
  integer, public :: compice = 5
  integer, public :: comprof = 6
  integer, public :: compwav = 7
  integer, public :: ncomps =  7 ! this will be incremented if the size of compglc is > 0
  integer, public, allocatable :: compglc(:)

  ! Generic component name (e.g. atm, ocn...)
  character(len=CS), public, allocatable :: compname(:)

  ! Specific component name (e.g. datm, mom6, etc...)
  character(len=CS), public :: med_name = ''
  character(len=CS), public :: atm_name = ''
  character(len=CS), public :: lnd_name = ''
  character(len=CS), public :: ocn_name = ''
  character(len=CS), public :: ice_name = ''
  character(len=CS), public :: rof_name = ''
  character(len=CS), public :: wav_name = ''
  character(len=CS), public :: glc_name = ''

  ! Coupling mode
  character(len=CS), public :: coupling_mode ! valid values are [cesm,ufs.nfrac,ufs.frac,ufs.nfrac.aoflux,ufs.frac.aoflux,hafs,hafs.mom6]

  ! Atmosphere-ocean flux algorithm
  character(len=CS), public :: aoflux_code   ! valid values are [cesm,ccpp]

  ! Atmosphere-ocean CCPP suite name
  character(len=CL), public :: aoflux_ccpp_suite

  ! Default src and destination masks for mapping
  integer, public, allocatable :: defaultMasks(:,:)

  ! Mapping
  integer , public, parameter :: mapunset          = 0
  integer , public, parameter :: mapbilnr          = 1
  integer , public, parameter :: mapconsf          = 2
  integer , public, parameter :: mapconsd          = 3
  integer , public, parameter :: mappatch          = 4
  integer , public, parameter :: mapfcopy          = 5
  integer , public, parameter :: mapnstod          = 6  ! nearest source to destination
  integer , public, parameter :: mapnstod_consd    = 7  ! nearest source to destination followed by conservative dst
  integer , public, parameter :: mapnstod_consf    = 8  ! nearest source to destination followed by conservative frac
  integer , public, parameter :: mappatch_uv3d     = 9  ! rotate u,v to 3d cartesian space, map from src->dest, then rotate back
  integer , public, parameter :: mapbilnr_uv3d     = 10 ! rotate u,v to 3d cartesian space, map from src->dest, then rotate back
  integer , public, parameter :: map_rof2ocn_ice   = 11 ! custom smoothing map to map ice from rof->ocn (cesm only)
  integer , public, parameter :: map_rof2ocn_liq   = 12 ! custom smoothing map to map liq from rof->ocn (cesm only)
  integer , public, parameter :: map_glc2ocn_liq   = 13 ! custom smoothing map to map liq from glc->ocn (cesm only)
  integer , public, parameter :: map_glc2ocn_ice   = 14 ! custom smoothing map to map ice from glc->ocn (cesm only)
  integer , public, parameter :: mapfillv_bilnr    = 15 ! fill value followed by bilinear
  integer , public, parameter :: mapbilnr_nstod    = 16 ! bilinear with nstod extrapolation
  integer , public, parameter :: mapconsf_aofrac   = 17 ! conservative with aofrac normalization (ufs only)
  integer , public, parameter :: nmappers          = 17
  character(len=*) , public, parameter :: mapnames(nmappers) = &
       (/'bilnr       ',&
         'consf       ',&
         'consd       ',&
         'patch       ',&
         'fcopy       ',&
         'nstod       ',&
         'nstod_consd ',&
         'nstod_consf ',&
         'patch_uv3d  ',&
         'bilnr_uv3d  ',&
         'rof2ocn_ice ',&
         'rof2ocn_liq ',&
         'glc2ocn_ice ',&
         'glc2ocn_liq ',&
         'fillv_bilnr ',&
         'bilnr_nstod ',&
         'consf_aofrac'/)

  type, public :: packed_data_type
     integer, allocatable :: fldindex(:) ! size of number of packed fields
     character(len=CS)    :: mapnorm     ! normalization for packed field
     type(ESMF_Field)     :: field_src    ! packed sourced field
     type(ESMF_Field)     :: field_dst    ! packed destination field
     type(ESMF_Field)     :: field_fracsrc
     type(ESMF_Field)     :: field_fracdst
  end type packed_data_type

  logical, public :: dststatus_print = .false.

  ! Mesh info
  type, public ::  mesh_info_type
     real(r8), pointer :: areas(:) => null()
     real(r8), pointer :: lats(:) => null()
     real(r8), pointer :: lons(:) => null()
  end type mesh_info_type

  logical          , public :: samegrid_atmlnd = .true. ! true=>atm and lnd are on the same grid
  character(len=CS), public :: mrg_fracname_lnd2atm_state
  character(len=CS), public :: mrg_fracname_lnd2atm_flux
  character(len=CS), public :: map_fracname_lnd2atm
  character(len=CS), public :: mrg_fracname_lnd2rof
  character(len=CS), public :: map_fracname_lnd2rof
  character(len=CS), public :: mrg_fracname_lnd2glc
  character(len=CS), public :: map_fracname_lnd2glc

  ! private internal state to keep instance data
  type InternalStateStruct

    ! Present/allowed coupling/active coupling logical flags
    logical, pointer :: comp_present(:)               ! comp present flag
    logical, pointer :: med_coupling_active(:,:)      ! computes the active coupling
    logical, pointer :: med_data_active(:,:)          ! uses stream data to provide background fill
    logical, pointer :: med_data_force_first(:)       ! force to use stream data for first coupling timestep
    integer          :: num_icesheets                 ! obtained from attribute
    logical          :: ocn2glc_coupling = .false.    ! obtained from attribute
    logical          :: lnd2glc_coupling = .false.
    logical          :: accum_lnd2glc = .false.

    ! Mediator vm
    type(ESMF_VM) :: vm

    ! Global nx,ny dimensions of input arrays (needed for mediator history output)
    integer, pointer   :: nx(:), ny(:)
    ! Number of nx*ny domains (needed for cubed-sphere and regional domains)
    integer, pointer   :: ntile(:)

    ! Import/Export Scalars
    character(len=CL) :: flds_scalar_name = ''
    integer           :: flds_scalar_num = 0
    integer           :: flds_scalar_index_nx = 0
    integer           :: flds_scalar_index_ny = 0
    integer           :: flds_scalar_index_ntile = 0
    integer           :: flds_scalar_index_nextsw_cday = 0
    integer           :: flds_scalar_index_precip_factor = 0
    real(r8)          :: flds_scalar_precip_factor = 1._r8  ! actual value of precip factor from ocn

    ! NState_Imp and NState_Exp are the standard NUOPC coupling datatypes
    ! FBImp and FBExp are the internal mediator datatypes
    ! NState_Exp(n) = FBExp(n), copied in the connector prep phase
    ! FBImp(n,n) = NState_Imp(n), copied in connector post phase
    ! FBImp(n,k) is the FBImp(n,n) interpolated to grid k
    ! Import/export States and field bundles (the field bundles have the scalar fields removed)
    type(ESMF_State)       , pointer :: NStateImp(:)   ! Import data from various component, on their grid
    type(ESMF_State)       , pointer :: NStateExp(:)   ! Export data to various component, on their grid
    type(ESMF_FieldBundle) , pointer :: FBImp(:,:)     ! Import data from various components interpolated to various grids
    type(ESMF_FieldBundle) , pointer :: FBExp(:)       ! Export data for various components, on their grid

    ! Mediator field bundles for ocean albedo
    type(ESMF_FieldBundle) :: FBMed_ocnalb_o            ! Ocn albedo on ocn grid
    type(ESMF_FieldBundle) :: FBMed_ocnalb_a            ! Ocn albedo on atm grid
    type(packed_data_type), pointer :: packed_data_ocnalb_o2a(:) ! packed data for mapping ocn->atm

    ! Mediator field bundles and other info for atm/ocn flux computation
    character(len=CS)      :: aoflux_grid                        ! 'ogrid', 'agrid' or 'xgrid'
    type(ESMF_Mesh)        :: aoflux_mesh                        ! Mesh used for atm/ocn flux computation
    type(ESMF_FieldBundle) :: FBMed_aoflux_a                     ! Ocn/Atm flux output fields on atm grid
    type(ESMF_FieldBundle) :: FBMed_aoflux_o                     ! Ocn/Atm flux output fields on ocn grid
    type(packed_data_type), pointer :: packed_data_aoflux_o2a(:) ! packed data for mapping ocn->atm

    ! Mapping
    ! RH(n,k,m) is a RH from grid n to grid k, map type m
    type(ESMF_RouteHandle) , pointer :: RH(:,:,:)            ! Routehandles for pairs of components and different mappers
    type(ESMF_Field)       , pointer :: field_NormOne(:,:,:) ! Unity static normalization
    type(packed_data_type) , pointer :: packed_data(:,:,:)   ! Packed data structure needed to efficiently map field bundles

    ! Fractions
    type(ESMF_FieldBundle), pointer :: FBfrac(:)     ! Fraction data for various components, on their grid

    ! Data
    type(ESMF_FieldBundle) , pointer :: FBData(:)    ! Background data for various components, on their grid, provided by CDEPS inline

    ! Accumulators for export field bundles
    type(ESMF_FieldBundle) :: FBExpAccumOcn      ! Accumulator for Ocn export on Ocn grid
    integer                :: ExpAccumOcnCnt = 0 ! Accumulator counter for FBExpAccumOcn
    type(ESMF_FieldBundle) :: FBExpAccumWav      ! Accumulator for Wav export on Wav grid
    integer                :: ExpAccumWavCnt = 0 ! Accumulator counter for FBExpAccumWav

    ! Component Mesh info
    type(mesh_info_type)   , pointer :: mesh_info(:)
    type(ESMF_FieldBundle) , pointer :: FBArea(:)     ! needed for mediator history writes

  end type InternalStateStruct

  type, public :: InternalState
    type(InternalStateStruct), pointer :: wrap
  end type InternalState

  character(len=*), parameter :: u_FILE_u  = &
       __FILE__

!=====================================================================
contains
!=====================================================================

  subroutine med_internalstate_init(gcomp, rc)

    use ESMF              , only : ESMF_LogFoundAllocError, ESMF_AttributeGet
    use NUOPC_Comp        , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    integer          , intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    logical                    :: ispresent, isset
    integer                    :: n, ns, n1
    character(len=8)           :: cnum
    character(len=CS)          :: cvalue
    character(len=ESMF_MAXSTR) :: mesh_glc
    character(len=CX)          :: msgString
    character(len=3)           :: name
    integer                    :: num_icesheets
    character(len=CL)          :: atm_mesh_name
    character(len=CL)          :: lnd_mesh_name
    logical                    :: isPresent_lnd, isSet_lnd
    logical                    :: isPresent_atm, isSet_atm
    character(len=*),parameter :: subname=' (internalstate init) '
    !-----------------------------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine if atm and lnd have the same mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=atm_mesh_name, &
         isPresent=isPresent_atm, isSet=isSet_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='mesh_lnd', value=lnd_mesh_name, &
         isPresent=isPresent_lnd, isSet=isSet_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ((isPresent_lnd .and. isSet_lnd) .and. (isPresent_atm .and. isSet_atm)) then
      if (trim(atm_mesh_name) == trim(lnd_mesh_name)) then
        samegrid_atmlnd = .true.
      else
        samegrid_atmlnd = .false.
      end if
    else
      samegrid_atmlnd = .true.
    end if

    ! See med_fraction_mod for the following definitions
    if (samegrid_atmlnd) then
      map_fracname_lnd2atm       = 'lfrin' ! in fraclist_a
      mrg_fracname_lnd2atm_state = 'lfrac' ! in fraclist_a
      mrg_fracname_lnd2atm_flux  = 'lfrac' ! in fraclist_a
      map_fracname_lnd2rof       = 'lfrac' ! in fraclist_r
      mrg_fracname_lnd2rof       = 'lfrac' ! in fraclist_r
      map_fracname_lnd2glc       = 'lfrac' ! in fraclist_g
      mrg_fracname_lnd2glc       = 'lfrac' ! in fraclist_g
    else
      map_fracname_lnd2atm       = 'lfrin' ! in fraclist_a
      mrg_fracname_lnd2atm_state = 'lfrac' ! in fraclist_a
      mrg_fracname_lnd2atm_flux  = 'lfrin' ! in fraclist_a
      map_fracname_lnd2rof       = 'lfrin' ! in fraclist_r
      mrg_fracname_lnd2rof       = 'lfrin' ! in fraclist_r
      map_fracname_lnd2glc       = 'lfrin' ! in fraclist_g
      mrg_fracname_lnd2rof       = 'lfrin' ! in fraclist_g
    endif

    if (maintask) then
      write(logunit,'(a,i8)') trim(subname)//'      map_fracname_lnd2atm       = '//trim(map_fracname_lnd2atm)      //' in fraclist_a'
      write(logunit,'(a,i8)') trim(subname)//'      mrg_fracname_lnd2atm_state = '//trim(mrg_fracname_lnd2atm_state)//' in fraclist_a'
      write(logunit,'(a,i8)') trim(subname)//'      mrg_fracname_lnd2atm_flux  = '//trim(mrg_fracname_lnd2atm_flux) //' in fraclist_a'
      write(logunit,'(a,i8)') trim(subname)//'      map_fracname_lnd2rof       = '//trim(map_fracname_lnd2rof)      //' in fraclist_r'
      write(logunit,'(a,i8)') trim(subname)//'      mrg_fracname_lnd2rof       = '//trim(mrg_fracname_lnd2rof)      //' in fraclist_r'
      write(logunit,'(a,i8)') trim(subname)//'      map_fracname_lnd2glc       = '//trim(map_fracname_lnd2glc)      //' in fraclist_g'
      write(logunit,'(a,i8)') trim(subname)//'      mrg_fracname_lnd2rof       = '//trim(mrg_fracname_lnd2rof)      //' in fraclist_g'
    end if

    ! Determine if glc is present
    call NUOPC_CompAttributeGet(gcomp, name='GLC_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    num_icesheets = 0
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sglc') then
          call NUOPC_CompAttributeGet(gcomp, name='mesh_glc', value=mesh_glc, isPresent=isPresent, isSet=isSet, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          glc_name = trim(cvalue)
          if (isPresent .and. isSet) then
             ! determine number of ice sheets - search in mesh_glc for colon deliminted strings
             if (len_trim(cvalue) > 0) then
                do n = 1, len_trim(mesh_glc)
                   if (mesh_glc(n:n) == ':') num_icesheets = num_icesheets + 1
                end do
                num_icesheets = num_icesheets + 1
             endif
             if (maintask) then
                write(logunit,'(a,i8)') trim(subname)//' number of ice sheets is ',num_icesheets
             end if
          end if
          ! now determing the number of multiple ice sheets and increment ncomps accordingly
          allocate(compglc(num_icesheets))
          compglc(:) = 0
          do ns = 1,num_icesheets
             ncomps = ncomps + 1
             compglc(ns) = ncomps
          end do
       end if
    end if

    ! Determine present flags starting with glc component
    allocate(is_local%wrap%comp_present(ncomps))
    is_local%wrap%comp_present(:) = .false.
    if (num_icesheets > 0) then
       do ns = 1,num_icesheets
          is_local%wrap%comp_present(compglc(ns)) = .true.
       end do
    end if
    is_local%wrap%num_icesheets = num_icesheets
    call NUOPC_CompAttributeGet(gcomp, name='mediator_present', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) is_local%wrap%comp_present(compmed)
    end if
    call NUOPC_CompAttributeGet(gcomp, name='MED_model', value=med_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', value=atm_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(atm_name) /= 'satm') is_local%wrap%comp_present(compatm) = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='LND_model', value=lnd_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(lnd_name) /= 'slnd') is_local%wrap%comp_present(complnd) = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='OCN_model', value=ocn_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(ocn_name) /= 'socn') is_local%wrap%comp_present(compocn) = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ICE_model', value=ice_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(ice_name) /= 'sice') is_local%wrap%comp_present(compice) = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ROF_model', value=rof_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(rof_name) /= 'srof') is_local%wrap%comp_present(comprof) = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='WAV_model', value=wav_name, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(wav_name) /= 'swav') is_local%wrap%comp_present(compwav) = .true.
    end if

    ! Allocate memory now that ncomps is determined
    allocate(is_local%wrap%med_coupling_active(ncomps,ncomps))
    allocate(is_local%wrap%med_data_active(ncomps,ncomps))
    allocate(is_local%wrap%med_data_force_first(ncomps))
    allocate(is_local%wrap%nx(ncomps))
    allocate(is_local%wrap%ny(ncomps))
    allocate(is_local%wrap%ntile(ncomps))
    allocate(is_local%wrap%NStateImp(ncomps))
    allocate(is_local%wrap%NStateExp(ncomps))
    allocate(is_local%wrap%FBImp(ncomps,ncomps))
    allocate(is_local%wrap%FBExp(ncomps))
    allocate(is_local%wrap%packed_data_ocnalb_o2a(nmappers))
    allocate(is_local%wrap%packed_data_aoflux_o2a(nmappers))
    allocate(is_local%wrap%RH(ncomps,ncomps,nmappers))
    allocate(is_local%wrap%field_NormOne(ncomps,ncomps,nmappers))
    allocate(is_local%wrap%packed_data(ncomps,ncomps,nmappers))
    allocate(is_local%wrap%FBfrac(ncomps))
    allocate(is_local%wrap%FBArea(ncomps))
    allocate(is_local%wrap%FBData(ncomps))
    allocate(is_local%wrap%mesh_info(ncomps))

    ! Determine component names
    allocate(compname(ncomps))
    compname(compmed) = 'med'
    compname(compatm) = 'atm'
    compname(complnd) = 'lnd'
    compname(compocn) = 'ocn'
    compname(compice) = 'ice'
    compname(comprof) = 'rof'
    compname(compwav) = 'wav'
    do ns = 1,is_local%wrap%num_icesheets
       write(cnum,'(i0)') ns
       compname(compglc(ns)) = 'glc' // trim(cnum)
    end do

    if (maintask) then
       ! Write out present flags
       write(logunit,*)
       do n1 = 1,ncomps
          name = trim(compname(n1))  ! this trims the ice sheets index from the glc name
          write(msgString,'(A,L4)') trim(subname)//' comp_present(comp'//name//') = ',&
               is_local%wrap%comp_present(n1)
          write(logunit,'(a)') trim(msgString)
       end do

       ! Write out model names if they are present
       write(logunit,*)
       if (is_local%wrap%comp_present(compatm)) write(logunit,'(a)') trim(subname) // " atm model= "//trim(atm_name)
       if (is_local%wrap%comp_present(complnd)) write(logunit,'(a)') trim(subname) // " lnd model= "//trim(lnd_name)
       if (is_local%wrap%comp_present(compocn)) write(logunit,'(a)') trim(subname) // " ocn model= "//trim(ocn_name)
       if (is_local%wrap%comp_present(compice)) write(logunit,'(a)') trim(subname) // " ice model= "//trim(ice_name)
       if (is_local%wrap%comp_present(comprof)) write(logunit,'(a)') trim(subname) // " rof model= "//trim(rof_name)
       if (is_local%wrap%comp_present(compwav)) write(logunit,'(a)') trim(subname) // " wav model= "//trim(wav_name)
       if (is_local%wrap%comp_present(compmed)) write(logunit,'(a)') trim(subname) // " med model= "//trim(med_name)
       if (is_local%wrap%num_icesheets > 0) then
          if (is_local%wrap%comp_present(compglc(1))) write(logunit,'(a)') trim(subname) // " glc model= "//trim(glc_name)
       end if
       write(logunit,*)
    end if

    ! Obtain dststatus_print setting if present
    call NUOPC_CompAttributeGet(gcomp, name='dststatus_print', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) dststatus_print=(trim(cvalue) == "true")
    write(msgString,*) trim(subname)//': Mediator dststatus_print is ',dststatus_print
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! Initialize flag for background fill using data
    is_local%wrap%med_data_active(:,:) = .false.
    is_local%wrap%med_data_active(compocn,compatm) = .true.
    is_local%wrap%med_data_active(compatm,compocn) = .true.
    is_local%wrap%med_data_active(compatm,compwav) = .true.

    ! Initialize flag to force using data in first coupling time step
    is_local%wrap%med_data_force_first(:) = .false.

  end subroutine med_internalstate_init

  !=====================================================================
  subroutine med_internalstate_coupling(gcomp, rc)

    !----------------------------------------------------------
    ! Check for active coupling interactions
    ! must be allowed, bundles created, and both sides have some fields
    ! This is called from med.F90 in the DataInitialize routine
    !----------------------------------------------------------

    use ESMF            , only : ESMF_StateIsCreated
    use NUOPC           , only : NUOPC_CompAttributeGet
    use med_methods_mod , only : State_getNumFields => med_methods_State_getNumFields

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)  :: is_local
    integer              :: n1, n2, ns
    integer              :: cntn1, cntn2
    logical, allocatable :: med_coupling_allowed(:,:)
    character(len=CL)    :: cvalue
    character(len=CX)    :: msgString
    logical              :: isPresent, isSet
    character(len=*),parameter :: subname=' (internalstate allowed coupling) '
    !-----------------------------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! This defines the med_coupling_allowed a starting point for what is
    ! allowed in this coupled system.  It will be revised further after the system
    ! starts, but any coupling set to false will never be allowed.
    ! are allowed, just update the table below.

    if (maintask) then
       write(logunit,'(a)') trim(subname) // "Initializing active coupling flags"
    end if

    ! Initialize med_coupling_allowed
    allocate(med_coupling_allowed(ncomps,ncomps))
    med_coupling_allowed(:,:) = .false.
    is_local%wrap%med_coupling_active(:,:) = .false.

    ! to atmosphere
    med_coupling_allowed(complnd,compatm) = .true.
    med_coupling_allowed(compice,compatm) = .true.
    med_coupling_allowed(compocn,compatm) = .true.
    med_coupling_allowed(compwav,compatm) = .true.

    ! to land
    med_coupling_allowed(compatm,complnd) = .true.
    med_coupling_allowed(comprof,complnd) = .true.
    do ns = 1,is_local%wrap%num_icesheets
       med_coupling_allowed(compglc(ns),complnd) = .true.
    end do

    ! to ocean
    med_coupling_allowed(compatm,compocn) = .true.
    med_coupling_allowed(compice,compocn) = .true.
    med_coupling_allowed(comprof,compocn) = .true.
    med_coupling_allowed(compwav,compocn) = .true.

    ! to ice
    med_coupling_allowed(compatm,compice) = .true.
    med_coupling_allowed(compocn,compice) = .true.
    med_coupling_allowed(comprof,compice) = .true.
    med_coupling_allowed(compwav,compice) = .true.
    do ns = 1,is_local%wrap%num_icesheets
       med_coupling_allowed(compglc(ns),compice) = .true.
    end do

    ! to river
    med_coupling_allowed(complnd,comprof) = .true.
    do ns = 1,is_local%wrap%num_icesheets
       med_coupling_allowed(compglc(ns),comprof) = .true.
    end do

    ! to wave
    med_coupling_allowed(compatm,compwav) = .true.
    med_coupling_allowed(compocn,compwav) = .true.
    med_coupling_allowed(compice,compwav) = .true.

    ! to land-ice
    call NUOPC_CompAttributeGet(gcomp, name='ocn2glc_coupling', value=cvalue, &
            isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       ! multiple ocean depths for temperature and salinity sent from the ocn to glc
       read(cvalue,*) is_local%wrap%ocn2glc_coupling
    else
       is_local%wrap%ocn2glc_coupling = .false.
    end if
    do ns = 1,is_local%wrap%num_icesheets
       med_coupling_allowed(complnd,compglc(ns)) = .true.
       med_coupling_allowed(compocn,compglc(ns)) = is_local%wrap%ocn2glc_coupling
    end do

    ! initialize med_coupling_active table
    is_local%wrap%med_coupling_active(:,:) = .false.
    do n1 = 1,ncomps
       if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call State_GetNumFields(is_local%wrap%NStateImp(n1), cntn1, rc=rc) ! Import Field Count
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (cntn1 > 0) then
             do n2 = 1,ncomps
                if (is_local%wrap%comp_present(n2) .and. ESMF_StateIsCreated(is_local%wrap%NStateExp(n2),rc=rc) .and. &
                     med_coupling_allowed(n1,n2)) then
                   call State_GetNumFields(is_local%wrap%NStateExp(n2), cntn2, rc=rc) ! Import Field Count
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   if (cntn2 > 0) is_local%wrap%med_coupling_active(n1,n2) = .true.
                endif
             enddo
          end if
       endif
    enddo

    ! create tables of allowed and active coupling flags
    ! - the rows are the destination of coupling
    ! - the columns are the source of coupling
    ! - So, the second column indicates which models the atm is coupled to.
    ! - And the second row indicates which models are coupled to the atm.
    if (maintask) then
       write(logunit,*) ' '
       write(logunit,'(A)') trim(subname)//' Allowed coupling flags'
       write(logunit,'(2x,A10,20(A5))') '|from to -> ',(compname(n2),n2=1,ncomps)
       do n1 = 1,ncomps
          write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)), &
               (med_coupling_allowed(n1,n2),n2=1,ncomps)
          do n2 = 1,len_trim(msgString)
             if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
          enddo
          write(logunit,'(A)') trim(msgString)
       enddo

       write(logunit,*) ' '
       write(logunit,'(A)') subname//' Active coupling flags'
       write(logunit,'(2x,A10,20(A5))') '|from to -> ',(compname(n2),n2=1,ncomps)
       do n1 = 1,ncomps
          write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)), &
               (is_local%wrap%med_coupling_active(n1,n2),n2=1,ncomps)
          do n2 = 1,len_trim(msgString)
             if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
          enddo
          write(logunit,'(A)') trim(msgString)
       enddo
       write(logunit,*) ' '
    endif

    ! Determine lnd2glc_coupling flag
    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(complnd,compglc(ns))) then
          is_local%wrap%lnd2glc_coupling = .true.
          exit
       end if
    end do

    ! Determine accum_lnd2glc flag
    if (is_local%wrap%lnd2glc_coupling) then
       is_local%wrap%accum_lnd2glc = .true.
    else
       ! Determine if will create auxiliary history file that contains
       ! lnd2glc data averaged over the year
       call NUOPC_CompAttributeGet(gcomp, name="histaux_l2x1yrg", value=cvalue, &
            isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          read(cvalue,*) is_local%wrap%accum_lnd2glc
       end if
    end if

    ! Determine ocn2glc_coupling flag
    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(compocn,compglc(ns))) then
          is_local%wrap%ocn2glc_coupling = .true.
          exit
       end if
    end do
    if (.not. is_local%wrap%ocn2glc_coupling) then
       ! Reset ocn2glc active coupling based in input attribute
       do ns = 1,is_local%wrap%num_icesheets
          is_local%wrap%med_coupling_active(compocn,compglc(ns)) = .false.
       end do
    end if

    ! Dealloate memory
    deallocate(med_coupling_allowed)

  end subroutine med_internalstate_coupling

  subroutine med_internalstate_defaultmasks(gcomp, rc)

    use med_constants_mod        , only : ispval_mask        => med_constants_ispval_mask

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    integer          , intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local

    !----------------------------------------------------------
    ! Default masking: for each component, the first element is
    ! when it is the src and the second element is when it is
    ! the destination
    !----------------------------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(defaultMasks(ncomps,2))
    defaultMasks(:,:) = ispval_mask
    if (is_local%wrap%comp_present(compocn)) defaultMasks(compocn,:) = 0
    if (is_local%wrap%comp_present(compice)) defaultMasks(compice,:) = 0
    if (is_local%wrap%comp_present(compwav)) defaultMasks(compwav,:) = 0
    if ( coupling_mode(1:3) == 'ufs') then
       if (is_local%wrap%comp_present(compatm)) defaultMasks(compatm,:) = 1
    endif
    if ( trim(coupling_mode) == 'hafs') then
       if (is_local%wrap%comp_present(compatm)) defaultMasks(compatm,1) = 1
    endif
    if ( trim(coupling_mode) /= 'cesm') then
       if (is_local%wrap%comp_present(compatm) .and. trim(atm_name(1:4)) == 'datm') then
          defaultMasks(compatm,1) = 0
       end if
    end if

  end subroutine med_internalstate_defaultmasks

end module med_internalstate_mod
