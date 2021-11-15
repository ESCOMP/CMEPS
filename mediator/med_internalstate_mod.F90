module med_internalstate_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use ESMF         , only : ESMF_RouteHandle, ESMF_FieldBundle, ESMF_State, ESMF_Field, ESMF_VM
  use ESMF         , only : ESMF_GridComp, ESMF_MAXSTR
  use med_kind_mod , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  public :: med_internalstate_init

  integer, public :: logunit            ! logunit for mediator log output
  integer, public :: diagunit           ! diagunit for budget output (med master only)
  integer, public :: loglevel           ! loglevel for mediator log output
  logical, public :: mastertask=.false. ! is this the mastertask
  integer, public :: med_id             ! needed currently in med_io_mod and set in esm.F90

  integer, public :: ncomps  = 0
  integer, public :: compmed = 0
  integer, public :: compatm = 0
  integer, public :: complnd = 0
  integer, public :: compocn = 0
  integer, public :: compice = 0
  integer, public :: comprof = 0
  integer, public :: compwav = 0
  integer, public, allocatable :: compglc(:)

  character(len=CS), public, allocatable :: compname(:)
  character(len=CS), public :: med_name = ''
  character(len=CS), public :: atm_name = ''
  character(len=CS), public :: lnd_name = ''
  character(len=CS), public :: ocn_name = ''
  character(len=CS), public :: ice_name = ''
  character(len=CS), public :: rof_name = ''
  character(len=CS), public :: wav_name = ''
  character(len=CS), public :: glc_name = ''

  integer, public :: num_icesheets     ! obtained from attribute
  logical, public :: ocn2glc_coupling  ! obtained from attribute
  logical, public :: lnd2glc_coupling  ! obtained in med.F90
  logical, public :: accum_lnd2glc     ! obtained in med.F90 (this can be true even if lnd2glc_coupling is false)
  logical, public :: dststatus_print = .false.

  ! Coupling mode
  character(len=CS), public :: coupling_mode ! valid values are [cesm,nems_orig,nems_frac,nems_orig_data,hafs]

  ! Active coupling definitions (will be initialize in med.F90)
  logical, public, allocatable :: med_coupling_allowed(:, :)

  ! Set mappers
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

  type, public ::  mesh_info_type
     real(r8), pointer :: areas(:) => null()
     real(r8), pointer :: lats(:) => null()
     real(r8), pointer :: lons(:) => null()
  end type mesh_info_type

  type, public :: packed_data_type
     integer, allocatable :: fldindex(:) ! size of number of packed fields
     character(len=CS)    :: mapnorm     ! normalization for packed field
     type(ESMF_Field)     :: field_src    ! packed sourced field
     type(ESMF_Field)     :: field_dst    ! packed destination field
     type(ESMF_Field)     :: field_fracsrc
     type(ESMF_Field)     :: field_fracdst
  end type packed_data_type

  ! private internal state to keep instance data
  type InternalStateStruct

    ! NState_Imp and NState_Exp are the standard NUOPC coupling datatypes
    ! FBImp and FBExp are the internal mediator datatypes
    ! NState_Exp(n) = FBExp(n), copied in the connector prep phase
    ! FBImp(n,n) = NState_Imp(n), copied in connector post phase
    ! FBImp(n,k) is the FBImp(n,n) interpolated to grid k
    ! RH(n,k,m) is a RH from grid n to grid k, map type m

    ! Present/Active logical flags
    logical, pointer   :: comp_present(:)          ! comp present flag
    logical, pointer   :: med_coupling_active(:,:) ! computes the active coupling

    ! Mediator vm
    type(ESMF_VM)          :: vm

    ! Global nx,ny dimensions of input arrays (needed for mediator history output)
    integer, pointer   :: nx(:), ny(:)

    ! Import/Export Scalars
    character(len=CL)      :: flds_scalar_name = ''
    integer                :: flds_scalar_num = 0
    integer                :: flds_scalar_index_nx = 0
    integer                :: flds_scalar_index_ny = 0
    integer                :: flds_scalar_index_nextsw_cday = 0
    integer                :: flds_scalar_index_precip_factor = 0
    real(r8)               :: flds_scalar_precip_factor = 1._r8  ! actual value of precip factor from ocn

    ! Import/export States and field bundles (the field bundles have the scalar fields removed)
    type(ESMF_State)       , pointer :: NStateImp(:) ! Import data from various component, on their grid
    type(ESMF_State)       , pointer :: NStateExp(:) ! Export data to various component, on their grid
    type(ESMF_FieldBundle) , pointer :: FBImp(:,:)   ! Import data from various components interpolated to various grids
    type(ESMF_FieldBundle) , pointer :: FBExp(:)     ! Export data for various components, on their grid

    ! Mediator field bundles for ocean albedo
    type(ESMF_FieldBundle) :: FBMed_ocnalb_o            ! Ocn albedo on ocn grid
    type(ESMF_FieldBundle) :: FBMed_ocnalb_a            ! Ocn albedo on atm grid
    type(packed_data_type), pointer :: packed_data_ocnalb_o2a(:) ! packed data for mapping ocn->atm

    ! Mediator field bundles and other info for atm/ocn flux computation
    character(len=CS)      :: aoflux_grid                        ! 'ogrid', 'agrid' or 'xgrid'
    type(ESMF_FieldBundle) :: FBMed_aoflux_a                     ! Ocn/Atm flux output fields on atm grid
    type(ESMF_FieldBundle) :: FBMed_aoflux_o                     ! Ocn/Atm flux output fields on ocn grid
    type(packed_data_type), pointer :: packed_data_aoflux_o2a(:) ! packed data for mapping ocn->atm

    ! Mapping
    type(ESMF_RouteHandle) , pointer :: RH(:,:,:)            ! Routehandles for pairs of components and different mappers
    type(ESMF_Field)       , pointer :: field_NormOne(:,:,:) ! Unity static normalization
    type(packed_data_type) , pointer :: packed_data(:,:,:)   ! Packed data structure needed to efficiently map field bundles

    ! Fractions
    type(ESMF_FieldBundle), pointer :: FBfrac(:)     ! Fraction data for various components, on their grid

    ! Accumulators for export field bundles
    type(ESMF_FieldBundle) :: FBExpAccumOcn      ! Accumulator for various components export on their grid
    integer                :: ExpAccumOcnCnt = 0 ! Accumulator counter for each FBExpAccum

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

  subroutine med_internalstate_init(gcomp, &
       med_present, atm_present, lnd_present, &
       ocn_present, ice_present, rof_present, wav_present, glc_present, rc)

    use ESMF         , only : ESMF_LogFoundAllocError
    use NUOPC_Comp   , only : NUOPC_CompAttributeAdd
    use NUOPC_Comp   , only : NUOPC_CompAttributeGet 
    use NUOPC_Comp   , only : NUOPC_CompAttributeSet
    use med_utils_mod, only : chkerr => med_utils_ChkErr

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    character(len=*), intent(out) :: med_present
    character(len=*), intent(out) :: atm_present
    character(len=*), intent(out) :: lnd_present
    character(len=*), intent(out) :: ocn_present
    character(len=*), intent(out) :: ice_present
    character(len=*), intent(out) :: rof_present
    character(len=*), intent(out) :: wav_present
    character(len=*), intent(out) :: glc_present
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    integer                    :: n, ns, ncomp 
    integer                    :: stat
    character(len=8)           :: cnum
    character(len=cs)          :: cvalue
    logical                    :: ispresent, isset
    character(len=ESMF_MAXSTR) :: mesh_glc
    character(len=CS)          :: attrList(8)
    character(len=*),parameter :: subname=' (med_internal_state_init)'
    !-----------------------------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeAdd(gcomp, &
         attrList=(/'atm_present','lnd_present','ocn_present','ice_present',&
                    'rof_present','wav_present','glc_present','med_present'/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    med_present = "false"
    atm_present = "false"
    lnd_present = "false"
    ocn_present = "false"
    ice_present = "false"
    rof_present = "false"
    wav_present = "false"
    glc_present = "false"

    ncomp = 0

    call NUOPC_CompAttributeGet(gcomp, name='mediator_present', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       med_present = trim(cvalue)
       if (med_present == 'true') then
          ncomp = ncomp + 1
          compmed = ncomp
       end if
    end if

    ! Note that the present flag is set to true if the component is not stub
    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'satm') then
          atm_present = "true"
       end if
       ncomp = ncomp + 1
       compatm = ncomp
       atm_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='LND_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'slnd') then
          lnd_present = "true"
       end if
       ncomp = ncomp + 1
       complnd = ncomp
       lnd_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='OCN_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'socn') then
          ocn_present = "true"
       end if
       ncomp = ncomp + 1
       compocn = ncomp
       ocn_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ICE_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sice') then
          ice_present = "true"
       end if
       ncomp = ncomp + 1
       compice = ncomp
       ice_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ROF_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'srof') then
          rof_present = "true"
       end if
       ncomp = ncomp + 1
       comprof = ncomp
       rof_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='WAV_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'swav') then
          wav_present = "true"
       end if
       ncomp = ncomp + 1
       compwav = ncomp
       wav_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='GLC_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sglc') then
          glc_present = "true"
          call NUOPC_CompAttributeGet(gcomp, name='mesh_glc', value=mesh_glc, isPresent=isPresent, isSet=isSet, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          num_icesheets = 0
          if (isPresent .and. isSet) then
             ! determine number of ice sheets - search in mesh_glc for colon deliminted strings
             if (len_trim(cvalue) > 0) then
                do n = 1, len_trim(mesh_glc)
                   if (mesh_glc(n:n) == ':') num_icesheets = num_icesheets + 1
                end do
                num_icesheets = num_icesheets + 1
             endif
             if (mastertask) then
                write(logunit,'(a,i8)') trim(subname)//' number of ice sheets is ',num_icesheets
             end if
          end if
          allocate(compglc(num_icesheets))
          do ns = 1,num_icesheets
             write(cnum,'(i0)') ns
             ncomp = ncomp + 1
             compglc(ns) = ncomp
          end do
       end if
       glc_name = trim(cvalue)
    end if

    ncomps = ncomp

    allocate(compname(ncomps))
    if (med_present == 'true') then
       compname(compmed) = 'med'
    end if
    if (atm_present == 'true') then
       compname(compatm) = 'atm'
    end if
    if (lnd_present == 'true') then
       compname(complnd) = 'lnd'
    end if
    if (ocn_present == 'true') then
       compname(compocn) = 'ocn'
    end if
    if (ice_present == 'true') then
       compname(compice) = 'ice'
    end if
    if (rof_present == 'true') then
       compname(comprof) = 'rof'
    end if
    if (wav_present == 'true') then
       compname(compwav) = 'wav'
    end if
    if (glc_present == 'true') then
       do ns = 1,num_icesheets
          write(cnum,'(i0)') ns
          compname(compglc(ns)) = 'glc' // trim(cnum)
       end do
    end if
    call NUOPC_CompAttributeSet(gcomp, name="atm_present", value=atm_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="lnd_present", value=lnd_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="ocn_present", value=ocn_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="ice_present", value=ice_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="rof_present", value=rof_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="wav_present", value=trim(wav_present), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="glc_present", value=trim(glc_present), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="med_present", value=med_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       write(logunit,*)
       if (trim(atm_present).eq."true") write(logunit,*) "atm_name="//trim(atm_name)
       if (trim(lnd_present).eq."true") write(logunit,*) "lnd_name="//trim(lnd_name)
       if (trim(ocn_present).eq."true") write(logunit,*) "ocn_name="//trim(ocn_name)
       if (trim(ice_present).eq."true") write(logunit,*) "ice_name="//trim(ice_name)
       if (trim(rof_present).eq."true") write(logunit,*) "rof_name="//trim(rof_name)
       if (trim(wav_present).eq."true") write(logunit,*) "wav_name="//trim(wav_name)
       if (trim(glc_present).eq."true") write(logunit,*) "glc_name="//trim(glc_name)
       if (trim(med_present).eq."true") write(logunit,*) "med_name="//trim(med_name)
       write(logunit,*)
    end if

    write(6,*)'DEBUG: ncomps = ',ncomps
    allocate(med_coupling_allowed(ncomps,ncomps))
    allocate(is_local%wrap%comp_present(ncomps))
    allocate(is_local%wrap%med_coupling_active(ncomps,ncomps))
    allocate(is_local%wrap%nx(ncomps))
    allocate(is_local%wrap%ny(ncomps))
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
    allocate(is_local%wrap%mesh_info(ncomps))
    allocate(is_local%wrap%FBArea(ncomps))

  end subroutine med_internalstate_init

end module med_internalstate_mod
