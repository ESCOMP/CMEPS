module esmFldsExchange_hafs_mod

  use ESMF
  use NUOPC
  use med_utils_mod, only : chkerr => med_utils_chkerr
  use med_kind_mod,  only : CX=>SHR_KIND_CX
  use med_kind_mod,  only : CS=>SHR_KIND_CS
  use med_kind_mod,  only : CL=>SHR_KIND_CL
  use med_kind_mod,  only : R8=>SHR_KIND_R8
  use esmflds,       only : compmed
  use esmflds,       only : compatm
  use esmflds,       only : compocn
  use esmflds,       only : compice
  use esmflds,       only : ncomps
  use esmflds,       only : fldListTo
  use esmflds,       only : fldListFr
  use esmFlds,       only : coupling_mode

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_hafs

  character(*), parameter :: u_FILE_u = &
       __FILE__

  type systemType
    sequence
    private
      integer :: system
  end type

  type(systemType), parameter ::  &
    SYS_ERR = systemType(-1), & ! Error code
    SYS_CDP = systemType(0),  & ! Community Data Models for Earth Prediction Sys
    SYS_UFS = systemType(1)     ! Unified Forecast System

  type gcomp_attr
    character(len=CX)   :: atm2ice_fmap='unset'
    character(len=CX)   :: atm2ice_smap='unset'
    character(len=CX)   :: atm2ice_vmap='unset'
    character(len=CX)   :: atm2ocn_fmap='unset'
    character(len=CX)   :: atm2ocn_smap='unset'
    character(len=CX)   :: atm2ocn_vmap='unset'
    character(len=CX)   :: ice2atm_fmap='unset'
    character(len=CX)   :: ice2atm_smap='unset'
    character(len=CX)   :: ocn2atm_fmap='unset'
    character(len=CX)   :: ocn2atm_smap='unset'
    character(len=CS)   :: mapnorm     ='one'
    type(systemType)    :: hafs_sysType=SYS_CDP
  end type

  interface operator (==)
    module procedure systemType_eq
  end interface

  interface assignment (=)
    module procedure systemType_tostring
    module procedure systemType_frstring
  end interface

!===============================================================================
contains
!===============================================================================

  subroutine esmFldsExchange_hafs(gcomp, phase, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    if (phase == 'advertise') then
      call esmFldsExchange_hafs_advt(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'fieldcheck') then
      call esmFldsExchange_hafs_fchk(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'initialize') then
      call esmFldsExchange_hafs_init(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogSetError(ESMF_FAILURE, &
         msg=trim(subname)//": Phase is set to "//trim(phase), &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_advt(gcomp, phase, rc)

    use esmFlds               , only : addfld => med_fldList_AddFld

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    integer             :: num, i, n
    logical             :: isPresent
    !character(len=5)    :: iso(2)
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: hafs_attr
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_advt)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !=====================================================================
    ! scalar information
    !=====================================================================

    call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", &
          value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld(fldListFr(n)%flds, trim(cvalue))
          call addfld(fldListTo(n)%flds, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------------------------------------------------------
    ! CDEPS Coupling
    ! ---------------------------------------------------------------------
    if (hafs_attr%hafs_sysType == SYS_CDP) then

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    call addfld(fldListFr(compocn)%flds, 'So_omask')
    call addfld(fldListFr(compice)%flds, 'Si_imask')

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    !----------------------------------------------------------
    ! to atm: Fractions
    !----------------------------------------------------------
    ! the following are computed in med_phases_prep_atm
    call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
    call addfld(fldListTo(compatm)%flds, 'So_ofrac')

    !----------------------------------------------------------
    ! to atm: surface temperatures from ocn
    !----------------------------------------------------------
    call addfld(fldListFr(compocn)%flds, 'So_t')
    call addfld(fldListTo(compatm)%flds, 'So_t')

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: req. fields to satisfy mediator (can be removed later)
    !----------------------------------------------------------
    allocate(flds(2))
    flds = (/'Faxa_snowc',&
             'Faxa_snowl'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
    end do
    deallocate(flds)

    allocate(flds(4))
    flds = (/'Sa_topo',&
             'Sa_z   ',&
             'Sa_ptem',&
             'Sa_pbot'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
    end do
    deallocate(flds)

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    call addfld(fldListFr(compice)%flds, 'Si_ifrac')
    call addfld(fldListTo(compocn)%flds, 'Si_ifrac')

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_lwdn ', &
             'Faxa_swndr', &
             'Faxa_swndf', &
             'Faxa_swvdr', &
             'Faxa_swvdf'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: longwave net heat flux
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds, 'Faxa_lwnet')
    call addfld(fldListTo(compocn)%flds, 'Foxx_lwnet')

    ! ---------------------------------------------------------------------
    ! to ocn: downward shortwave heat flux
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds, 'Faxa_swdn')
    call addfld(fldListTo(compocn)%flds, 'Faxa_swdn')

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from atm
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    call addfld(fldListTo(compocn)%flds, 'Foxx_swnet')

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate from atm
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
    call addfld(fldListFr(compatm)%flds, 'Faxa_rainl')
    call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
    call addfld(fldListTo(compocn)%flds, 'Faxa_rain' )

    ! ---------------------------------------------------------------------
    ! to ocn: sensible heat flux from atm
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds , 'Faxa_sen')
    call addfld(fldListTo(compocn)%flds , 'Foxx_sen')

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    call addfld(fldListFr(compatm)%flds , 'Faxa_lat')
    call addfld(fldListTo(compocn)%flds , 'Foxx_lat')

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! to ocn: zonal wind at the lowest model level from atm
    ! to ocn: meridional wind at the lowest model level from atm
    ! to ocn: wind speed at the lowest model level from atm
    ! to ocn: temperature at the lowest model level from atm
    ! to ocn: sea surface skin temperature
    ! to ocn: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(7))
    flds = (/'Sa_pslv', &
             'Sa_u   ', &
             'Sa_v   ', &
             'Sa_wspd', &
             'Sa_tbot', &
             'Sa_tskn', &
             'Sa_shum'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: zonal and meridional surface stress from atm
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       call addfld(fldListFr(compatm)%flds , 'Faxa_'//trim(suffix(n)))
       call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(suffix(n)))
    end do
    deallocate(suffix)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: density at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(1))
    flds = (/'Sa_dens'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compice)%flds, trim(fldname))
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: zonal sea water velocity from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(2))
    flds = (/'So_u   ', 'So_v   '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compocn)%flds, trim(fldname))
       call addfld(fldListTo(compice)%flds, trim(fldname))
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! UFS Coupling
    ! ---------------------------------------------------------------------
    elseif (hafs_attr%hafs_sysType == SYS_UFS) then
       !----------------------------------------------------------
       ! to med: masks from components
       !----------------------------------------------------------
       call addfld(fldListFr(compocn)%flds, 'So_omask')
       !----------------------------------------------------------
       ! to med: frac from components
       !----------------------------------------------------------
       call addfld(fldListTo(compatm)%flds, 'So_ofrac')
       ! ---------------------------------------------------------------------
       ! from atm to ocn
       ! ---------------------------------------------------------------------
       ! state fields
       allocate(S_flds(6))
       S_flds = (/'Sa_u10m', & ! inst_zonal_wind_height10m
                  'Sa_v10m', & ! inst_merid_wind_height10m
                  'Sa_t2m ', & ! inst_temp_height2m
                  'Sa_q2m ', & ! inst_spec_humid_height2m
                  'Sa_pslv', & ! inst_pres_height_surface
                  'Sa_tskn' /) ! inst_temp_height_surface
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       end do
       deallocate(S_flds)
       ! flux fields
       allocate(F_flds(5,2))
       F_flds(1,:) = (/'Faxa_taux ','Faxa_taux '/) ! mean_zonal_moment_flx_atm
       F_flds(2,:) = (/'Faxa_tauy ','Faxa_tauy '/) ! mean_merid_moment_flx_atm
       F_flds(3,:) = (/'Faxa_rain ','Faxa_rain '/) ! mean_prec_rate
       F_flds(4,:) = (/'Faxa_swnet','Faxa_swnet'/) ! mean_net_sw_flx
       F_flds(5,:) = (/'Faxa_lwnet','Faxa_lwnet'/) ! mean_net_lw_flx
       do n = 1,size(F_flds,1)
          fldname1 = trim(F_flds(n,1))
          fldname2 = trim(F_flds(n,2))
          call addfld(fldListFr(compatm)%flds, trim(fldname1))
          call addfld(fldListTo(compocn)%flds, trim(fldname2))
       end do
       deallocate(F_flds)
       ! ---------------------------------------------------------------------
       ! from ocn to atm
       ! ---------------------------------------------------------------------
       ! state fields
       allocate(S_flds(1))
       S_flds = (/'So_t'/) ! sea_surface_temperature
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          call addfld(fldListFr(compocn)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       end do
       deallocate(S_flds)
    endif ! hafs_sysType

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_advt

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_fchk(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_fchk)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fldchk(is_local%wrap%FBImp(compocn,compocn),'So_omask',rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": Field connected "//"So_omask", &
          ESMF_LOGMSG_INFO)
    else
       call ESMF_LogSetError(ESMF_FAILURE, &
          msg=trim(subname)//": Field is not connected "//"So_omask", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_fchk

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_init(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd
    use esmflds               , only : mapfillv_bilnr
    use esmflds               , only : mapnstod_consf

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: num, i, n
    integer             :: n1, n2, n3, n4
    !character(len=5)    :: iso(2)
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: hafs_attr
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_init)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------
    ! Merging arguments:
    ! mrg_fromN = source component index that for the field to be merged
    ! mrg_fldN  = souce field name to be merged
    ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum',
    !                         'sum_with_weights', 'merge')
    ! NOTE:
    ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn
    ! fluxes or for ocn albedos
    !
    ! NOTE:
    ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in
    ! the atm/ocn flux calculations. Input fields required from either the atm
    ! or the ocn for these computation will use the logical 'use_med_aoflux'
    ! below. This is used to determine mappings between the atm and ocn needed
    ! for these computations.
    !--------------------------------------

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------------------------------------------------------
    ! CDEPS Coupling
    ! ---------------------------------------------------------------------
    if (hafs_attr%hafs_sysType == SYS_CDP) then

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    call addmap(fldListFr(compocn)%flds, 'So_omask', compice, &
         mapfcopy, 'unset', 'unset')

    ! ---------------------------------------------------------------------
    ! to med: atm and ocn fields required for atm/ocn flux calculation
    ! ---------------------------------------------------------------------
    allocate(flds(6))
    flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', 'Sa_shum'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (trim(fldname) == 'Sa_u' .or. trim(fldname) == 'Sa_v') then
          !call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
          !     mappatch, 'one', hafs_attr%atm2ocn_vmap)
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapbilnr, 'one', hafs_attr%atm2ocn_smap)
       else
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapbilnr, 'one', hafs_attr%atm2ocn_smap)
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to med: unused fields needed by the atm/ocn flux computation
    ! ---------------------------------------------------------------------
    !call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, &
    !     mappatch, 'one', hafs_attr%atm2ocn_vmap)
    call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_vmap)
    !call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, &
    !     mappatch, 'one', hafs_attr%atm2ocn_vmap)
    call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_vmap)
    call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    call addmap(fldListFr(compatm)%flds, 'Sa_shum', compocn, &
         mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    if (fldchk(is_local%wrap%FBImp(compatm,compatm),'Sa_ptem',rc=rc)) then
       call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, &
            mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    end if
    if (fldchk(is_local%wrap%FBImp(compatm,compatm),'Sa_dens',rc=rc)) then
       call addmap(fldListFr(compatm)%flds, 'Sa_dens', compocn, &
            mapbilnr, 'one', hafs_attr%atm2ocn_smap)
    end if

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, &
         mapconsf, 'one', hafs_attr%atm2ocn_fmap)

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    !----------------------------------------------------------
    ! to atm: sea surface temperature
    !----------------------------------------------------------
    allocate(flds(1))
    flds = (/'So_t'/) ! sea_surface_temperature
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBExp(compatm),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compocn)%flds, trim(fldname), compatm, &
               mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%ocn2atm_smap)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: req. fields to satisfy mediator (can be removed later)
    !----------------------------------------------------------
    allocate(flds(2))
    flds = (/'Faxa_snowc', 'Faxa_snowl'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBexp(compocn),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapconsf, 'one', hafs_attr%atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    allocate(flds(4))
    flds = (/'Sa_topo', 'Sa_z   ', 'Sa_ptem', 'Sa_pbot'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBexp(compocn),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapbilnr, 'one', hafs_attr%atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn, &
         mapfcopy, 'unset', 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', &
         mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_lwdn ', 'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', &
             'Faxa_swvdf'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapconsf, 'one', hafs_attr%atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), &
               mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: longwave net heat flux
    ! ---------------------------------------------------------------------
    call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, &
         mapconsf, 'one', hafs_attr%atm2ocn_fmap)
    call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
         mrg_from=compatm, mrg_fld='Faxa_lwnet', mrg_type='copy')

    ! ---------------------------------------------------------------------
    ! to ocn: downward shortwave heat flux
    ! ---------------------------------------------------------------------
    if (fldchk(is_local%wrap%FBImp(compatm,compatm),'Faxa_swdn',rc=rc) .and. &
        fldchk(is_local%wrap%FBExp(compocn),'Faxa_swdn',rc=rc) &
       ) then
       call addmap(fldListFr(compatm)%flds, 'Faxa_swdn', compocn, &
            mapconsf, 'one', hafs_attr%atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Faxa_swdn', &
            mrg_from=compatm, mrg_fld='Faxa_swdn', mrg_type='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from atm
    ! ---------------------------------------------------------------------
    call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, &
         mapconsf, 'one', hafs_attr%atm2ocn_fmap)
    call addmrg(fldListTo(compocn)%flds, 'Foxx_swnet', &
         mrg_from=compatm, mrg_fld='Faxa_swnet', mrg_type='copy')

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate from atm
    ! ---------------------------------------------------------------------
    if (fldchk(is_local%wrap%FBImp(compatm,compatm),'Faxa_rainl',rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(compatm,compatm),'Faxa_rainc',rc=rc) .and. &
        fldchk(is_local%wrap%FBExp(compocn),'Faxa_rain',rc=rc) &
       ) then
        call addmap(fldListFr(compatm)%flds, 'Faxa_rainl', compocn, &
             mapconsf, 'one', hafs_attr%atm2ocn_fmap)
        call addmap(fldListFr(compatm)%flds, 'Faxa_rainc', compocn, &
             mapconsf, 'one', hafs_attr%atm2ocn_fmap)
        call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', &
             mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
             mrg_type='sum_with_weights', mrg_fracname='ofrac')
    else if (fldchk(is_local%wrap%FBExp(compocn),'Faxa_rain',rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm,compatm),'Faxa_rain',rc=rc) &
            ) then
        call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compocn, &
             mapconsf, 'one', hafs_attr%atm2ocn_fmap)
        call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', &
             mrg_from=compatm, mrg_fld='Faxa_rain', mrg_type='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sensible heat flux from atm
    ! ---------------------------------------------------------------------
    call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, &
         mapconsf, 'one', hafs_attr%atm2ocn_fmap)
    call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
         mrg_from=compatm, mrg_fld='Faxa_sen', mrg_type='copy')

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, &
         mapconsf, 'one', hafs_attr%atm2ocn_fmap)
    call addmrg(fldListTo(compocn)%flds, 'Foxx_lat', &
         mrg_from=compatm, mrg_fld='Faxa_lat', mrg_type='copy')

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! to ocn: zonal wind at the lowest model level from atm
    ! to ocn: meridional wind at the lowest model level from atm
    ! to ocn: wind speed at the lowest model level from atm
    ! to ocn: temperature at the lowest model level from atm
    ! to ocn: sea surface skin temperature
    ! to ocn: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(7))
    flds = (/'Sa_pslv', 'Sa_u   ', 'Sa_v   ', 'Sa_wspd', 'Sa_tbot', 'Sa_tskn', &
             'Sa_shum'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapbilnr, 'one', hafs_attr%atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: zonal and meridional surface stress from atm
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)), compocn, &
            mapconsf, 'one', hafs_attr%atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)), &
            mrg_from=compatm, mrg_fld='Faxa_'//trim(suffix(n)), &
            mrg_type='copy')
    end do
    deallocate(suffix)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: density at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(1))
    flds = (/'Sa_dens'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBexp(compice),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compice, &
               mapbilnr, 'one', hafs_attr%atm2ice_smap)
          call addmrg(fldListTo(compice)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: zonal sea water velocity from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(2))
    flds = (/'So_u   ', 'So_v   '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (fldchk(is_local%wrap%FBexp(compice),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compocn)%flds, trim(fldname), compice, &
               mapfcopy , 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, trim(fldname), &
               mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! UFS Coupling
    ! ---------------------------------------------------------------------
    elseif (hafs_attr%hafs_sysType == SYS_UFS) then
       ! ---------------------------------------------------------------------
       ! from atm to ocn
       ! ---------------------------------------------------------------------
       ! state fields
       allocate(S_flds(6))
       S_flds = (/'Sa_u   ', & ! inst_zonal_wind_height10m
                'Sa_v   ', & ! inst_merid_wind_height10m
                'Sa_tbot', & ! inst_temp_height2m
                'Sa_shum', & ! inst_spec_humid_height2m
                'Sa_pslv', & ! inst_pres_height_surface
                'Sa_tskn' /) ! inst_temp_height_surface
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
             ) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
                  mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%atm2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end do
       deallocate(S_flds)
       ! flux fields
       allocate(F_flds(5,2))
       F_flds(1,:) = (/'Faxa_taux ','Faxa_taux '/) ! mean_zonal_moment_flx_atm
       F_flds(2,:) = (/'Faxa_tauy ','Faxa_tauy '/) ! mean_merid_moment_flx_atm
       F_flds(3,:) = (/'Faxa_rain ','Faxa_rain '/) ! mean_prec_rate
       F_flds(4,:) = (/'Faxa_swnet','Faxa_swnet'/) ! mean_net_sw_flx
       F_flds(5,:) = (/'Faxa_lwnet','Faxa_lwnet'/) ! mean_net_lw_flx
       do n = 1,size(F_flds,1)
          fldname1 = trim(F_flds(n,1))
          fldname2 = trim(F_flds(n,2))
          if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname2),rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname1),rc=rc) &
            ) then
             call addmap(fldListFr(compatm)%flds, trim(fldname1), compocn, &
                  mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%atm2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname2), &
                  mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')
          end if
       end do
       deallocate(F_flds)
       ! ---------------------------------------------------------------------
       ! from ocn to atm
       ! ---------------------------------------------------------------------
       ! state fields
       allocate(S_flds(1))
       S_flds = (/'So_t'/) ! sea_surface_temperature
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          if (fldchk(is_local%wrap%FBExp(compatm),trim(fldname),rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
             ) then
             call addmap(fldListFr(compocn)%flds, trim(fldname), compatm, &
                  mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%ocn2atm_smap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end do
       deallocate(S_flds)
    endif ! hafs_sysType

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_init

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    type(gcomp_attr) , intent(inout) :: hafs_attr
    integer          , intent(inout) :: rc

    ! local variables:
    character(32)       :: cname
    integer             :: verbosity, diagnostic
    character(len=CL)   :: cvalue
    logical             :: isPresent
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_attr)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
#if ESMF_VERSION_MAJOR >= 8
    call NUOPC_CompGet(gcomp, name=cname, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
#else
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=cvalue, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    verbosity = ESMF_UtilString2Int(cvalue, &
      specialStringList=(/"off ","low ","high","max "/), &
      specialValueList=(/0,9985,32513,131071/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=cvalue, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    diagnostic = ESMF_UtilString2Int(cvalue, &
      specialStringList=(/"off ","max "/), &
      specialValueList=(/0,131071/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
#endif

    !----------------------------------------------------------
    ! Initialize system type
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='system_type', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='system_type', &
          value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       hafs_attr%hafs_sysType = cvalue
    end if

    !----------------------------------------------------------
    ! Normalization type
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='normalization', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='normalization', &
          value=hafs_attr%mapnorm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm
    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmapname', &
          value=hafs_attr%ice2atm_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smapname', &
          value=hafs_attr%ice2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
          value=hafs_attr%ocn2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
          value=hafs_attr%ocn2atm_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ice
    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmapname', &
          value=hafs_attr%atm2ice_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', &
          value=hafs_attr%atm2ice_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmapname', &
          value=hafs_attr%atm2ice_vmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ocn
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
          value=hafs_attr%atm2ocn_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       value=hafs_attr%atm2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
          value=hafs_attr%atm2ocn_vmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Log Attribute Settings
    if (btest(verbosity,16)) then
       write(cvalue,"(I0)") verbosity
       call ESMF_LogWrite(trim(subname)//': Verbosity        = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       write(cvalue,"(I0)") diagnostic
       call ESMF_LogWrite(trim(subname)//': Diagnostic       = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       cvalue = hafs_attr%hafs_sysType
       call ESMF_LogWrite(trim(subname)//': system_type      = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': normalization    = '// &
          trim(hafs_attr%mapnorm), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ice2atm_fmapname = '// &
          trim(hafs_attr%ice2atm_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ice2atm_smapname = '// &
          trim(hafs_attr%ice2atm_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_smapname = '// &
          trim(hafs_attr%ocn2atm_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_fmapname = '// &
          trim(hafs_attr%ocn2atm_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ice_fmapname = '// &
          trim(hafs_attr%atm2ice_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ice_smapname = '// &
          trim(hafs_attr%atm2ice_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ice_vmapname = '// &
          trim(hafs_attr%atm2ice_vmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_fmapname = '// &
          trim(hafs_attr%atm2ocn_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_smapname = '// &
          trim(hafs_attr%atm2ocn_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_vmapname = '// &
          trim(hafs_attr%atm2ocn_vmap), ESMF_LOGMSG_INFO)
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_attr

  !-----------------------------------------------------------------------------

  function systemType_eq(type1, type2)
    logical systemType_eq
    type(systemType), intent(in) :: type1, type2
    systemType_eq = (type1%system == type2%system)
  end function

  !-----------------------------------------------------------------------------

  subroutine systemType_tostring(string, tval)
    character(len=*), intent(out) :: string
    type(systemType), intent(in)  :: tval
    select case (tval%system)
      case(SYS_CDP%system)
        string = 'CDEPS'
      case(SYS_UFS%system)
        string = 'UFS'
      case default
        string = 'ERROR'
    end select
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine systemType_frstring(tval, string)
    type(systemType), intent(out)  :: tval
    character(len=*), intent(in)   :: string
    select case (ESMF_UtilStringUpperCase(string))
      case ('CDEPS')
        tval = SYS_CDP
      case ('UFS')
        tval = SYS_UFS
      case default
        tval = SYS_ERR
    end select
  end subroutine

end module esmFldsExchange_hafs_mod
