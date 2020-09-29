module med_fraction_mod

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  ! Sets fractions on all component grids
  !  the fractions fields are now ifrac, ofrac, lfrac
  !    lfrac = fraction of lnd on a grid
  !    ifrac = fraction of ice on a grid
  !    ofrac = fraction of ocn on a grid
  !    ifrad = fraction of ocn on a grid at last radiation time
  !    ofrad = fraction of ice on a grid at last radiation time
  !
  !   lfrac, ifrac, and ofrac:
  !       are the self-consistent values in the system
  !    ifrad and ofrad:
  !       are needed for the swnet calculation.
  !
  !  the fractions fields are defined for each grid in the fraction bundles as
  !    needed as follows.
  !    character(*),parameter :: fraclist_a = 'ifrac:ofrac:lfrac
  !    character(*),parameter :: fraclist_o = 'ifrac:ofrac:ifrad:ofrad'
  !    character(*),parameter :: fraclist_i = 'ifrac:ofrac'
  !    character(*),parameter :: fraclist_l = 'lfrac'
  !    character(*),parameter :: fraclist_g = 'gfrac:lfrac'
  !    character(*),parameter :: fraclist_r = 'lfrac:rfrac'
  !
  !  we assume ocean and ice are on the same grids, same masks
  !  we assume ocn2atm and ice2atm are masked maps
  !  we assume lnd2atm is a global map
  !  we assume that the ice fraction evolves in time but that
  !    the land model fraction does not.  the ocean fraction then
  !    is just the complement of the ice fraction over the region
  !    of the ocean/ice mask.
  !  we assume that component fractions sent at runtime
  !    are always the relative fraction covered.
  !    for example, if an ice cell can be up to 50% covered in
  !    ice and 50% land, then the ice should have a fraction
  !    value of 0.5 at that grid cell.  at run time though, the ice
  !    fraction will be between 0.0 and 1.0 meaning that grid cells
  !    is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
  !    sent at run-time are corrected by the model to be total fractions
  !    such that in general, on every grid,
  !       fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
  !  where fractions_* are a bundle of fractions on a particular grid and
  !    *frac is the fraction of a particular component in the bundle.
  !
  !  the fractions are computed fundamentally as follows (although the
  !    detailed implementation might be slightly different)
  !
  !  initialization:
  !    initially assume ifrac on all grids is zero
  !      fractions_*(ifrac) = 0.0
  !    fractions/masks provided by surface components
  !      fractions_o(ofrac) = ocean "mask" provided by ocean
  !    then mapped to the atm model
  !      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
  !    and a few things are then derived
  !      fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
  !           this is truncated to zero for very small values (< 0.001)
  !           to attempt to preserve non-land gridcells.
  !      fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
  !      fractions_r(lfrac) = mapl2r(fractions_l(lfrac))
  !      fractions_g(lfrac) = mapl2g(fractions_l(lfrac))
  !
  !  run-time (frac_set):
  !    update fractions on ice grid
  !      fractions_i(ifrac) = i2x_i(Si_ifrac)  ! ice frac from ice model
  !      fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
  !        note: the relative fractions are corrected to total fractions
  !      fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
  !      fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
  !      fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
  !      fractions_a(ofrac) = mapi2a(fractions_i(ofrac))
  !
  !  fractions used in merging are as follows
  !  merge to atm   uses fractions_a(lfrac,ofrac,ifrac)
  !  merge to ocean uses fractions_o(ofrac,ifrac) normalized to one
  !
  !  fraction corrections in mapping are as follows
  !    mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
  !    mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
  !    mapl2a uses *fractions_l(lfrac)
  !    mapl2g weights by fractions_l(lfrac) with normalization and multiplies by fractions_g(lfrac)
  !
  !  run time:
  !      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
  !      0.0-eps < fractions_*(*) < 1.0+eps
  !
  ! Note that the following FBImp field names are current hard-wired below
  ! TODO: this needs to be generalized - these names should be set dynamically at run time in the
  ! source component
  !    is_local%wrap%FBImp(compglc,compglc) => 'frac'
  !    is_local%wrap%FBImp(complnd,complnd) => 'Sl_lfrin'
  !    is_local%wrap%FBImp(compice,compice) => 'Si_imask'
  !    is_local%wrap%FBImp(compocn,compocn) => 'So_omask'
  !    is_local%wrap%FBImp(compice,compice) => 'Si_ifrac' (runtime)
  !
  !-----------------------------------------------------------------------------

  use med_kind_mod      , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod , only : dbug_flag      => med_constants_dbug_flag
  use med_constants_mod , only : czero          => med_constants_czero
  use med_utils_mod     , only : chkErr         => med_utils_ChkErr
  use med_methods_mod   , only : FB_init        => med_methods_FB_init
  use med_methods_mod   , only : FB_reset       => med_methods_FB_reset
  use med_methods_mod   , only : FB_getFldPtr   => med_methods_FB_getFldPtr
  use med_methods_mod   , only : FB_diagnose    => med_methods_FB_diagnose
  use med_methods_mod   , only : FB_fldChk      => med_methods_FB_fldChk
  use med_map_mod       , only : FB_FieldRegrid => med_map_FB_Field_Regrid
  use esmFlds           , only : ncomps
  use ESMF              , only : ESMF_RouteHandle

  implicit none
  private

  ! Note - everything is private in this module other than these routines
  public med_fraction_init
  public med_fraction_set

  integer, parameter                      :: nfracs = 5
  character(len=5)                        :: fraclist(nfracs,ncomps)
  character(len=5),parameter,dimension(3) :: fraclist_a = (/'ifrac','ofrac','lfrac'/)
  character(len=5),parameter,dimension(4) :: fraclist_o = (/'ifrac','ofrac','ifrad','ofrad'/)
  character(len=5),parameter,dimension(2) :: fraclist_i = (/'ifrac','ofrac'/)
  character(len=5),parameter,dimension(1) :: fraclist_l = (/'lfrac'/)
  character(len=5),parameter,dimension(2) :: fraclist_g = (/'gfrac','lfrac'/)
  character(len=5),parameter,dimension(2) :: fraclist_r = (/'rfrac','lfrac'/)
  character(len=5),parameter,dimension(1) :: fraclist_w = (/'wfrac'/)

  !--- standard ---
  real(R8)    , parameter :: eps_fraclim = 1.0e-03      ! truncation limit in fractions_a(lfrac)
  character(*), parameter :: u_FILE_u =  &
       __FILE__

  type(ESMF_RouteHandle) :: rh_ice2atm

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_fraction_init(gcomp, rc)

    ! Initialize FBFrac(:) field bundles

    use ESMF                  , only : ESMF_GridComp, ESMF_Field
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_GridCompGet, ESMF_StateIsCreated
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleDestroy
    use esmFlds               , only : coupling_mode
    use esmFlds               , only : compatm, compocn, compice, complnd
    use esmFlds               , only : comprof, compglc, compwav, compname
    use esmFlds               , only : mapfcopy, mapconsd, mapnstod_consd
    use med_map_mod           , only : med_map_routehandles_init, med_map_rh_is_created
    use med_internalstate_mod , only : InternalState, logunit, mastertask
    use perf_mod              , only : t_startf, t_stopf

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_MAXSTR
    use ESMF                  , only : ESMF_Field, ESMF_FieldRegrid
    use ESMF                  , only : ESMF_FieldRegridStore, ESMF_FieldRegridRelease
    use ESMF                  , only : ESMF_FieldRedistStore, ESMF_FieldRedistRelease
    use ESMF                  , only : ESMF_TERMORDER_SRCSEQ, ESMF_Region_Flag, ESMF_REGION_TOTAL, ESMF_REGION_SELECT
    use ESMF                  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_FRACAREA
    use ESMF                  , only : ESMF_NORMTYPE_DSTAREA, ESMF_REGRIDMETHOD_PATCH, ESMF_RouteHandlePrint
    use med_methods_mod       , only : FB_GetFieldByName => med_methods_FB_GetFieldByName

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    real(R8), pointer          :: frac(:)
    real(R8), pointer          :: ofrac(:)
    real(R8), pointer          :: lfrac(:)
    real(R8), pointer          :: ifrac(:)
    real(R8), pointer          :: gfrac(:)
    real(R8), pointer          :: rfrac(:)
    real(R8), pointer          :: wfrac(:)
    real(R8), pointer          :: Sl_lfrin(:)
    real(R8), pointer          :: Si_imask(:)
    real(R8), pointer          :: So_omask(:)
    integer                    :: i,j,n,n1
    integer                    :: maptype
    type(ESMF_RouteHandle)     :: rh_ocn2atm
    type(ESMF_RouteHandle)     :: rh_lnd2atm
    type(ESMF_RouteHandle)     :: rh_atm2lnd
    integer                    :: srcTermProcessing_Value = 0
    logical, save              :: first_call = .true.
    character(len=*),parameter :: subname='(med_fraction_init)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then

       !---------------------------------------
       ! Initialize the fraclist arrays
       !---------------------------------------

       fraclist(:,:) = ' '
       fraclist(1:size(fraclist_a),compatm) = fraclist_a
       fraclist(1:size(fraclist_o),compocn) = fraclist_o
       fraclist(1:size(fraclist_i),compice) = fraclist_i
       fraclist(1:size(fraclist_l),complnd) = fraclist_l
       fraclist(1:size(fraclist_r),comprof) = fraclist_r
       fraclist(1:size(fraclist_w),compwav) = fraclist_w
       fraclist(1:size(fraclist_g),compglc) = fraclist_g

       !---------------------------------------
       ! Create field bundles and initialize them to zero
       !---------------------------------------

       ! Note - must use import state here - since export state might not
       ! contain anything other than scalar data if the component is not prognostic
       do n1 = 1,ncomps
          if ( is_local%wrap%comp_present(n1) .and. &
               ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
             ! create FBFrac and zero out FBfrac(n1)
             call FB_init(is_local%wrap%FBfrac(n1), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(n1), fieldNameList=fraclist(:,n1), &
                  name='FBfrac'//trim(compname(n1)), rc=rc)
             call FB_reset(is_local%wrap%FBfrac(n1), value=czero, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
       first_call = .false.

    endif

    !---------------------------------------
    ! Set 'lfrac' for FBFrac(complnd) - this might be overwritten later
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then
       call FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd) , 'Sl_lfrin' , Sl_lfrin, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrac', lfrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       lfrac(:) = Sl_lfrin(:)
    end if

    !---------------------------------------
    ! Set 'ifrac' in FBFrac(compice) and FBFrac(compatm)
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then

       ! Set 'ifrac' FBFrac(compice)
       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , Si_imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ifrac(:) = Si_imask(:)

       ! Set 'ifrac' in  FBFrac(compatm) - at this point this is the ice mask mapped to the atm mesh
       ! This maps the ice mask (which is the same as the ocean mask) to the atm mesh
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%med_coupling_active(compice,compatm)) then

          if (med_map_RH_is_created(is_local%wrap%RH(compice,compatm,:),mapfcopy, rc=rc)) then
             ! If ice and atm are on the same mesh - a redist route handle has already been created
             maptype = mapfcopy
          else
             if (trim(coupling_mode) == 'nems_orig' ) then
                maptype = mapnstod_consd
             else
                maptype = mapconsd
             end if
             if (.not. med_map_RH_is_created(is_local%wrap%RH(compice,compatm,:),maptype, rc=rc)) then
                call med_map_routehandles_init( compice, compatm, &
                     FBSrc=is_local%wrap%FBImp(compice,compice), &
                     FBDst=is_local%wrap%FBImp(compice,compatm), &
                     mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compice), 'ifrac', &
               is_local%wrap%FBfrac(compatm), 'ifrac', &
               is_local%wrap%RH(compice,compatm,:), maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    !---------------------------------------
    ! Set 'ofrac' in FBFrac(compocn) and 'ofrac' in FBFrac(compatm)
    !---------------------------------------

    if (is_local%wrap%comp_present(compocn)) then

       ! Set 'ofrac' in FBFrac(compocn)
       call FB_getFldPtr(is_local%wrap%FBImp(compocn,compocn) , 'So_omask', So_omask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ofrac(:) = So_omask(:)

       ! Set 'ofrac' in FBFrac(compatm) - at this point this is the ocean mask mapped to the atm grid
       ! This is mapping the ocean mask to the atm grid - so in effect it is (1-land fraction) on the atm grid
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%med_coupling_active(compocn,compatm)) then

          if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:),mapfcopy, rc=rc)) then
             ! If ocn and atm are on the same mesh - a redist route handle has already been created
             maptype = mapfcopy
          else
             if (trim(coupling_mode) == 'nems_orig' ) then
                maptype = mapnstod_consd
             else
                maptype = mapconsd
             end if
             if (.not. med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:),maptype, rc=rc)) then
                call med_map_routehandles_init( compocn, compatm, &
                     FBSrc=is_local%wrap%FBImp(compocn,compocn), &
                     FBDst=is_local%wrap%FBImp(compocn,compatm), &
                     mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compocn), 'ofrac', &
               is_local%wrap%FBfrac(compatm), 'ofrac', &
               is_local%wrap%RH(compocn,compatm,:), maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       end if
    end if

    !---------------------------------------
    ! Reset 'lfrac' in FBFrac(complnd) by mapping FBFrac(compatm) if appropriate
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%med_coupling_active(complnd,compatm)) then

       ! Reset 'lfrac' in FBFrac(complnd) by mapping the Frac
       ! If lnd -> atm coupling is active - map 'lfrac' from FBFrac(compatm) to FBFrac(complnd)

       if (med_map_RH_is_created(is_local%wrap%RH(compatm,complnd,:),mapfcopy, rc=rc)) then
          maptype = mapfcopy
       else
          maptype = mapconsd
          if (.not. med_map_RH_is_created(is_local%wrap%RH(compatm,complnd,:),maptype, rc=rc)) then
             call med_map_routehandles_init( compatm, complnd, &
                  FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                  FBDst=is_local%wrap%FBImp(compatm,complnd), &
                  mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
       call FB_FieldRegrid(&
            is_local%wrap%FBfrac(compatm), 'lfrac', &
            is_local%wrap%FBfrac(complnd), 'lfrac', &
            is_local%wrap%RH(compatm,complnd,:),maptype, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    !---------------------------------------
    ! Set 'lfrac' in FBFrac(compatm) and correct 'ofrac' in FBFrac(compatm)
    ! ---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then

       if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compice)) then

          ! Ocean is present
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)

          if (is_local%wrap%comp_present(complnd)) then
             do n = 1,size(lfrac)
                lfrac(n) = 1.0_R8 - ofrac(n)
                if (abs(lfrac(n)) < eps_fraclim) then
                   lfrac(n) = 0.0_R8
                end if
             end do
          else
             lfrac(:) = 0.0_R8
          end if

       else if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%med_coupling_active(complnd,compatm)) then

          ! If the ocean or ice are absent, regrid 'lfrac' from FBFrac(complnd) -> FBFrac(compatm)
          if (med_map_RH_is_created(is_local%wrap%RH(complnd,compatm,:),mapfcopy, rc=rc)) then
             maptype = mapfcopy
          else
             maptype = mapconsd
             if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,compatm,:),maptype, rc=rc)) then
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,compatm))) then
                   call med_map_routehandles_init( complnd, compatm, &
                        FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                        FBDst=is_local%wrap%FBImp(complnd,compatm), &
                        mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrac', &
               is_local%wrap%FBfrac(compatm), 'lfrac', &
               is_local%wrap%RH(complnd,compatm,:),maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
          do n = 1,size(lfrac)
             ofrac(n) = 1.0_R8 - lfrac(n)
             if (abs(ofrac(n)) < eps_fraclim) then
                ofrac(n) = 0.0_R8
             end if
          end do

       end if
    end if

    !---------------------------------------
    ! Set 'rfrac' and 'lfrac' for FBFrac(comprof)
    !---------------------------------------

    if (is_local%wrap%comp_present(comprof)) then

       ! Set 'rfrac' in FBFrac(comprof)
       if ( FB_FldChk(is_local%wrap%FBfrac(comprof)       , 'rfrac', rc=rc) .and. &
            FB_FldChk(is_local%wrap%FBImp(comprof,comprof), 'frac' , rc=rc)) then
          call FB_getFldPtr(is_local%wrap%FBfrac(comprof)       , 'rfrac', rfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBImp(comprof,comprof), 'frac' , frac, rc=rc)
          rfrac(:) = frac(:)
       else
          ! Set 'rfrac' in FBfrac(comprof) to 1.
          call FB_getFldPtr(is_local%wrap%FBfrac(comprof), 'rfrac', rfrac, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          rfrac(:) = 1.0_R8
       endif

       ! Set 'lfrac' in FBFrac(comprof)
       if (is_local%wrap%comp_present(complnd)) then
          maptype = mapconsd
          if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,comprof,:),maptype, rc=rc)) then
             call med_map_routehandles_init( complnd, comprof, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,comprof), &
                  mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrac', &
               is_local%wrap%FBfrac(comprof), 'lfrac', &
               is_local%wrap%RH(complnd,comprof,:),maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    !---------------------------------------
    ! Set 'gfrac' and 'lfrac' for FBFrac(compglc)
    !---------------------------------------

    if (is_local%wrap%comp_present(compglc)) then
       ! Set 'gfrac' in FBFrac(compglc)
       if ( FB_FldChk(is_local%wrap%FBfrac(compglc)        , 'gfrac', rc=rc) .and. &
            FB_FldChk(is_local%wrap%FBImp(compglc, compglc), 'frac' , rc=rc)) then
          call FB_getFldPtr(is_local%wrap%FBfrac(compglc)       , 'gfrac', gfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), 'frac' , frac, rc=rc)
          gfrac(:) = frac(:)
       else
          ! Set 'gfrac' in FBfrac(compglc) to 1.
          call FB_getFldPtr(is_local%wrap%FBfrac(compglc), 'gfrac', gfrac, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          gfrac(:) = 1.0_R8
       endif

       ! Set 'lfrac' in FBFrac(compglc)
       if ( is_local%wrap%comp_present(complnd) .and. is_local%wrap%med_coupling_active(complnd,compglc)) then
          maptype = mapconsd
          if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,compglc,:),maptype, rc=rc)) then
             call med_map_routehandles_init( complnd, compglc, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,compglc), &
                  mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrac', &
               is_local%wrap%FBfrac(compglc), 'lfrac', &
               is_local%wrap%RH(complnd,compglc,:),maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    !---------------------------------------
    ! Set 'wfrac' for FBFrac(compwav)
    !---------------------------------------

    if (is_local%wrap%comp_present(compwav)) then
       ! Set 'wfrac' in FBfrac(compwav) to 1.
       call FB_getFldPtr(is_local%wrap%FBfrac(compwav), 'wfrac', wfrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       wfrac(:) = 1.0_R8
    endif

    !---------------------------------------
    ! Create route handles ocn<->ice if not created
    !---------------------------------------

    if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
       if (.not. med_map_RH_is_created(is_local%wrap%RH(compice,compocn,:),mapfcopy, rc=rc)) then
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compice,compocn))) then
             call FB_init(is_local%wrap%FBImp(compice,compocn), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(compocn), &
                  STflds=is_local%wrap%NStateImp(compice), &
                  name='FBImp'//trim(compname(compice))//'_'//trim(compname(compocn)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call med_map_routehandles_init(compice, compocn, &
               FBSrc=is_local%wrap%FBImp(compice,compice), &
               FBDst=is_local%wrap%FBImp(compice,compocn), &
               mapindex=mapfcopy, RouteHandle=is_local%wrap%RH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (.not. med_map_RH_is_created(is_local%wrap%RH(compocn,compice,:),mapfcopy, rc=rc)) then
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compocn,compice))) then
             call FB_init(is_local%wrap%FBImp(compocn,compice), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(compice), &
                  STflds=is_local%wrap%NStateImp(compocn), &
                  name='FBImp'//trim(compname(compocn))//'_'//trim(compname(compice)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call med_map_routehandles_init( compocn, compice, &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compocn,compice), &
               mapindex=mapfcopy, RouteHandle=is_local%wrap%RH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    !---------------------------------------
    ! Diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       do n = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call FB_diagnose(is_local%wrap%FBfrac(n), &
                  trim(subname) // trim(compname(n)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_fraction_init

  !-----------------------------------------------------------------------------

  subroutine med_fraction_set(gcomp, rc)

    ! Update time varying fractions

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_Field, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use esmFlds               , only : compatm, compocn, compice, compname
    use esmFlds               , only : mapfcopy, mapconsd, mapnstod_consd
    use esmFlds               , only : coupling_mode
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_RH_is_created
    use perf_mod              , only : t_startf, t_stopf
    use med_methods_mod       , only : FB_GetFieldByName => med_methods_FB_GetFieldByName

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    type(ESMF_Field)           :: fldsrc
    type(ESMF_Field)           :: flddst
    real(r8), pointer          :: lfrac(:)
    real(r8), pointer          :: ifrac(:)
    real(r8), pointer          :: ofrac(:)
    real(r8), pointer          :: Si_ifrac(:)
    real(r8), pointer          :: Si_imask(:)
    integer                    :: n
    integer                    :: maptype
    character(len=*),parameter :: subname='(med_fraction_set)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Update FBFrac(compice), FBFrac(compocn) and FBFrac(compatm) field bundles
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then

       ! -------------------------------------------
       ! Set FBfrac(compice)
       ! -------------------------------------------

       ! Si_imask is the ice domain mask which is constant over time
       ! Si_ifrac is the time evolving ice fraction on the ice grid

       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_ifrac', Si_ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , Si_imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! set ifrac = Si_ifrac * Si_imask
       ifrac(:) = Si_ifrac(:) * Si_imask(:)

       ! set ofrac = Si_imask - ifrac
       ofrac(:) = Si_imask(:) - ifrac(:)

       ! -------------------------------------------
       ! Set FBfrac(compocn)
       ! -------------------------------------------

       ! The following is just a redistribution from FBFrac(compice)

       if (is_local%wrap%comp_present(compocn)) then
          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compocn)
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compice), 'ifrac', &
               is_local%wrap%FBfrac(compocn), 'ifrac', &
               is_local%wrap%RH(compice,compocn,:),mapfcopy, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compocn)
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compice), 'ofrac', &
               is_local%wrap%FBfrac(compocn), 'ofrac', &
               is_local%wrap%RH(compice,compocn,:),mapfcopy, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! -------------------------------------------
       ! Set FBfrac(compatm)
       ! -------------------------------------------

       if (is_local%wrap%comp_present(compatm)) then

          ! Determine maptype
          if (trim(coupling_mode) == 'nems_orig' ) then
             maptype = mapnstod_consd
          else
             if (med_map_RH_is_created(is_local%wrap%RH(compice,compatm,:),mapfcopy, rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsd
             end if
          end if

          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compatm)
          if (is_local%wrap%med_coupling_active(compice,compatm)) then
             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compatm), 'ifrac', &
                  is_local%wrap%RH(compice,compatm,:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compatm)
          if (is_local%wrap%med_coupling_active(compocn,compatm)) then
             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compocn), 'ofrac', &
                  is_local%wrap%FBfrac(compatm), 'ofrac', &
                  is_local%wrap%RH(compocn,compatm,:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if ! end of if present compatm

    end if ! end of if present compice

    !---------------------------------------
    ! Diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       do n = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call FB_diagnose(is_local%wrap%FBfrac(n), trim(subname) // trim(compname(n))//' frac', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       enddo
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_fraction_set

end module med_fraction_mod
