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
  !    ofrad = fraction of ice on a grid at last radiation time
  !      afrac, lfrac, ifrac, and ofrac are the self-consistent values in the
  !      system.  lfrin is the fraction on the land grid and is allowed to
  !      vary from the self-consistent value as descibed below.  ifrad
  !      and ofrad are needed for the swnet calculation.
  !   lfrac, ifrac, and ofrac:
  !       are the self-consistent values in the system
  !    ifrad and ofrad:
  !       are needed for the swnet calculation.
  !
  !  the fractions fields are defined for each grid in the fraction bundles as
  !    needed as follows.
  !    character(*),parameter :: fraclist_a = 'ifrac:ofrac:lfrac:lfrin:aofrac
  !    character(*),parameter :: fraclist_o = 'ifrac:ofrac:ifrad:ofrad'
  !    character(*),parameter :: fraclist_i = 'ifrac:ofrac'
  !    character(*),parameter :: fraclist_l = 'lfrac:lfrin'
  !    character(*),parameter :: fraclist_g = 'gfrac:lfrac:lfrin'
  !    character(*),parameter :: fraclist_r = 'rfrac:lfrac:lfrin'
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
  !  in general, on every grid,
  !    fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
  !
  !  the fractions are computed fundamentally as follows (although the
  !    detailed implementation might be slightly different)
  !
  !  initialization:
  !    initially assume ifrac on all grids is zero
  !      fractions_*(ifrac) = 0.0
  !    fractions/masks provided by surface components
  !      fractions_o(ofrac) = ocean "mask" provided by ocean
  !      fractions_l(lfrin) = Sl_lfrin  ! land model fraction computed as
  !                                       map of ocean mask to land grid
  !    then mapped to the atm model
  !      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
  !      fractions_a(lfrin) = mapl2a(fractions_l(lfrin))
  !
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
  !    mapl2a uses *fractions_l(lfrin) and /fractions_a(lfrin)
  !    mapl2g weights by fractions_l(lfrin) with normalization and multiplies by fractions_g(lfrin) ???
  !
  !  run time:
  !      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
  !      0.0-eps < fractions_*(*) < 1.0+eps
  !
  ! Note that the following FBImp field names are current hard-wired below
  ! TODO: this needs to be generalized - these names should be set dynamically at run time in the
  ! source component
  !    is_local%wrap%FBImp(compglc,compglc(:)) => 'frac'
  !    is_local%wrap%FBImp(complnd,complnd)    => 'Sl_lfrin'
  !    is_local%wrap%FBImp(compice,compice)    => 'Si_imask'
  !    is_local%wrap%FBImp(compocn,compocn)    => 'So_omask'
  !    is_local%wrap%FBImp(compice,compice)    => 'Si_ifrac' (runtime)
  !
  !  NOTE: In trigrid configurations, lfrin MUST be defined as the
  !  conservative o2l mapping of the complement of the ocean mask.
  !  In non-trigrid configurations, lfrin is generally associated with
  !  the fraction of land grid defined by the surface dataset and might
  !  be 1 everywhere for instance.  In many cases, the non-trigrid
  !  lfrin is defined to be the conservative o2a mapping of the complement
  !  of the ocean mask.  In this case, it is defined the same as the
  !  trigrid.  But to support all cases,
  !  for trigrid:
  !    mapping from the land grid should use the lfrin field (same in non-trigrid)
  !    budget diagnostics should use lfrin (lfrac in non-trigrid)
  !    merges in the atm should use lfrac (same in non-trigrid)
  !    the runoff should use the lfrin fraction in the runoff merge (lfrac in non-trigrid)
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX =>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_constants_mod     , only : czero            => med_constants_czero
  use med_utils_mod         , only : chkErr           => med_utils_ChkErr
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_fldchk    => med_methods_FB_fldchk
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use med_methods_mod       , only : fldbun_getdata2d => med_methods_FB_getdata2d
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_init      => med_methods_FB_init
  use med_methods_mod       , only : fldbun_reset     => med_methods_FB_reset
  use med_map_mod           , only : med_map_field
  use med_internalstate_mod , only : ncomps, samegrid_atmlnd

  implicit none
  private

  ! Note - everything is private in this module other than these routines
  public med_fraction_init
  public med_fraction_set

  integer, parameter           :: nfracs = 5
  character(len=6),allocatable :: fraclist(:,:)
  character(len=6),parameter   :: fraclist_a(5) = (/'ifrac ','ofrac ','lfrac ','lfrin ','aofrac'/)
  character(len=6),parameter   :: fraclist_o(4) = (/'ifrac ','ofrac ','ifrad ','ofrad '/)
  character(len=6),parameter   :: fraclist_i(2) = (/'ifrac ','ofrac '/)
  character(len=6),parameter   :: fraclist_l(2) = (/'lfrac ','lfrin '/)
  character(len=6),parameter   :: fraclist_g(3) = (/'gfrac ','lfrac ','lfrin '/)
  character(len=6),parameter   :: fraclist_r(3) = (/'rfrac ','lfrac ','lfrin '/)
  character(len=6),parameter   :: fraclist_w(1) = (/'wfrac '/)

  !--- standard ---
  real(R8)    , parameter :: eps_fraclim = 1.0e-03      ! truncation limit in fractions_a(lfrac)
  character(*), parameter :: u_FILE_u =  &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_fraction_init(gcomp, rc)

    ! Initialize FBFrac(:) field bundles

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_LogSetError, ESMF_RC_NOT_VALID
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_StateIsCreated
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleDestroy
    use ESMF                  , only : ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet
    use med_internalstate_mod , only : coupling_mode
    use med_internalstate_mod , only : compatm, compocn, compice, complnd
    use med_internalstate_mod , only : comprof, compglc, compwav, compname
    use med_internalstate_mod , only : mapfcopy, mapconsd, mapnstod_consd
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_routehandles_init, med_map_rh_is_created
    use med_methods_mod       , only : State_getNumFields => med_methods_State_getNumFields
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    real(R8), pointer   :: frac(:)
    real(R8), pointer   :: ofrac(:)
    real(R8), pointer   :: aofrac(:)
    real(R8), pointer   :: lfrac(:)
    real(R8), pointer   :: lfrin(:)
    real(R8), pointer   :: ifrac(:)
    real(R8), pointer   :: gfrac(:)
    real(R8), pointer   :: rfrac(:)
    real(R8), pointer   :: wfrac(:)
    real(R8), pointer   :: Sl_lfrin(:)
    real(R8), pointer   :: Si_imask(:)
    real(R8), pointer   :: So_omask(:)
    real(R8), pointer   :: Sa_ofrac(:)
    integer             :: n,n1,ns
    integer             :: maptype
    integer             :: fieldCount
    logical, save       :: first_call = .true.
    character(len=*),parameter :: subname=' (med_fraction_init)'
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

       ! allocate module variable
       allocate(fraclist(nfracs,ncomps))

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
       do ns = 1,is_local%wrap%num_icesheets
          fraclist(1:size(fraclist_g),compglc(ns)) = fraclist_g
       end do

       !---------------------------------------
       ! Create field bundles and initialize them to zero
       !---------------------------------------

       ! Note - must use import state here - since export state might not
       ! contain anything other than scalar data if the component is not prognostic
       do n1 = 1,ncomps
          if ( is_local%wrap%comp_present(n1) .and. &
              (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .or. &
               ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc))) then
             ! Check number of fields in the state
             call State_GetNumFields(is_local%wrap%NStateImp(n1), fieldCount, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! create FBFrac
             if (fieldCount == 0) then
               call fldbun_init(is_local%wrap%FBfrac(n1), is_local%wrap%flds_scalar_name, &
                    STgeom=is_local%wrap%NStateExp(n1), fieldNameList=fraclist(:,n1), &
                    name='FBfrac'//trim(compname(n1)), rc=rc)
             else
               call fldbun_init(is_local%wrap%FBfrac(n1), is_local%wrap%flds_scalar_name, &
                    STgeom=is_local%wrap%NStateImp(n1), fieldNameList=fraclist(:,n1), &
                    name='FBfrac'//trim(compname(n1)), rc=rc)
             end if
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! zero out FBfrac(n1)
             call fldbun_reset(is_local%wrap%FBfrac(n1), value=czero, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
       first_call = .false.

    endif

    !---------------------------------------
    ! Set 'lfrac' in FBFrac(complnd) - this might be overwritten later
    ! Set 'lfrin' in FBFrac(complnd)
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then
       call fldbun_getdata1d(is_local%wrap%FBImp(complnd,complnd) , 'Sl_lfrin', Sl_lfrin, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBFrac(complnd) , 'lfrac', lfrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (associated(lfrac)) then
          lfrac(:) = Sl_lfrin(:)
       end if
       call fldbun_getdata1d(is_local%wrap%FBFrac(complnd) , 'lfrin', lfrin, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (associated(lfrin)) then
          lfrin(:) = Sl_lfrin(:)
       end if
    end if

    !---------------------------------------
    ! Set 'ifrac' in FBFrac(compice)
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then
       ! Set 'ifrac' FBFrac(compice)
       call fldbun_getdata1d(is_local%wrap%FBImp(compice,compice), 'Si_imask', Si_imask, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBFrac(compice), 'ifrac', ifrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (associated(ifrac)) then
          ifrac(:) = Si_imask(:)
       end if
    end if

    !---------------------------------------
    ! Set 'ifrac' in FBFrac(compatm)
    !---------------------------------------

    if ( is_local%wrap%comp_present(compice) .and. &
         is_local%wrap%comp_present(compatm) .and. &
         is_local%wrap%med_coupling_active(compice,compatm)) then

       ! Set 'ifrac' in  FBFrac(compatm) - at this point this is the ice mask mapped to the atm mesh
       ! This maps the ice mask (which is the same as the ocean mask) to the atm mesh
       if (med_map_RH_is_created(is_local%wrap%RH(compice,compatm,:),mapfcopy, rc=rc)) then
          ! If ice and atm are on the same mesh - a redist route handle has already been created
          maptype = mapfcopy
       else
          if (coupling_mode(1:9) == 'ufs.nfrac' ) then
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
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compice), 'ifrac', field=field_src, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'ifrac', field=field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_map_field(field_src, field_dst, is_local%wrap%RH(compice,compatm,:), maptype, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Set 'ofrac' in FBFrac(compocn)
    !---------------------------------------

    if (is_local%wrap%comp_present(compocn)) then
       ! Set 'ofrac' in FBFrac(compocn)
       call fldbun_getdata1d(is_local%wrap%FBImp(compocn,compocn), 'So_omask', So_omask, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBFrac(compocn), 'ofrac', ofrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (associated(ofrac)) then
          ofrac(:) = So_omask(:)
       end if
    end if

    !---------------------------------------
    ! Set 'ofrac' in FBFrac(compatm)
    !---------------------------------------

    if ( is_local%wrap%comp_present(compocn) .and. &
         is_local%wrap%comp_present(compatm) .and. &
         is_local%wrap%med_coupling_active(compocn,compatm)) then

       ! Set 'ofrac' in FBFrac(compatm) - at this point this is the
       ! ocean mask mapped to the atm grid This is mapping the ocean mask to
       ! the atm grid - so in effect it is (1-land fraction) on the atm grid

       if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:),mapfcopy, rc=rc)) then
          ! If ocn and atm are on the same mesh - a redist route handle has already been created
          maptype = mapfcopy
       else
          if (coupling_mode(1:9) == 'ufs.nfrac' ) then
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
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compocn), fieldname='ofrac', field=field_src, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), fieldname='ofrac', field=field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_map_field(field_src, field_dst, is_local%wrap%RH(compocn,compatm,:), maptype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Set 'aofrac' in FBfrac(compatm) if available
       if ( fldbun_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ofrac', rc=rc)) then
          call fldbun_getdata1d(is_local%wrap%FBImp(compatm,compatm), 'Sa_ofrac', Sa_ofrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call fldbun_getdata1d(is_local%wrap%FBFrac(compatm), 'aofrac', aofrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (associated(ofrac)) then
             aofrac(:) = Sa_ofrac(:)
          end if
       end if
    end if

    !---------------------------------------
    ! Set 'lfrin' in FBFrac(compatm)
    ! ---------------------------------------

    if ( is_local%wrap%comp_present(compatm) .and. &
         is_local%wrap%comp_present(complnd) .and. &
         is_local%wrap%med_coupling_active(complnd,compatm)) then

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

      call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrin', field=field_src, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'lfrin', field=field_dst, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,compatm,:), maptype, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !---------------------------------------
    ! Set 'lfrac' in FBFrac(compatm)
    ! Reset 'ofrac' in FBFrac(compatm) if appropriate
    ! ---------------------------------------
    ! These should actually be mapo2a of ofrac and lfrac but we can't
    ! map lfrac from o2a due to masked mapping weights.  So we have to
    ! settle for a residual calculation that is truncated to zero to
    ! try to preserve "all ocean" cells.

    if (is_local%wrap%comp_present(compatm)) then

       if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compice)) then

          ! Ocean or ice is present
          call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (associated(lfrac)) then
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
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrac', field=field_src, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'lfrac', field=field_dst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,compatm,:), maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Reset ofrac in FBFrac(compatm)
          if (samegrid_atmlnd) then
            call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
            call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'lfrin', lfrac, rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (associated(ofrac)) then
            do n = 1,size(lfrac)
              ofrac(n) = 1.0_R8 - lfrac(n)
              if (abs(ofrac(n)) < eps_fraclim) then
                ofrac(n) = 0.0_R8
              end if
            end do
          end if

        end if
    end if

    !---------------------------------------
    ! Reset 'lfrac' in FBFrac(complnd) if appropriate
    !---------------------------------------

    if ( is_local%wrap%comp_present(complnd) .and. &
         is_local%wrap%med_coupling_active(complnd,compatm)) then

       ! If lnd -> atm coupling is active - map 'lfrac' from FBFrac(compatm) to FBFrac(complnd)
       ! Note that if the atmosphere is absent, then simply set fractions_l(lfrac) = fractions_l(lfrin)
       ! from above

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
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'lfrac', field=field_src, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrac', field=field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_map_field(field_src, field_dst, is_local%wrap%RH(compatm,complnd,:), maptype, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Set 'rfrac' and 'lfrac' for FBFrac(comprof)
    !---------------------------------------

    if (is_local%wrap%comp_present(comprof)) then

       ! Set 'rfrac' in FBFrac(comprof)
       if ( fldbun_fldchk(is_local%wrap%FBfrac(comprof)       , 'rfrac', rc=rc) .and. &
            fldbun_fldchk(is_local%wrap%FBImp(comprof,comprof), 'frac' , rc=rc)) then
          call fldbun_getdata1d(is_local%wrap%FBImp(comprof,comprof), 'frac', frac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call fldbun_getdata1d(is_local%wrap%FBfrac(comprof), 'rfrac', rfrac, rc)
          if (associated(rfrac)) then
             rfrac(:) = frac(:)
          endif
       else
          ! Set 'rfrac' in FBfrac(comprof) to 1.
          call fldbun_getdata1d(is_local%wrap%FBfrac(comprof), 'rfrac', rfrac, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (associated(rfrac)) then
             rfrac(:) = 1.0_R8
          endif
       endif

       ! Set 'lfrac' and 'lfrin' in FBFrac(comprof)
       if (is_local%wrap%comp_present(complnd)) then
          maptype = mapconsd
          if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,comprof,:),maptype, rc=rc)) then
             call med_map_routehandles_init( complnd, comprof, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,comprof), &
                  mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrac', field=field_src, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(comprof), 'lfrac', field=field_dst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,comprof,:), maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrin', field=field_src, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(comprof), 'lfrin', field=field_dst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,comprof,:), maptype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    !---------------------------------------
    ! Set 'gfrac', 'lfrac' and 'lfrin' in FBFrac(compglc)
    !---------------------------------------

    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%comp_present(compglc(ns))) then

          ! Set 'gfrac' in FBFrac(compglc(ns))
          if ( fldbun_fldchk(is_local%wrap%FBfrac(compglc(ns))            , 'gfrac', rc=rc) .and. &
               fldbun_fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'frac' , rc=rc)) then
             call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), 'frac', frac, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call fldbun_getdata1d(is_local%wrap%FBfrac(compglc(ns)), 'gfrac', gfrac, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (associated(gfrac)) then
                gfrac(:) = frac(:)
             endif
          else
             ! Set 'gfrac' in FBfrac(compglc(ns)) to 1.
             call fldbun_getdata1d(is_local%wrap%FBfrac(compglc(ns)), 'gfrac', gfrac, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (associated(gfrac)) then
                gfrac(:) = 1.0_R8
             endif
          endif

          ! Set 'lfrac' and 'lfrin' in FBFrac(compglc(ns))
          if ( is_local%wrap%comp_present(complnd) .and. is_local%wrap%med_coupling_active(complnd,compglc(ns))) then
             maptype = mapconsd
             if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,compglc(ns),:),maptype, rc=rc)) then
                call med_map_routehandles_init( complnd, compglc(ns), &
                     FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                     FBDst=is_local%wrap%FBImp(complnd,compglc(ns)), &
                     mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrac', field=field_src, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compglc(ns)), 'lfrac', field=field_dst, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,compglc(ns),:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(complnd), 'lfrin', field=field_src, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compglc(ns)), 'lfrin', field=field_dst, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_map_field(field_src, field_dst, is_local%wrap%RH(complnd,compglc(ns),:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif
    end do

    !---------------------------------------
    ! Set 'wfrac' for FBFrac(compwav)
    !---------------------------------------

    if (is_local%wrap%comp_present(compwav)) then
       ! Set 'wfrac' in FBfrac(compwav) to 1.
       call fldbun_getdata1d(is_local%wrap%FBfrac(compwav), 'wfrac', wfrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (associated(wfrac)) then
          wfrac(:) = 1.0_R8
       endif
    endif

    !---------------------------------------
    ! Create route handles ocn<->ice if not created
    !---------------------------------------

    if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
       if (.not. med_map_RH_is_created(is_local%wrap%RH(compice,compocn,:),mapfcopy, rc=rc)) then
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compice,compocn))) then
             call fldbun_init(is_local%wrap%FBImp(compice,compocn), is_local%wrap%flds_scalar_name, &
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
             call fldbun_init(is_local%wrap%FBImp(compocn,compice), is_local%wrap%flds_scalar_name, &
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
             call fldbun_diagnose(is_local%wrap%FBfrac(n), trim(subname) // trim(compname(n)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_fraction_init

  !================================================================================================
  subroutine med_fraction_set(gcomp, rc)

    ! Update time varying fractions

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_internalstate_mod , only : compatm, compocn, compice, compname
    use med_internalstate_mod , only : mapfcopy, mapconsd, mapnstod_consd
    use med_internalstate_mod , only : coupling_mode
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_RH_is_created
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    real(r8), pointer          :: ifrac(:)
    real(r8), pointer          :: ofrac(:)
    real(r8), pointer          :: aofrac(:)
    real(r8), pointer          :: Si_ifrac(:)
    real(r8), pointer          :: Si_imask(:)
    real(r8), pointer          :: Sa_ofrac(:)
    type(ESMF_Field)           :: field_src
    type(ESMF_Field)           :: field_dst
    integer                    :: n
    integer                    :: maptype
    character(len=*),parameter :: subname=' (med_fraction_set)'
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

       call fldbun_getdata1d(is_local%wrap%FBImp(compice,compice), 'Si_ifrac', Si_ifrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBImp(compice,compice), 'Si_imask', Si_imask, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBFrac(compice), 'ifrac', ifrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBFrac(compice), 'ofrac', ofrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Note that Si_imask and Si_ifrac are the same when the ocn/ice mask is either 0 or 1
       ! However, in the case where the atm/lnd/ice/ocn are all on the atm grid, Si_imask is
       ! Si_imask is the the model mask mapped to the atm grid (e.g. a gx1v7 mask mapped to the atm grid)
       ! The model mask is normally assumed to be an selected ocean mask from a fully coupled run
       ! So in it is (1-land fraction) on the atm grid

       ! set ifrac
       if (associated(ifrac)) then
          ifrac(:) = Si_ifrac(:) * Si_imask(:)
       endif

       ! set ofrac = Si_imask - ifrac
       if (associated(ofrac)) then
          ofrac(:) = Si_imask(:) - ifrac(:)
       endif

       ! -------------------------------------------
       ! Set FBfrac(compocn)
       ! -------------------------------------------

       ! The following is just a redistribution from FBFrac(compice)

       call t_startf('MED:'//trim(subname)//' fbfrac(compocn)')
       if (is_local%wrap%comp_present(compocn)) then
          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compocn)
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compice), 'ifrac', field=field_src, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compocn), 'ifrac', field=field_dst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_map_field(field_src, field_dst, is_local%wrap%RH(compice,compocn,:), mapfcopy, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compocn)
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compice), 'ofrac', field=field_src, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compocn), 'ofrac', field=field_dst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_map_field(field_src, field_dst, is_local%wrap%RH(compice,compocn,:), mapfcopy, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       call t_stopf('MED:'//trim(subname)//' fbfrac(compocn)')

       ! -------------------------------------------
       ! Set FBfrac(compatm)
       ! -------------------------------------------

       if (is_local%wrap%comp_present(compatm)) then

          call t_startf('MED:'//trim(subname)//' fbfrac(compatm)')
          ! Determine maptype
          if (coupling_mode(1:9) == 'ufs.nfrac' ) then
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
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compice), 'ifrac', field=field_src, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'ifrac', field=field_dst, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_map_field(field_src, field_dst, is_local%wrap%RH(compice,compatm,:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compatm)
          if (is_local%wrap%med_coupling_active(compocn,compatm)) then
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compice), 'ofrac', field=field_src, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(is_local%wrap%FBfrac(compatm), 'ofrac', field=field_dst, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_map_field(field_src, field_dst, is_local%wrap%RH(compice,compatm,:), maptype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          ! Set 'aofrac' from FBImp(compatm) to FBfrac(compatm) if available
          if ( fldbun_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ofrac', rc=rc)) then
             call fldbun_getdata1d(is_local%wrap%FBImp(compatm,compatm), 'Sa_ofrac', Sa_ofrac, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call fldbun_getdata1d(is_local%wrap%FBFrac(compatm), 'aofrac', aofrac, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             aofrac(:) = Sa_ofrac(:)
          end if
       end if ! end of if present compatm
       call t_stopf('MED:'//trim(subname)//' fbfrac(compatm)')

    end if ! end of if present compatm

    !---------------------------------------
    ! Diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       do n = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call fldbun_diagnose(is_local%wrap%FBfrac(n), trim(subname) // trim(compname(n))//' frac', rc=rc)
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
