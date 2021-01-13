module med_phases_prep_rof_mod

  !-----------------------------------------------------------------------------
  ! Create rof export fields
  ! - accumulate import lnd fields on the land grid that are sent to rof
  !   this will be done in med_phases_prep_rof_accum
  ! - time avergage accumulated import lnd fields when necessary
  !   map the time averaged accumulated lnd fields to the rof grid
  !   merge the mapped lnd fields to create FBExp(comprof)
  !   this will be done in med_phases_prep_rof_avg
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_FieldBundle, ESMF_Field
  use esmFlds               , only : ncomps, complnd, comprof, compname, mapconsf, mapconsd, mapfcopy
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_constants_mod     , only : czero            => med_constants_czero
  use med_utils_mod         , only : chkerr           => med_utils_chkerr
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use med_methods_mod       , only : fldbun_getdata2d => med_methods_FB_getdata2d
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_reset     => med_methods_FB_reset
  use med_methods_mod       , only : fldbun_average   => med_methods_FB_average
  use med_methods_mod       , only : field_getdata2d  => med_methods_Field_getdata2d
  use med_methods_mod       , only : field_getdata1d  => med_methods_Field_getdata1d
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_rof        ! called by run sequence
  public  :: med_phases_prep_rof_accum  ! called by med_phases_post_lnd

  private :: med_phases_prep_rof_irrig

  ! the following are needed for lnd2rof irrigation
  type(ESMF_Field) :: field_lndVolr
  type(ESMF_Field) :: field_rofVolr
  type(ESMF_Field) :: field_lndIrrig
  type(ESMF_Field) :: field_rofIrrig
  type(ESMF_Field) :: field_lndIrrig0
  type(ESMF_Field) :: field_rofIrrig0
  type(ESMF_Field) :: field_lfrac_rof

  character(len=*), parameter :: volr_field             = 'Flrr_volrmch'
  character(len=*), parameter :: irrig_flux_field       = 'Flrl_irrig'
  character(len=*), parameter :: irrig_normalized_field = 'Flrl_irrig_normalized'
  character(len=*), parameter :: irrig_volr0_field      = 'Flrl_irrig_volr0     '

  ! the following are the fields that will be accumulated from the land
  character(CS) :: lnd2rof_flds(6) = (/'Flrl_rofsur','Flrl_rofgwl','Flrl_rofsub', &
                                       'Flrl_rofdto','Flrl_rofi  ','Flrl_irrig '/)

  integer :: maptype_lnd2rof
  integer :: maptype_rof2lnd

  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_prep_rof_accum(gcomp, rc)

    !------------------------------------
    ! Carry out fast accumulation for the river (rof) component
    ! Accumulation and averaging is done on the land input on the land grid for the fields that will
    ! will be sent to the river component
    ! Mapping from the land to the rof grid is then done with the time averaged fields
    !------------------------------------

    use NUOPC , only : NUOPC_IsConnected
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_FieldBundleGet, ESMF_StateIsCreated, ESMF_StateGet
    use ESMF  , only : ESMF_FieldBundleIsCreated, ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp

    integer, intent(out) :: rc

    ! local variables
    type(InternalState)       :: is_local
    integer                   :: i,j,n,ncnt
    integer                   :: fieldCount
    integer                   :: ungriddedUBound(1)
    logical                   :: exists
    real(r8), pointer         :: dataptr1d(:) => null()
    real(r8), pointer         :: dataptr2d(:,:) => null()
    real(r8), pointer         :: dataptr1d_accum(:) => null()
    real(r8), pointer         :: dataptr2d_accum(:,:) => null()
    type(ESMF_Field)          :: lfield
    type(ESMF_Field)          :: lfield_accum
    type(ESMF_Field), pointer :: fieldlist(:) => null()
    type(ESMF_Field), pointer :: fieldlist_accum(:) => null()
    character(CL), pointer    :: lfieldnamelist(:) => null()
    character(len=*), parameter :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Accumulate lnd input on lnd grid for fields that will be sent to rof
    do n = 1,size(lnd2rof_flds)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
            isPresent=exists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBImpaccum(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
               field=lfield_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (ungriddedUBound(1) > 0) then
             call field_getdata2d(lfield, dataptr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call field_getdata2d(lfield_accum, dataptr2d_accum, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr2d_accum(:,:) = dataptr2d_accum(:,:) + dataptr2d(:,:)
          else
             call field_getdata1d(lfield, dataptr1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call field_getdata1d(lfield_accum, dataptr1d_accum, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr1d_accum(:) = dataptr1d_accum(:) + dataptr1d(:)
          end if
       end if
    end do

    ! Accumulate counter
    is_local%wrap%FBImpAccumCnt(complnd) = is_local%wrap%FBImpAccumCnt(complnd) + 1

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBImpAccum(complnd,complnd), &
            string=trim(subname)//' FBImpAccum(complnd,complnd) ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof_accum

  !===============================================================================
  subroutine med_phases_prep_rof(gcomp, rc)

    !------------------------------------
    ! Prepare the ROF export Fields from the mediator
    !------------------------------------

    use NUOPC             , only : NUOPC_IsConnected
    use ESMF              , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF              , only : ESMF_FieldBundleGet, ESMF_FieldGet
    use ESMF              , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use esmFlds           , only : fldListTo
    use med_map_mod       , only : med_map_field_packed
    use med_merge_mod     , only : med_merge_auto
    use med_constants_mod , only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)       :: is_local
    integer                   :: i,j,n,n1,ncnt
    integer                   :: count
    logical                   :: exists
    real(r8), pointer         :: dataptr(:) => null()
    real(r8), pointer         :: dataptr1d(:) => null()
    real(r8), pointer         :: dataptr2d(:,:) => null()
    type(ESMF_Field)          :: field_irrig_flux
    integer                   :: fieldcount
    type(ESMF_Field)          :: lfield
    type(ESMF_Field), pointer :: fieldlist(:) => null()
    integer                   :: ungriddedUBound(1)
    character(CL), pointer    :: lfieldnamelist(:) => null()
    character(len=*),parameter  :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Average import from land accumuled FB
    !---------------------------------------

    count = is_local%wrap%FBImpAccumCnt(complnd)
    if (count == 0) then
       if (mastertask) then
          write(logunit,'(a)')trim(subname)//'accumulation count for land input averging to river is 0 '// &
               ' accumulation field is set to zero'
       end if
    end if

    do n = 1,size(lnd2rof_flds)
       call ESMF_FieldBundleGet(is_local%wrap%FBImpAccum(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
            isPresent=exists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          call ESMF_FieldBundleGet(is_local%wrap%FBImpAccum(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (ungriddedUBound(1) > 0) then
             call field_getdata2d(lfield, dataptr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (count == 0) then
                dataptr2d(:,:) = czero
             else
                dataptr2d(:,:) = dataptr2d(:,:) / real(count, r8)
             end if
          else
             call field_getdata1d(lfield, dataptr1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (count == 0) then
                dataptr1d(:) = czero
             else
                dataptr1d(:) = dataptr1d(:) / real(count, r8)
             end if
          end if
       end if
    end do

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBImpAccum(complnd,complnd), &
            string=trim(subname)//' FBImpAccum(complnd,complnd) after avg ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Map to create FBImpAccum(complnd,comprof)
    !---------------------------------------

    ! The following assumes that only land import fields are needed to create the
    ! export fields for the river component and that ALL mappings are done with mapconsf

    if (is_local%wrap%med_coupling_active(complnd,comprof)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImpAccum(complnd,complnd), &
            FBDst=is_local%wrap%FBImpAccum(complnd,comprof), &
            FBFracSrc=is_local%wrap%FBFrac(complnd), &
            field_normOne=is_local%wrap%field_normOne(complnd,comprof,:), &
            packed_data=is_local%wrap%packed_data(complnd,comprof,:), &
            routehandles=is_local%wrap%RH(complnd,comprof,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call fldbun_diagnose(is_local%wrap%FBImpAccum(complnd,comprof), &
               string=trim(subname)//' FBImpAccum(complnd,comprof) after map ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Reset the irrig_flux_field with the map_lnd2rof_irrig calculation below if appropriate
       if ( NUOPC_IsConnected(is_local%wrap%NStateImp(complnd), fieldname=trim(irrig_flux_field))) then
          call med_phases_prep_rof_irrig( gcomp, rc=rc )
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          ! This will ensure that no irrig is sent from the land
          call fldbun_getdata1d(is_local%wrap%FBImpAccum(complnd,comprof), irrig_flux_field, dataptr, rc)
          dataptr(:) = czero
       end if
    endif

    !---------------------------------------
    ! auto merges to create FBExp(comprof)
    !---------------------------------------

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBFrac(comprof), &
            string=trim(subname)//' FBFrac(comprof) before merge ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call med_merge_auto(comprof, &
         is_local%wrap%med_coupling_active(:,comprof), &
         is_local%wrap%FBExp(comprof), &
         is_local%wrap%FBFrac(comprof), &
         is_local%wrap%FBImpAccum(:,comprof), &
         fldListTo(comprof), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBExp(comprof), &
            string=trim(subname)//' FBexp(comprof) ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! zero accumulator and FBAccum
    !---------------------------------------

    ! zero counter
    is_local%wrap%FBImpAccumCnt(complnd) = 0

    ! zero lnd2rof fields in FBImpAccum
    do n = 1,size(lnd2rof_flds)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
            isPresent=exists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          call ESMF_FieldBundleGet(is_local%wrap%FBImpaccum(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (ungriddedUBound(1) > 0) then
             call field_getdata2d(lfield, dataptr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr2d(:,:) = czero
          else
             call field_getdata1d(lfield, dataptr1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr1d(:) = czero
          end if
       end if
    end do

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof

  !===============================================================================
  subroutine med_phases_prep_rof_irrig(gcomp, rc)

    !---------------------------------------------------------------
    ! Description
    ! Do custom mapping for the irrigation flux, from land -> rof.
    !
    ! The basic idea is that we want to pull irrigation out of ROF cells proportionally to
    ! the river volume (volr) in each cell. This is important in cases where the various
    ! ROF cells overlapping a CTSM cell have very different volr: If we didn't do this
    ! volr-normalized remapping, we'd try to extract the same amount of water from each
    ! of the ROF cells, which would be more likely to have withdrawals exceeding
    ! available volr.
    !
    ! (Both RTM and MOSART have code to handle excess withdrawals by pulling the excess
    ! directly out of the ocean. We'd like to avoid resorting to this if possible.
    !
    ! This mapping works by:
    ! (1) Normalizing the land's irrigation flux by volr
    ! (2) Mapping this volr-normalized flux to the rof grid
    ! (3) Converting the mapped, volr-normalized flux back to a normal
    !     (non-volr-normalized) flux on the rof grid.
    !---------------------------------------------------------------

    use ESMF        , only : ESMF_GridComp, ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
    use ESMF        , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldIsCreated
    use ESMF        , only : ESMF_Mesh, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT
    use ESMF        , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF        , only : ESMF_LOGMSG_INFO, ESMF_LogWrite, ESMF_LOGMSG_ERROR
    use med_map_mod , only : med_map_rh_is_created, med_map_field, med_map_field_normalized

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                   :: r,l
    type(InternalState)       :: is_local
    integer                   :: fieldcount
    type(ESMF_Field)          :: field_import_rof
    type(ESMF_Field)          :: field_import_lnd
    type(ESMF_Field)          :: field_irrig_flux
    type(ESMF_Field)          :: field_lfrac_lnd
    type(ESMF_Field), pointer :: fieldlist_lnd(:) => null()
    type(ESMF_Field), pointer :: fieldlist_rof(:) => null()
    type(ESMF_Mesh)           :: lmesh_lnd
    type(ESMF_Mesh)           :: lmesh_rof
    real(r8), pointer         :: volr_l(:) => null()
    real(r8), pointer         :: volr_r(:) => null()
    real(r8), pointer         :: volr_r_import(:) => null()
    real(r8), pointer         :: irrig_normalized_l(:) => null()
    real(r8), pointer         :: irrig_normalized_r(:) => null()
    real(r8), pointer         :: irrig_volr0_l(:) => null()
    real(r8), pointer         :: irrig_volr0_r(:) => null()
    real(r8), pointer         :: irrig_flux_l(:) => null()
    real(r8), pointer         :: irrig_flux_r(:) => null()
    character(len=*), parameter :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof_irrig)'
    !---------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (med_map_RH_is_created(is_local%wrap%RH(complnd,comprof,:),mapconsf, rc=rc)) then
       maptype_lnd2rof = mapconsf
    else if ( med_map_RH_is_created(is_local%wrap%RH(complnd,comprof,:),mapfcopy, rc=rc)) then
       maptype_lnd2rof = mapfcopy
    else
       call ESMF_LogWrite(trim(subname)//&
            ": ERROR conservative or redist route handles not created for lnd->rof mapping", &
            ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    if (med_map_RH_is_created(is_local%wrap%RH(comprof,complnd,:),mapconsf, rc=rc)) then
       maptype_rof2lnd = mapconsf
    else if ( med_map_RH_is_created(is_local%wrap%RH(comprof,complnd,:),mapfcopy, rc=rc)) then
       maptype_rof2lnd = mapfcopy
    else
       call ESMF_LogWrite(trim(subname)//&
            ": ERROR conservative or redist route handles not created for rof->lnd mapping", &
            ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! ------------------------------------------------------------------------
    ! Initialize module field bundles if not already initialized
    ! ------------------------------------------------------------------------

    if (.not. ESMF_FieldIsCreated(field_lndVolr)   .and. &
        .not. ESMF_FieldIsCreated(field_rofVolr)   .and. &
        .not. ESMF_FieldIsCreated(field_lndIrrig)  .and. &
        .not. ESMF_FieldIsCreated(field_rofIrrig)  .and. &
        .not. ESMF_FieldIsCreated(field_lndIrrig0) .and. &
        .not. ESMF_FieldIsCreated(field_rofIrrig0) .and. &
        .not. ESMF_FieldIsCreated(field_lfrac_rof)) then

       ! get fields in source and destination field bundles
       call fldbun_getmesh(is_local%wrap%FBImp(complnd,complnd), lmesh_lnd, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getmesh(is_local%wrap%FBImp(comprof,comprof), lmesh_rof, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       field_lndVolr = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       field_rofVolr = ESMF_FieldCreate(lmesh_rof, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       field_lndIrrig = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       field_rofIrrig = ESMF_FieldCreate(lmesh_rof, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       field_lndIrrig0 = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       field_rofIrrig0 = ESMF_FieldCreate(lmesh_rof, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       field_lfrac_rof = ESMF_FieldCreate(lmesh_rof, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    ! ------------------------------------------------------------------------
    ! 1) Create volr_l: Adjust volr_r, and map it to the land grid
    ! ------------------------------------------------------------------------

    ! Treat any rof point with volr < 0 as if it had volr = 0. Negative volr values can
    ! arise in RTM. This fix is needed to avoid mapping negative irrigation to those
    ! cells: while conservative, this would be unphysical (it would mean that irrigation
    ! actually adds water to those cells).

    ! Create volr_r
    call fldbun_getdata1d(is_local%wrap%FBImp(comprof,comprof), trim(volr_field), volr_r_import, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(field_rofVolr, volr_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do r = 1, size(volr_r)
       if (volr_r_import(r) < 0._r8) then
          volr_r(r) = 0._r8
       else
          volr_r(r) = volr_r_import(r)
       end if
    end do

    ! Map volr_r to volr_l (rof->lnd) using conservative mapping without any fractional weighting
    call med_map_field( &
         field_src=field_rofVolr, &
         field_dst=field_lndVolr, &
         routehandles=is_local%wrap%RH(comprof,complnd,:), &
         maptype=maptype_rof2lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get volr_l
    call field_getdata1d(field_lndVolr, volr_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! (2) Determine irrigation from land on land grid normalized by volr_l
    ! ------------------------------------------------------------------------

    ! In order to avoid possible divide by 0, as well as to handle non-sensical negative
    ! volr on the land grid, we divide the land's irrigation flux into two separate flux
    ! components:
    ! - a component where we have positive volr on the land grid (put in
    !   irrig_normalized_l, which is mapped using volr-normalization)
    ! - a component where we have zero or negative volr on the land
    !   grid (put in irrig_volr0_l, which is mapped as a standard flux).
    ! We then remap both of these components to the rof grid, and then
    ! finally add the two components to determine the total irrigation
    ! flux on the rof grid.

    ! First extract accumulated irrigation flux from land
    call fldbun_getdata1d(is_local%wrap%FBImpAccum(complnd,complnd), trim(irrig_flux_field), irrig_flux_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Fill in values for irrig_normalized_l and irrig_volr0_l
    call field_getdata1d(field_lndIrrig, irrig_normalized_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(field_lndIrrig0, irrig_volr0_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do l = 1, size(volr_l)
       if (volr_l(l) > 0._r8) then
          irrig_normalized_l(l) = irrig_flux_l(l) / volr_l(l)
          irrig_volr0_l(l)      = 0._r8
       else
          irrig_normalized_l(l) = 0._r8
          irrig_volr0_l(l)      = irrig_flux_l(l)
       end if
    end do

    ! ------------------------------------------------------------------------
    ! (3) Map normalized irrigation from land to rof grid and
    !     convert to a total irrigation flux on the ROF grid
    ! ------------------------------------------------------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(complnd), 'lfrac', field=field_lfrac_lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call med_map_field_normalized(  &
         field_src=field_lndIrrig, &
         field_dst=field_rofIrrig, &
         routehandles=is_local%wrap%RH(complnd,comprof,:), &
         maptype=maptype_lnd2rof, &
         field_normsrc=field_lfrac_lnd, &
         field_normdst=field_lfrac_rof, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call med_map_field_normalized(  &
         field_src=field_lndIrrig0, &
         field_dst=field_rofIrrig0, &
         routehandles=is_local%wrap%RH(complnd,comprof,:), &
         maptype=maptype_lnd2rof, &
         field_normsrc=field_lfrac_lnd, &
         field_normdst=field_lfrac_rof, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Convert to a total irrigation flux on the ROF grid, and put this in the pre-merge FBImpAccum(complnd,comprof)
    call field_getdata1d(field_rofIrrig, irrig_normalized_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(field_rofIrrig0, irrig_volr0_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(is_local%wrap%FBImpAccum(complnd,comprof), trim(irrig_flux_field), irrig_flux_r, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(field_rofIrrig0, irrig_volr0_r, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do r = 1, size(irrig_flux_r)
       irrig_flux_r(r) = (irrig_normalized_r(r) * volr_r(r)) + irrig_volr0_r(r)
    end do

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof_irrig

end module med_phases_prep_rof_mod
