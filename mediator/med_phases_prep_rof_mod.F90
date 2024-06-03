module med_phases_prep_rof_mod

  !-----------------------------------------------------------------------------
  ! Create rof export fields
  ! - accumulate import lnd fields on the land grid that are sent to rof
  !   - done in med_phases_prep_rof_accum
  ! - time avergage accumulated import lnd fields on lnd grid when necessary and
  !   then map the time averaged accumulated lnd fields to the rof grid
  !   and then merge the mapped lnd fields to create FBExp(comprof)
  !   - done in med_phases_prep_rof_avg
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_FieldBundle, ESMF_Field
  use med_internalstate_mod , only : complnd, compglc, comprof, mapconsf, mapfcopy
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_constants_mod     , only : czero            => med_constants_czero
  use med_utils_mod         , only : chkerr           => med_utils_chkerr
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_reset     => med_methods_FB_reset
  use med_methods_mod       , only : fldbun_average   => med_methods_FB_average
  use med_methods_mod       , only : field_getdata1d  => med_methods_Field_getdata1d
  use med_methods_mod       , only : fldbun_fldchk    => med_methods_FB_fldchk
  use med_methods_mod       , only : FB_check_for_nans => med_methods_FB_check_for_nans
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_rof_init   ! called from med.F90
  public  :: med_phases_prep_rof_accum  ! called by med_phases_post_lnd.F90
  public  :: med_phases_prep_rof        ! called by run sequence

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

  ! the following are the fields that will be accumulated from the land and are derived from fldlistTo(comprof)
  character(CS), allocatable :: lnd2rof_flds(:)

  integer :: maptype_lnd2rof
  integer :: maptype_rof2lnd

  ! Accumulation to river field bundles - accumulation is done on the land mesh and then averaged and mapped to the
  ! rof mesh
  integer               , public :: lndAccum2rof_cnt
  type(ESMF_FieldBundle), public :: FBlndAccum2rof_l
  type(ESMF_FieldBundle), public :: FBlndAccum2rof_r

  character(len=9) :: fldnames_fr_glc(2) = (/'Fgrg_rofl', 'Fgrg_rofi'/)

  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_prep_rof_init(gcomp, rc)

    !---------------------------------------
    ! Create module field bundles FBlndAccum2rof_l and FBlndAccum2rof_r
    ! land accumulation on both complnd and comprof meshes
    !---------------------------------------

    use ESMF        , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF        , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
    use ESMF        , only : ESMF_Mesh, ESMF_MESHLOC_ELEMENT
    use ESMF        , only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleGet, ESMF_FieldBundleAdd
    use ESMF        , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF        , only : ESMF_TYPEKIND_R8
    use esmFlds     , only : med_fldList_GetfldListFr, med_fldList_GetfldlistTo, med_fldlist_GetNumFlds, med_fld_getFldInfo
    use esmFlds     , only : med_fldList_type, med_fldList_entry_type
    use med_map_mod , only : med_map_packed_field_create

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n, nflds
    type(ESMF_Mesh)     :: mesh_l
    type(ESMF_Mesh)     :: mesh_r
    type(ESMF_Field)    :: lfield
    type(med_fldList_type), pointer :: fldList
    type(med_fldList_entry_type), pointer :: fldptr
    character(len=CS)  :: fldname
    character(len=CS), allocatable  :: fldnames_temp(:)
    character(len=*),parameter  :: subname=' (med_phases_prep_rof_init) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lnd2rof_flds (module variable) - note that fldListTo is set in esmFldsExchange_cesm.F90
    ! Remove scalar field from lnd2rof_flds
    fldList => med_fldList_GetfldlistTo(comprof)
    nflds = med_fldlist_getnumflds(fldList)
    allocate(fldnames_temp(nflds))
    fldptr => fldList%fields
    n = 0
    do while(associated(fldptr))
       call med_fld_GetFldInfo(fldptr, stdname=fldname)
       if (trim(fldname) .ne. trim(is_local%wrap%flds_scalar_name)) then
          n = n+1
          fldnames_temp(n) = fldname
       endif
       fldptr => fldptr%next
    enddo
    allocate(lnd2rof_flds(n))
    lnd2rof_flds = fldnames_temp(1:n)
    deallocate(fldnames_temp)

    ! Get lnd and rof meshes
    call fldbun_getmesh(is_local%wrap%FBImp(complnd,complnd), mesh_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getmesh(is_local%wrap%FBImp(complnd,comprof), mesh_r, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create module field bundle FBlndAccum2rof_l on land mesh and  FBlndAccum2rof_r on rof mesh
    FBlndAccum2rof_l = ESMF_FieldBundleCreate(name='FBlndAccum2rof_l', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    FBlndAccum2rof_r = ESMF_FieldBundleCreate(name='FBlndAccum2rof_r', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,size(lnd2rof_flds)
       lfield = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, name=lnd2rof_flds(n), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBlndAccum2rof_l, (/lfield/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(trim(subname)//' adding field '//trim(lnd2rof_flds(n))//' to FBLndAccum2rof_l', &
            ESMF_LOGMSG_INFO)
       lfield = ESMF_FieldCreate(mesh_r, ESMF_TYPEKIND_R8, name=lnd2rof_flds(n), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBlndAccum2rof_r, (/lfield/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(trim(subname)//' adding field '//trim(lnd2rof_flds(n))//' to FBLndAccum2rof_r', &
            ESMF_LOGMSG_INFO)
    end do

    ! Initialize field bundles and accumulation count

    call fldbun_reset(FBlndAccum2rof_l, value=0.0_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call fldbun_reset(FBlndAccum2rof_r, value=0.0_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lndAccum2rof_cnt = 0

    fldList => med_fldList_GetFldListFr(complnd)
    ! Create packed mapping from rof->lnd

    call med_map_packed_field_create(destcomp=comprof, &
         flds_scalar_name=is_local%wrap%flds_scalar_name, &
         fieldsSrc=fldList, &
         FBSrc=FBLndAccum2rof_l, FBDst=FBLndAccum2rof_r, &
         packed_data=is_local%wrap%packed_data(complnd,comprof,:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_prep_rof_init

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
    integer                   :: n
    logical                   :: exists
    real(r8), pointer         :: dataptr1d(:)
    real(r8), pointer         :: dataptr1d_accum(:)
    type(ESMF_Field)          :: lfield
    type(ESMF_Field)          :: lfield_accum
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
          call ESMF_FieldBundleGet(FBlndAccum2rof_l, fieldName=trim(lnd2rof_flds(n)), &
               field=lfield_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(lfield, dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(lfield_accum, dataptr1d_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_accum(:) = dataptr1d_accum(:) + dataptr1d(:)
       end if
    end do

    ! Accumulate counter
    lndAccum2rof_cnt =  lndAccum2rof_cnt + 1

    if (dbug_flag > 1) then
       call fldbun_diagnose(FBlndAccum2rof_l, string=trim(subname)//' FBlndAccum2rof_l accum', rc=rc)
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
    use esmFlds           , only : med_fldList_GetfldListTo, med_fldList_type
    use med_map_mod       , only : med_map_field_packed
    use med_merge_mod     , only : med_merge_auto
    use med_constants_mod , only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)       :: is_local
    integer                   :: n,ns,nf
    integer                   :: count
    logical                   :: exists
    real(r8), pointer         :: dataptr_in(:)
    real(r8), pointer         :: dataptr_out(:)
    type(ESMF_Field)          :: lfield
    type(med_fldList_type), pointer :: fldList
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

    fldList => med_fldList_GetfldListTo(comprof)
    !---------------------------------------
    ! Average import from land accumuled FB
    !---------------------------------------

    count = lndAccum2rof_cnt
    if (count == 0) then
       if (maintask) then
          write(logunit,'(a)')trim(subname)//'accumulation count for land input averging to river is 0 '// &
               ' accumulation field is set to zero'
       end if
    end if

    do n = 1,size(lnd2rof_flds)
       call ESMF_FieldBundleGet(FBlndAccum2rof_l, fieldName=trim(lnd2rof_flds(n)), isPresent=exists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          call ESMF_FieldBundleGet(FBlndAccum2rof_l, fieldName=trim(lnd2rof_flds(n)), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(lfield, dataptr_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (count == 0) then
             dataptr_out(:) = czero
          else
             dataptr_out(:) = dataptr_out(:) / real(count, r8)
          end if
       end if
    end do

    if (dbug_flag > 1) then
       call fldbun_diagnose(FBlndAccum2rof_l, string=trim(subname)//' FBlndAccum2rof_l after avg ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Map to create FBlndAccum2rof_r
    !---------------------------------------

    ! The following assumes that only land import fields are needed to create the
    ! export fields for the river component and that ALL mappings are done with mapconsf

    call med_map_field_packed( FBSrc=FBlndAccum2rof_l, FBDst=FBlndAccum2rof_r, &
         FBFracSrc=is_local%wrap%FBFrac(complnd), &
         field_normOne=is_local%wrap%field_normOne(complnd,comprof,:), &
         packed_data=is_local%wrap%packed_data(complnd,comprof,:), &
         routehandles=is_local%wrap%RH(complnd,comprof,:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call fldbun_diagnose(FBlndAccum2rof_r, string=trim(subname)//' FBlndAccum2rof_r after map ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Reset the irrig_flux_field with the map_lnd2rof_irrig calculation below if appropriate
    if ( NUOPC_IsConnected(is_local%wrap%NStateImp(complnd), fieldname=trim(irrig_flux_field))) then
       call med_phases_prep_rof_irrig( gcomp, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       ! This will ensure that no irrig is sent from the land
       call fldbun_getdata1d(FBlndAccum2rof_r, irrig_flux_field, dataptr_out, rc)
       dataptr_out(:) = czero
    end if

    !---------------------------------------
    ! create FBExp(comprof)
    !---------------------------------------

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBFrac(comprof), &
            string=trim(subname)//' FBFrac(comprof) before merge ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! data coming from FBlndAccum2rof_r
    call med_merge_auto(compsrc=complnd, FBout=is_local%wrap%FBExp(comprof), &
         FBfrac=is_local%wrap%FBFrac(comprof), FBin=FBlndAccum2rof_r, fldListTo=fldList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! custom merge for glc->rof
    ! glc->rof is mapped in med_phases_post_glc
    do ns = 1,is_local%wrap%num_icesheets
      if (is_local%wrap%med_coupling_active(compglc(ns),comprof)) then
        do nf = 1,size(fldnames_fr_glc)
          if ( fldbun_fldchk(is_local%wrap%FBImp(compglc(ns),comprof), fldnames_fr_glc(nf), rc=rc) .and. &
               fldbun_fldchk(is_local%wrap%FBExp(comprof), fldnames_fr_glc(nf), rc=rc) ) then
            call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),comprof), &
                 trim(fldnames_fr_glc(nf)), dataptr_in, rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return
            call fldbun_getdata1d(is_local%wrap%FBExp(comprof), &
                 trim(fldnames_fr_glc(nf)), dataptr_out , rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return
            ! Determine export data
            if (ns == 1) then
              dataptr_out(:) = dataptr_in(:)
            else
              dataptr_out(:) = dataptr_out(:) + dataptr_in(:)
            end if
          end if
        end do
      end if
    end do

    ! Check for nans in fields export to rof
    call FB_check_for_nans(is_local%wrap%FBExp(comprof), maintask, logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call fldbun_diagnose(is_local%wrap%FBExp(comprof), &
            string=trim(subname)//' FBexp(comprof) ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! zero accumulator and FBAccum
    !---------------------------------------

    ! zero counter
    lndAccum2rof_cnt = 0

    ! zero lnd2rof fields in FBlndAccum2rof_l
    do n = 1,size(lnd2rof_flds)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldName=trim(lnd2rof_flds(n)), &
            isPresent=exists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          call ESMF_FieldBundleGet(FBlndAccum2rof_l, fieldName=trim(lnd2rof_flds(n)), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(lfield, dataptr_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr_out(:) = czero
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
    type(ESMF_Field)          :: field_lfrac_lnd
    type(ESMF_Mesh)           :: lmesh_lnd
    type(ESMF_Mesh)           :: lmesh_rof
    real(r8), pointer         :: volr_l(:)
    real(r8), pointer         :: volr_r(:)
    real(r8), pointer         :: volr_r_import(:)
    real(r8), pointer         :: irrig_normalized_l(:)
    real(r8), pointer         :: irrig_normalized_r(:)
    real(r8), pointer         :: irrig_volr0_l(:)
    real(r8), pointer         :: irrig_volr0_r(:)
    real(r8), pointer         :: irrig_flux_l(:)
    real(r8), pointer         :: irrig_flux_r(:)
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
    call fldbun_getdata1d(FBlndAccum2rof_l, trim(irrig_flux_field), irrig_flux_l, rc)
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

    ! Convert to a total irrigation flux on the ROF grid, and put this in the pre-merge FBlndAccum2rof_r
    call field_getdata1d(field_rofIrrig, irrig_normalized_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(field_rofIrrig0, irrig_volr0_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(FBlndAccum2rof_r, trim(irrig_flux_field), irrig_flux_r, rc)
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
