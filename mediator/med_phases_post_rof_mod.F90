module med_phases_post_rof_mod

  ! Post rof phase, if appropriate, map initial rof->lnd, rof->ocn, rof->ice

  use NUOPC_Mediator        , only : NUOPC_MediatorGet
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_Mesh, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_Field, ESMF_FieldCreate
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleCreate
  use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_VM, ESMF_VMAllreduce, ESMF_REDUCE_SUM
  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : complnd, compocn, compice, comprof
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_utils_mod         , only : chkerr    => med_utils_ChkErr
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_phases_history_mod, only : med_phases_history_write_comp
  use med_map_mod           , only : med_map_field_packed
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use perf_mod              , only : t_startf, t_stopf
  use shr_sys_mod           , only : shr_sys_abort

  implicit none
  private

  public  :: med_phases_post_rof_init
  public  :: med_phases_post_rof
  private :: med_phases_post_rof_create_rof_field_bundle
  private :: med_phases_post_rof_remove_negative_runoff

  ! A local FieldBundle to hold a copy of rof fields, so that when we modify them, we
  ! aren't modifying the import fields in-place.
  type(ESMF_FieldBundle) :: FBrof_r
  integer :: num_rof_fields
  character(len=CS), allocatable :: rof_field_names(:)

  logical :: remove_negative_runoff

  character(len=13), parameter :: fields_to_remove_negative_runoff(4) = &
       ['Forr_rofl    ', &
        'Forr_rofi    ', &
        'Forr_rofl_glc', &
        'Forr_rofi_glc']

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_post_rof_init(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(CL) :: cvalue
    logical       :: isPresent, isSet
    logical       :: flds_wiso

    character(len=*), parameter :: subname='(med_phases_post_rof_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    call med_phases_post_rof_create_rof_field_bundle(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='remove_negative_runoff', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) remove_negative_runoff
    else
      remove_negative_runoff = .false.
    end if

    ! remove_negative_runoff isn't yet set up to handle isotope fields, so ensure that
    ! this isn't set along with flds_wiso
    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_wiso
    else
       flds_wiso = .false.
    end if
    if (remove_negative_runoff .and. flds_wiso) then
      call shr_sys_abort('remove_negative_runoff must be set to false when flds_wiso is true')
    end if

    if (maintask) then
      write(logunit,'(a,l7)') trim(subname)//' remove_negative_runoff = ', remove_negative_runoff
    end if

    if (dbug_flag > 20) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)
  end subroutine med_phases_post_rof_init

  subroutine med_phases_post_rof(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    real(r8), pointer   :: data_orig(:)
    real(r8), pointer   :: data_copy(:)
    integer             :: n
    logical             :: exists
    character(len=*), parameter :: subname='(med_phases_post_rof)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, num_rof_fields
      call fldbun_getdata1d(is_local%wrap%FBImp(comprof,comprof), trim(rof_field_names(n)), data_orig, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call fldbun_getdata1d(FBrof_r, trim(rof_field_names(n)), data_copy, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      data_copy(:) = data_orig(:)
    end do

    if (remove_negative_runoff) then
      do n = 1, size(fields_to_remove_negative_runoff)
        call ESMF_FieldBundleGet(FBrof_r, fieldName=trim(fields_to_remove_negative_runoff(n)), isPresent=exists, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        if (exists) then
          call med_phases_post_rof_remove_negative_runoff(gcomp, fields_to_remove_negative_runoff(n), rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
        end if
      end do
    end if

    ! map rof to lnd
    if (is_local%wrap%med_coupling_active(comprof,complnd)) then
       call t_startf('MED:'//trim(subname)//' map_rof2lnd')
       call med_map_field_packed( &
            FBSrc=FBrof_r, &
            FBDst=is_local%wrap%FBImp(comprof,complnd), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,complnd,:), &
            packed_data=is_local%wrap%packed_data(comprof,complnd,:), &
            routehandles=is_local%wrap%RH(comprof,complnd,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2lnd')
    end if
    ! map rof to ocn
    if (is_local%wrap%med_coupling_active(comprof,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_rof2ocn')
       call med_map_field_packed( &
            FBSrc=FBrof_r, &
            FBDst=is_local%wrap%FBImp(comprof,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,compocn,:), &
            packed_data=is_local%wrap%packed_data(comprof,compocn,:), &
            routehandles=is_local%wrap%RH(comprof,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2ocn')
    end if
    ! map rof to ice
    if (is_local%wrap%med_coupling_active(comprof,compice)) then
       call t_startf('MED:'//trim(subname)//' map_rof2ice')
       call med_map_field_packed( &
            FBSrc=FBrof_r, &
            FBDst=is_local%wrap%FBImp(comprof,compice), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,compice,:), &
            packed_data=is_local%wrap%packed_data(comprof,compice,:), &
            routehandles=is_local%wrap%RH(comprof,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2ice')
    end if

    ! Write rof inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_comp(gcomp, comprof, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_rof

  subroutine med_phases_post_rof_create_rof_field_bundle(gcomp, rc)
    !---------------------------------------------------------------
    ! Create FBrof_r

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer :: n
    type(ESMF_Mesh)  :: mesh
    type(ESMF_Field) :: field
    integer, parameter :: dbug_threshold = 20 ! threshold for writing debug information in this subroutine
    character(len=*), parameter :: subname='(med_phases_post_rof_mod: med_phases_post_rof_create_rof_field_bundle)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > dbug_threshold) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbun_getmesh(is_local%wrap%FBImp(comprof,comprof), mesh, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(comprof,comprof), fieldCount=num_rof_fields, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(rof_field_names(num_rof_fields))
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(comprof,comprof), fieldNameList=rof_field_names, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note that, for simplicity, we'll add all rof fields to this local FieldBundle, even
    ! though we only need to modify a subset of the fields.
    FBrof_r = ESMF_FieldBundleCreate(name='FBrof_r', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1, num_rof_fields
      field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=rof_field_names(n), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleAdd(FBrof_r, (/field/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    if (dbug_flag > dbug_threshold) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_rof_create_rof_field_bundle

  subroutine med_phases_post_rof_remove_negative_runoff(gcomp, field_name, rc)
    !---------------------------------------------------------------
    ! For one runoff field, remove negative runoff by downweighting all positive runoff to
    ! spread the negative runoff globally.

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    character(len=*), intent(in) :: field_name  ! name of runoff flux field to process
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    real(r8), pointer   :: runoff_flux(:)  ! temporary 1d pointer
    real(r8), pointer   :: areas(:)
    real(r8)            :: local_positive(1), global_positive(1)
    real(r8)            :: local_negative(1), global_negative(1)
    real(r8)            :: global_sum
    real(r8)            :: multiplier
    real(r8)            :: local_positive_final(1), global_positive_final(1)
    real(r8)            :: local_negative_final(1), global_negative_final(1)
    real(r8)            :: global_sum_final
    integer :: n

    integer, parameter :: dbug_threshold = 20 ! threshold for writing debug information in this subroutine
    character(len=*), parameter :: subname='(med_phases_post_rof_mod: med_phases_post_rof_remove_negative_runoff)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > dbug_threshold) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note that we don't use rof fractions in the global sum. This is consistent with the
    ! global budget calculations in med_diag_mod and is because the rof fractions are 1
    ! everywhere.
    areas => is_local%wrap%mesh_info(comprof)%areas

    call fldbun_getdata1d(FBrof_r, trim(field_name), runoff_flux, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    local_positive(1) = 0.0_r8
    local_negative(1) = 0.0_r8
    do n = 1, size(runoff_flux)
      if (runoff_flux(n) >= 0.0_r8) then
        local_positive(1) = local_positive(1) + areas(n) * runoff_flux(n)
      else
        local_negative(1) = local_negative(1) + areas(n) * runoff_flux(n)
      end if
    end do

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMAllreduce(vm, senddata=local_positive, recvdata=global_positive, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMAllreduce(vm, senddata=local_negative, recvdata=global_negative, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    global_sum = global_positive(1) + global_negative(1)
    if (maintask .and. dbug_flag > dbug_threshold) then
      write(logunit,'(a)') subname//' Before correction: '//trim(field_name)
      write(logunit,'(a,e27.17)') subname//' global_positive = ', global_positive(1)
      write(logunit,'(a,e27.17)') subname//' global_negative = ', global_negative(1)
      write(logunit,'(a,e27.17)') subname//' global_sum      = ', global_sum
    end if

    if (global_sum > 0.0_r8) then
      ! There is enough positive runoff to absorb all of the negative runoff; so set
      ! negative runoff to 0 and downweight positive runoff to conserve.
      multiplier = global_sum/global_positive(1)
      do n = 1, size(runoff_flux)
        if (runoff_flux(n) > 0.0_r8) then
          runoff_flux(n) = runoff_flux(n) * multiplier
        else
          runoff_flux(n) = 0.0_r8
        end if
      end do
    else if (global_sum < 0.0_r8) then
      ! There is more negative than positive runoff. Hopefully this happens rarely, if
      ! ever; so set positive runoff to 0 and downweight negative runoff to minimize
      ! negative runoff and conserve.
      multiplier = global_sum/global_negative(1)
      do n = 1, size(runoff_flux)
        if (runoff_flux(n) < 0.0_r8) then
          runoff_flux(n) = runoff_flux(n) * multiplier
        else
          runoff_flux(n) = 0.0_r8
        end if
      end do
    else
      ! global_sum == 0 - i.e., positive and negative exactly balance (very rare, unless
      ! the fluxes are already 0 everywhere!); set all fluxes to 0 in this case.
      do n = 1, size(runoff_flux)
        runoff_flux(n) = 0.0_r8
      end do
    end if

    if (dbug_flag > dbug_threshold) then
      ! Recompute positives, negatives and total sum for output diagnostic purposes
      local_positive_final(1) = 0.0_r8
      local_negative_final(1) = 0.0_r8
      do n = 1, size(runoff_flux)
        if (runoff_flux(n) >= 0.0_r8) then
          local_positive_final(1) = local_positive_final(1) + areas(n) * runoff_flux(n)
        else
          local_negative_final(1) = local_negative_final(1) + areas(n) * runoff_flux(n)
        end if
      end do
      call ESMF_VMAllreduce(vm, senddata=local_positive_final, recvdata=global_positive_final, count=1, &
           reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_VMAllreduce(vm, senddata=local_negative_final, recvdata=global_negative_final, count=1, &
           reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      global_sum_final = global_positive_final(1) + global_negative_final(1)
      if (maintask) then
        write(logunit,'(a)') subname//' After correction: '//trim(field_name)
        write(logunit,'(a,e27.17)') subname//' global_positive_final = ', global_positive_final(1)
        write(logunit,'(a,e27.17)') subname//' global_negative_final = ', global_negative_final(1)
        write(logunit,'(a,e27.17)') subname//' global_sum_final      = ', global_sum_final
      end if
    end if

    if (dbug_flag > dbug_threshold) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_rof_remove_negative_runoff

end module med_phases_post_rof_mod
