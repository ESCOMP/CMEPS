module med_phases_prep_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ice export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_utils_mod         , only : chkerr            => med_utils_ChkErr
  use med_methods_mod       , only : fldchk            => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_diagnose       => med_methods_FB_diagnose
  use med_methods_mod       , only : State_GetScalar   => med_methods_State_GetScalar
  use med_methods_mod       , only : State_SetScalar   => med_methods_State_SetScalar
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_merge_mod         , only : med_merge_auto
  use med_map_mod           , only : med_map_field_packed
  use med_internalstate_mod , only : InternalState, logunit, mastertask
  use esmFlds               , only : compatm, compice, compocn, comprof, compglc, ncomps, compname
  use esmFlds               , only : fldListFr, fldListTo
  use esmFlds               , only : mapnames
  use esmFlds               , only : coupling_mode
  use esmFlds               , only : med_fldList_GetFldInfo
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_ice

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ice(gcomp, rc)

    use ESMF  , only : operator(/=)
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_StateGet 
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_Field
    use ESMF  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use ESMF  , only : ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
    use NUOPC , only : NUOPC_IsConnected

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_StateItem_Flag)      :: itemType
    type(InternalState)            :: is_local
    type(ESMF_Field)               :: lfield
    integer                        :: i,n,n1,ncnt
    character(len=CS)              :: fldname
    integer                        :: fldnum
    integer                        :: mapindex
    real(R8), pointer              :: dataptr(:) => null()
    real(R8)                       :: precip_fact
    character(len=CS)              :: cvalue
    character(len=64), allocatable :: fldnames(:)
    real(r8)                       :: nextsw_cday
    logical                        :: first_precip_fact_call = .true.
    character(len=*),parameter     :: subname='(med_phases_prep_ice)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compice), fieldCount=ncnt, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       ! map atm->ice
       if (is_local%wrap%med_coupling_active(compatm,compice)) then
          call t_startf('MED:'//trim(subname)//' map_atm2ice')
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBImp(compatm,compatm), &
               FBDst=is_local%wrap%FBImp(compatm,compice), &
               FBFracSrc=is_local%wrap%FBFrac(compatm), &
               field_normOne=is_local%wrap%field_normOne(compatm,compice,:), &
               packed_data=is_local%wrap%packed_data(compatm,compice,:), &
               routehandles=is_local%wrap%RH(compatm,compice,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call t_stopf('MED:'//trim(subname)//' map_atm2ice')
       end if
       ! map ocn->ice
       if (is_local%wrap%med_coupling_active(compocn,compice)) then
          call t_startf('MED:'//trim(subname)//' map_ocn2ice')
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compocn,compice), &
               FBFracSrc=is_local%wrap%FBFrac(compocn), &
               field_normOne=is_local%wrap%field_normOne(compocn,compice,:), &
               packed_data=is_local%wrap%packed_data(compocn,compice,:), &
               routehandles=is_local%wrap%RH(compocn,compice,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call t_stopf('MED:'//trim(subname)//' map_ocn2ice')
       end if
       ! rof->ice is mapped in med_phases_post_rof
       ! glc->ice is mapped in med_phases_post_glc
       
       !---------------------------------------
       !--- auto merges to create FBExp(compice)
       !---------------------------------------

       call med_merge_auto(compice, &
            is_local%wrap%med_coupling_active(:,compice), &
            is_local%wrap%FBExp(compice), &
            is_local%wrap%FBFrac(compice), &
            is_local%wrap%FBImp(:,compice), &
            fldListTo(compice), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       if (trim(coupling_mode) == 'cesm') then

          ! application of precipitation factor from ocean
          ! TODO (mvertens, 2019-03-18): precip_fact here is not valid if
          ! the component does not send it - hardwire it to 1 until this is resolved
          precip_fact = 1.0_R8
          if (precip_fact /= 1.0_R8) then
             if (first_precip_fact_call .and. mastertask) then
                write(logunit,'(a)')'(merge_to_ice): Scaling rain, snow, liquid and ice runoff by precip_fact '
                first_precip_fact_call = .false.
             end if
             write(cvalue,*) precip_fact
             call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO)

             allocate(fldnames(3))
             fldnames = (/'Faxa_rain', 'Faxa_snow', 'Fixx_rofi'/)
             do n = 1,size(fldnames)
                if (fldchk(is_local%wrap%FBExp(compice), trim(fldnames(n)), rc=rc)) then
                   call ESMF_FieldBundleGet(is_local%wrap%FBExp(compice), fieldname=trim(fldnames(n)), &
                        field=lfield, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   call ESMF_FieldGet(lfield, farrayptr=dataptr, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   dataptr(:) = dataptr(:) * precip_fact
                end if
             end do
             deallocate(fldnames)
          end if
       end if

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compice), string=trim(subname)//' FBexp(compice) ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- update scalar data
       !---------------------------------------

       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), trim(is_local%wrap%flds_scalar_name), itemType, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (itemType /= ESMF_STATEITEM_NOTFOUND) then
          if (is_local%wrap%flds_scalar_index_nextsw_cday .ne. 0) then
             ! send nextsw_cday to ice - first obtain it from atm import 
             call State_GetScalar(&
                  scalar_value=nextsw_cday, &
                  scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
                  state=is_local%wrap%NstateImp(compatm), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call State_SetScalar(&
                  scalar_value=nextsw_cday, &
                  scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
                  state=is_local%wrap%NstateExp(compice), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ice

end module med_phases_prep_ice_mod
