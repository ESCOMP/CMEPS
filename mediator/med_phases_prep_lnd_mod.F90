module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use ESMF                  , only : operator(/=)
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use esmFlds               , only : complnd, compatm, comprof, ncomps, compname
  use esmFlds               , only : fldListTo
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_fldchk    => med_methods_FB_fldchk
  use med_methods_mod       , only : State_GetScalar  => med_methods_State_GetScalar
  use med_methods_mod       , only : State_SetScalar  => med_methods_State_SetScalar
  use med_utils_mod         , only : chkerr           => med_utils_ChkErr
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_map_mod           , only : med_map_field_packed, med_map_field_normalized, med_map_field
  use med_merge_mod         , only : med_merge_auto
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_lnd

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_lnd(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    type(InternalState)       :: is_local
    integer                   :: n1,ncnt,ns
    real(r8)                  :: nextsw_cday
    logical                   :: first_call = .true.
    logical                   :: isPresent
    character(CL)             :: cvalue
    character(len=*), parameter :: subname='(med_phases_prep_lnd)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Count the number of fields outside of scalar data, if zero, then return
    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldCount=ncnt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       ! map to create FBimp(:,complnd) - other than glc->lnd
       !---------------------------------------

       call t_startf('MED:'//trim(subname)//' map')
       do n1 = 1,ncomps
          ! Skip glc here and handle it below
          if (n1 == compatm .or. n1 == comprof) then
             if (is_local%wrap%med_coupling_active(n1,complnd)) then
                call med_map_field_packed( &
                     FBSrc=is_local%wrap%FBImp(n1,n1), &
                     FBDst=is_local%wrap%FBImp(n1,complnd), &
                     FBFracSrc=is_local%wrap%FBFrac(n1), &
                     field_normOne=is_local%wrap%field_normOne(n1,complnd,:), &
                     packed_data=is_local%wrap%packed_data(n1,complnd,:), &
                     routehandles=is_local%wrap%RH(n1,complnd,:), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
       end do
       call t_stopf('MED:'//trim(subname)//' map')

       !---------------------------------------
       ! auto merges to create FBExp(complnd) - other than glc->lnd
       !---------------------------------------

       ! The following will merge all fields in fldsSrc

       call t_startf('MED:'//trim(subname)//' merge')
       call med_merge_auto(complnd, &
            is_local%wrap%med_coupling_active(:,complnd), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' merge')

       !---------------------------------------
       ! update scalar data
       !---------------------------------------

       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), trim(is_local%wrap%flds_scalar_name), itemType, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (itemType /= ESMF_STATEITEM_NOTFOUND) then
          call t_startf('MED:'//trim(subname)//' nextsw_cday')
          ! send nextsw_cday to land - first obtain it from atm import
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
               state=is_local%wrap%NstateExp(complnd), &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call t_stopf('MED:'//trim(subname)//' nextsw_cday')
       end if

       ! diagnose
       if (dbug_flag > 1) then
          call fldbun_diagnose(is_local%wrap%FBExp(complnd), &
               string=trim(subname)//' FBexp(complnd) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Set first call logical to false
    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_lnd

end module med_phases_prep_lnd_mod
