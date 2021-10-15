module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  public  :: med_phases_prep_lnd

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_lnd(gcomp, rc)

    use NUOPC                 , only : NUOPC_CompAttributeGet
    use ESMF                  , only : operator(/=)
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
    use esmFlds               , only : complnd, compatm, ncomps
    use esmFlds               , only : fldListTo
    use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
    use med_utils_mod         , only : chkerr           => med_utils_ChkErr
    use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_merge_mod         , only : med_merge_auto
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemType
    type(InternalState)         :: is_local
    type(ESMF_Field)            :: lfield
    integer                     :: ncnt,ns
    real(r8)                    :: nextsw_cday
    integer                     :: scalar_id
    real(r8)                    :: tmp(1)
    real(r8), pointer           :: dataptr2d(:,:)
    logical                     :: first_call = .true.
    real(r8), pointer           :: dataptr_scalar_lnd(:,:)
    real(r8), pointer           :: dataptr_scalar_atm(:,:)
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

       ! atm2lnd is done in med_phases_post_atm
       ! rof2lnd is done in med_phases_post_rof
       ! glc2lnd is done in med_phases_post_glc

       ! auto merges to create FBExp(complnd) - other than glc->lnd
       ! The following will merge all fields in fldsSrc
       call t_startf('MED:'//trim(subname)//' merge')
       call med_merge_auto(&
            is_local%wrap%med_coupling_active(:,complnd), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' merge')

       ! obtain nextsw_cday from atm if it is in the import state and send it to lnd
       scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday
       if (scalar_id > 0 .and. mastertask) then
          call ESMF_StateGet(is_local%wrap%NstateImp(compatm), &
               itemName=trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_atm, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_StateGet(is_local%wrap%NStateExp(complnd), &
               trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_lnd, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday
          dataptr_scalar_lnd(scalar_id,1) = dataptr_scalar_atm(scalar_id,1)
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
