module med_phases_prep_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ice export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  public  :: med_phases_prep_ice

  real(r8), pointer :: dataptr_scalar_ice(:,:) => null()
  real(r8), pointer :: dataptr_scalar_atm(:,:) => null()

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ice(gcomp, rc)

    use ESMF                  , only : operator(/=)
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_StateGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_Field
    use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use ESMF                  , only : ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
    use ESMF                  , only : ESMF_VMBroadCast
    use med_utils_mod         , only : chkerr       => med_utils_ChkErr
    use med_methods_mod       , only : FB_fldchk    => med_methods_FB_FldChk
    use med_methods_mod       , only : FB_diagnose  => med_methods_FB_diagnose
    use med_methods_mod       , only : FB_GetFldPtr => med_methods_FB_GetFldPtr
    use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
    use med_merge_mod         , only : med_merge_auto
    use med_internalstate_mod , only : InternalState, logunit, mastertask
    use esmFlds               , only : compatm, compice, compocn, comprof, compglc, ncomps, compname
    use esmFlds               , only : fldListTo
    use esmFlds               , only : coupling_mode
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)            :: is_local
    type(ESMF_Field)               :: lfield
    integer                        :: i,n
    real(R8), pointer              :: dataptr(:)
    real(R8), pointer              :: dataptr_scalar_ocn(:,:)
    real(R8)                       :: precip_fact(1)
    character(len=CS)              :: cvalue
    character(len=64), allocatable :: fldnames(:)
    real(r8)                       :: nextsw_cday
    integer                        :: scalar_id
    real(r8)                       :: tmp(1)
    logical                        :: first_precip_fact_call = .true.
    character(len=*),parameter     :: subname='(med_phases_prep_ice)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! atm->ice is mapped in med_phases_post_atm
    ! glc->ice is mapped in med_phases_post_glc
    ! rof->ice is mapped in med_phases_post_rof
    ! ocn->ice is mapped in med_phases_post_ocn

    ! auto merges to create FBExp(compice)
    call med_merge_auto(&
         is_local%wrap%med_coupling_active(:,compice), &
         is_local%wrap%FBExp(compice), &
         is_local%wrap%FBFrac(compice), &
         is_local%wrap%FBImp(:,compice), &
         fldListTo(compice), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Apply precipitation factor from ocean (that scales atm rain and snow to ice) if appropriate
    if (trim(coupling_mode) == 'cesm' .and. is_local%wrap%flds_scalar_index_precip_factor /= 0) then

       ! Note that in med_internal_mod.F90 all is_local%wrap%flds_scalar_index_precip_factor
       ! is initialized to 0.
       ! In addition, in med.F90, if this attribute is not present as a mediator component attribute,
       ! it is set to 0.
       if (mastertask) then
          call ESMF_StateGet(is_local%wrap%NstateImp(compocn), &
               itemName=trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_ocn, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          scalar_id=is_local%wrap%flds_scalar_index_precip_factor
          precip_fact(1) = dataptr_scalar_ocn(scalar_id,1)
          if (precip_fact(1) /= 1._r8) then
             write(logunit,'(a,f21.13)')&
                  '(merge_to_ice): Scaling rain, snow, liquid and ice runoff by non-unity precip_fact ',&
                  precip_fact(1)
          end if
       end if
       call ESMF_VMBroadCast(is_local%wrap%vm, precip_fact, 1, 0, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       is_local%wrap%flds_scalar_precip_factor = precip_fact(1)
       if (dbug_flag > 5) then
          write(cvalue,*) precip_fact(1)
          call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO)
       end if

       ! Scale rain and snow to ice from atm by the precipitation factor received from the ocean
       allocate(fldnames(3))
       fldnames = (/'Faxa_rain', 'Faxa_snow', 'Fixx_rofi'/)
       do n = 1,size(fldnames)
          if (FB_fldchk(is_local%wrap%FBExp(compice), trim(fldnames(n)), rc=rc)) then
             call FB_GetFldPtr(is_local%wrap%FBExp(compice), trim(fldnames(n)), dataptr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             dataptr(:) = dataptr(:) * is_local%wrap%flds_scalar_precip_factor
          end if
       end do
       deallocate(fldnames)
    end if

    ! obtain nextsw_cday from atm if it is in the import state and send it to ice
    scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday
    if (scalar_id > 0 .and. mastertask) then
       call ESMF_StateGet(is_local%wrap%NstateImp(compatm), &
            itemName=trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_atm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_StateGet(is_local%wrap%NStateExp(compice), &
            trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_ice, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       dataptr_scalar_ice(scalar_id,1) = dataptr_scalar_atm(scalar_id,1)
    end if

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBExp(compice), string=trim(subname)//' FBexp(compice) ', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ice

end module med_phases_prep_ice_mod
