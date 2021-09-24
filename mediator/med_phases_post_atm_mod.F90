module med_phases_post_atm_mod

  !-----------------------------------------------------------------------------
  ! Mediator phase for post atm calculations, maps atm->ice, atm->lnd and atm->ocn
  !-----------------------------------------------------------------------------

  implicit none
  private

  public :: med_phases_post_atm

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_atm(gcomp, rc)

    !---------------------------------------
    ! map atm to ocn and atm to ice and atm to land
    !---------------------------------------

    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_phases_history_mod, only : med_phases_history_write_atm
    use med_map_mod           , only : med_map_field_packed
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr    => med_utils_ChkErr
    use esmFlds               , only : compocn, compatm, compice, complnd
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*), parameter :: subname='(med_phases_post_atm)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map atm to ocn
    if (is_local%wrap%med_coupling_active(compatm,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compocn,:), &
            packed_data=is_local%wrap%packed_data(compatm,compocn,:), &
            routehandles=is_local%wrap%RH(compatm,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ocn')
    end if
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
    ! map atm->lnd
    if (is_local%wrap%med_coupling_active(compatm,complnd)) then
       call t_startf('MED:'//trim(subname)//' map_atm2lnd')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,complnd), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,complnd,:), &
            packed_data=is_local%wrap%packed_data(compatm,complnd,:), &
            routehandles=is_local%wrap%RH(compatm,complnd,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2lnd')
    end if

    ! Write atm inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_atm(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_atm

end module med_phases_post_atm_mod
