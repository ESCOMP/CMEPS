module med_phases_post_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases updating fractions and mapping ice->atm and ice->ocn
  !-----------------------------------------------------------------------------

  implicit none
  private

  public :: med_phases_post_ice

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_ice(gcomp, rc)

    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridComp
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_constants_mod     , only : dbug_flag   => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr      => med_utils_ChkErr
    use med_methods_mod       , only : FB_diagnose => med_methods_FB_diagnose
    use med_map_mod           , only : med_map_field_packed
    use med_fraction_mod      , only : med_fraction_set
    use med_internalstate_mod , only : InternalState, mastertask
    use med_phases_history_mod, only : med_phases_history_write_ice
    use esmFlds               , only : compice, compatm, compocn, compwav
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*),parameter :: subname='(med_phases_post_ice)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! update ice fraction
    call med_fraction_set(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map ice to atm - scaling by updated ice fraction
    if (is_local%wrap%med_coupling_active(compice,compatm)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compice,compice), &
            FBDst=is_local%wrap%FBImp(compice,compatm), &
            FBFracSrc=is_local%wrap%FBFrac(compice), &
            field_NormOne=is_local%wrap%field_normOne(compice,compatm,:), &
            packed_data=is_local%wrap%packed_data(compice,compatm,:), &
            routehandles=is_local%wrap%RH(compice,compatm,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if
    ! map ice to ocn
    if (is_local%wrap%med_coupling_active(compice,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_ice2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compice,compice), &
            FBDst=is_local%wrap%FBImp(compice,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compice), &
            field_normOne=is_local%wrap%field_normOne(compice,compocn,:), &
            packed_data=is_local%wrap%packed_data(compice,compocn,:), &
            routehandles=is_local%wrap%RH(compice,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_ice2ocn')
    end if
    ! map ice to wav
    if (is_local%wrap%med_coupling_active(compice,compwav)) then
       call t_startf('MED:'//trim(subname)//' map_ice2wav')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compice,compice), &
            FBDst=is_local%wrap%FBImp(compice,compwav), &
            FBFracSrc=is_local%wrap%FBFrac(compice), &
            field_normOne=is_local%wrap%field_normOne(compice,compwav,:), &
            packed_data=is_local%wrap%packed_data(compice,compwav,:), &
            routehandles=is_local%wrap%RH(compice,compwav,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_ice2wav')
    end if

    ! Write ice inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_ice(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if

  end subroutine med_phases_post_ice

end module med_phases_post_ice_mod
