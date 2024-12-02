module med_phases_post_ocn_mod

  !-----------------------------------------------------------------------------
  ! Mediator post ocn phase - maps ocn->ice, accumulate glc input from ocn
  !-----------------------------------------------------------------------------

  implicit none
  private

  public  :: med_phases_post_ocn

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_ocn(gcomp, rc)

    use NUOPC_Mediator          , only : NUOPC_MediatorGet
    use ESMF                    , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                    , only : ESMF_GridComp
    use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_kind_mod            , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod           , only : chkerr      => med_utils_ChkErr
    use med_constants_mod       , only : dbug_flag   => med_constants_dbug_flag
    use med_map_mod             , only : med_map_field_packed
    use med_internalstate_mod   , only : InternalState
    use med_internalstate_mod   , only : compice, compocn, compwav
    use med_phases_history_mod  , only : med_phases_history_write_comp
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_accum_ocn
    use perf_mod                , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*),parameter :: subname='(med_phases_post_ocn)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Map ocn->ice
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
    ! Map ocn->wav
    if (is_local%wrap%med_coupling_active(compocn,compwav)) then
       call t_startf('MED:'//trim(subname)//' map_ocn2wav')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compocn,compocn), &
            FBDst=is_local%wrap%FBImp(compocn,compwav), &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            field_normOne=is_local%wrap%field_normOne(compocn,compwav,:), &
            packed_data=is_local%wrap%packed_data(compocn,compwav,:), &
            routehandles=is_local%wrap%RH(compocn,compwav,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_ocn2wav')
    end if

    ! Accumulate ocn input for glc if there is ocn->glc coupling
    if (is_local%wrap%ocn2glc_coupling) then
       call ESMF_LogWrite(subname//' DEBUG: calling med_phases_prep_glc_accum_ocn', ESMF_LOGMSG_INFO)
       call med_phases_prep_glc_accum_ocn(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Write ocn inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_comp(gcomp, compocn, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_ocn

end module med_phases_post_ocn_mod
