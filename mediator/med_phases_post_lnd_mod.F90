module med_phases_post_lnd_mod

  implicit none
  private

  public :: med_phases_post_lnd

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_lnd(gcomp, rc)

    use NUOPC_Mediator          , only : NUOPC_MediatorGet
    use ESMF                    , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                    , only : ESMF_GridComp
    use med_kind_mod            , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_constants_mod       , only : dbug_flag   => med_constants_dbug_flag
    use med_utils_mod           , only : chkerr      => med_utils_ChkErr
    use med_methods_mod         , only : FB_diagnose => med_methods_FB_diagnose
    use med_map_mod             , only : med_map_field_packed
    use med_internalstate_mod   , only : InternalState, mastertask
    use med_phases_prep_rof_mod , only : med_phases_prep_rof_accum
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_accum_lnd, med_phases_prep_glc_avg
    use med_phases_history_mod  , only : med_phases_history_write_comp
    use esmFlds                 , only : complnd, compatm, comprof, compglc, num_icesheets
    use esmFlds                 , only : lnd2glc_coupling, accum_lnd2glc
    use perf_mod                , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*),parameter :: subname='(med_phases_post_lnd)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! If driver clock is created then are in the run phase otherwise are in the initialization phase
    if (ESMF_ClockIsCreated(dclock)) then

       ! map lnd to atm
       if (is_local%wrap%med_coupling_active(complnd,compatm)) then
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBImp(complnd,complnd), &
               FBDst=is_local%wrap%FBImp(complnd,compatm), &
               FBFracSrc=is_local%wrap%FBFrac(complnd), &
               field_NormOne=is_local%wrap%field_normOne(complnd,compatm,:), &
               packed_data=is_local%wrap%packed_data(complnd,compatm,:), &
               routehandles=is_local%wrap%RH(complnd,compatm,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! accumulate lnd input for rof
       if (is_local%wrap%med_coupling_active(complnd,comprof)) then
          call med_phases_prep_rof_accum(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! accumulate lnd input for glc (note that lnd2glc_coupling and accum_lnd2glc is determined in med.F90)
       if (lnd2glc_coupling) then
          call med_phases_prep_glc_accum_lnd(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          ! Note that in this case med_phases_prep_glc_avg is called 
          ! from med_phases_prep_glc in the run sequence
       else if (accum_lnd2glc) then
          call med_phases_prep_glc_accum_lnd(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_phases_prep_glc_avg(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Write lnd inst, avg or aux if requested in mediator attributes
       call med_phases_history_write_comp(gcomp, complnd, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    else

       ! initialization phase - map lnd to atm
       if (is_local%wrap%med_coupling_active(complnd,compatm)) then
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBImp(complnd,complnd), &
               FBDst=is_local%wrap%FBImp(complnd,compatm), &
               FBFracSrc=is_local%wrap%FBFrac(complnd), &
               field_NormOne=is_local%wrap%field_normOne(complnd,compatm,:), &
               packed_data=is_local%wrap%packed_data(complnd,compatm,:), &
               routehandles=is_local%wrap%RH(complnd,compatm,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_lnd

end module med_phases_post_lnd_mod
