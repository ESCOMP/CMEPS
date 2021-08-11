module med_phases_post_wav_mod

  implicit none
  private

  public :: med_phases_post_wav

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_wav(gcomp, rc)

    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridComp
    use med_constants_mod     , only : dbug_flag   => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr      => med_utils_ChkErr
    use med_methods_mod       , only : FB_diagnose => med_methods_FB_diagnose
    use med_map_mod           , only : med_map_field_packed
    use med_internalstate_mod , only : InternalState, mastertask
    use esmFlds               , only : compwav, compatm, compocn, compice
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    character(len=*),parameter :: subname='(med_phases_post_wav)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map wav to atm
    if (is_local%wrap%med_coupling_active(compwav,compatm)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compwav,compwav), &
            FBDst=is_local%wrap%FBImp(compwav,compatm), &
            FBFracSrc=is_local%wrap%FBFrac(compwav), &
            field_NormOne=is_local%wrap%field_normOne(compwav,compatm,:), &
            packed_data=is_local%wrap%packed_data(compwav,compatm,:), &
            routehandles=is_local%wrap%RH(compwav,compatm,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    ! map wav to ocn
    if (is_local%wrap%med_coupling_active(compwav,compocn)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compwav,compwav), &
            FBDst=is_local%wrap%FBImp(compwav,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compwav), &
            field_NormOne=is_local%wrap%field_normOne(compwav,compocn,:), &
            packed_data=is_local%wrap%packed_data(compwav,compocn,:), &
            routehandles=is_local%wrap%RH(compwav,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    ! map wav to ice
    if (is_local%wrap%med_coupling_active(compwav,compice)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compwav,compwav), &
            FBDst=is_local%wrap%FBImp(compwav,compice), &
            FBFracSrc=is_local%wrap%FBFrac(compwav), &
            field_NormOne=is_local%wrap%field_normOne(compwav,compice,:), &
            packed_data=is_local%wrap%packed_data(compwav,compice,:), &
            routehandles=is_local%wrap%RH(compwav,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if

  end subroutine med_phases_post_wav

end module med_phases_post_wav_mod
