module med_phases_post_lnd_mod

  implicit none
  private

  public :: med_phases_post_lnd_init ! does not accumulate input to rof
  public :: med_phases_post_lnd

  logical :: lnd2glc_coupling

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_lnd(gcomp, rc)

    use med_kind_mod            , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                    , only : ESMF_GridComp
    use med_constants_mod       , only : dbug_flag   => med_constants_dbug_flag
    use med_utils_mod           , only : chkerr      => med_utils_ChkErr
    use med_methods_mod         , only : FB_diagnose => med_methods_FB_diagnose
    use med_map_mod             , only : med_map_field_packed
    use med_internalstate_mod   , only : InternalState, mastertask
    use med_phases_prep_rof_mod , only : med_phases_prep_rof_accum
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_accum
    use esmFlds                 , only : complnd, compatm, comprof, compglc, num_icesheets
    use perf_mod                , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: ns
    logical             :: first_call = .true.
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

    ! first determine if there will be any lnd to glc coupling
    if (first_call) then
       do ns = 1,num_icesheets
          if (is_local%wrap%med_coupling_active(complnd,compglc(ns))) then
             lnd2glc_coupling = .true.
             exit
          end if
       end do
       first_call = .false.
    end if

    ! accumulate lnd input for glc
    if (lnd2glc_coupling) then
       call med_phases_prep_glc_accum(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_lnd

  !===============================================================================
  subroutine med_phases_post_lnd_init(gcomp, rc)

    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridComp
    use med_constants_mod     , only : dbug_flag   => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr      => med_utils_ChkErr
    use med_methods_mod       , only : FB_diagnose => med_methods_FB_diagnose
    use med_map_mod           , only : med_map_field_packed
    use med_internalstate_mod , only : InternalState, mastertask
    use esmFlds               , only : complnd, compatm
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
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

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_lnd_init

end module med_phases_post_lnd_mod
