module med_phases_post_rof_mod

  ! Post rof phase, if appropriate, map initial rof->lnd, rof->ocn, rof->ice

  implicit none
  private

  public :: med_phases_post_rof

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_post_rof(gcomp, rc)

    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use esmFlds               , only : complnd, compocn, compice, compatm, comprof, ncomps, compname
    use med_utils_mod         , only : chkerr    => med_utils_ChkErr
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_map_mod           , only : med_map_field_packed
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    character(len=*), parameter :: subname='(med_phases_post_rof)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map rof to lnd
    if (is_local%wrap%med_coupling_active(comprof,complnd)) then
       call t_startf('MED:'//trim(subname)//' map_rof2lnd')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(comprof,comprof), &
            FBDst=is_local%wrap%FBImp(comprof,complnd), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,complnd,:), &
            packed_data=is_local%wrap%packed_data(comprof,complnd,:), &
            routehandles=is_local%wrap%RH(comprof,complnd,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2lnd')
    end if
    ! map rof to ocn
    if (is_local%wrap%med_coupling_active(comprof,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_rof2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(comprof,comprof), &
            FBDst=is_local%wrap%FBImp(comprof,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,compocn,:), &
            packed_data=is_local%wrap%packed_data(comprof,compocn,:), &
            routehandles=is_local%wrap%RH(comprof,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2ocn')
    end if
    ! map rof to ice
    if (is_local%wrap%med_coupling_active(comprof,compice)) then
       call t_startf('MED:'//trim(subname)//' map_rof2ice')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(comprof,comprof), &
            FBDst=is_local%wrap%FBImp(comprof,compice), &
            FBFracSrc=is_local%wrap%FBFrac(comprof), &
            field_normOne=is_local%wrap%field_normOne(comprof,compice,:), &
            packed_data=is_local%wrap%packed_data(comprof,compice,:), &
            routehandles=is_local%wrap%RH(comprof,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_rof2ice')
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_rof

end module med_phases_post_rof_mod
