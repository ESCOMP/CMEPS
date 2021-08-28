module med_phases_prep_atm_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing atm export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use med_constants_mod     , only : dbug_flag   => med_constants_dbug_flag
  use med_utils_mod         , only : memcheck    => med_memcheck
  use med_utils_mod         , only : chkerr      => med_utils_ChkErr
  use med_methods_mod       , only : FB_diagnose => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_fldchk   => med_methods_FB_FldChk
  use med_merge_mod         , only : med_merge_auto
  use med_map_mod           , only : med_map_field_packed
  use med_internalstate_mod , only : InternalState, mastertask
  use esmFlds               , only : compatm, compocn, compice, ncomps, compname
  use esmFlds               , only : fldListTo, fldListMed_aoflux, coupling_mode
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_prep_atm

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_atm(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Field)           :: lfield
    character(len=64)          :: timestr
    type(InternalState)        :: is_local
    real(R8), pointer          :: dataPtr1(:) => null()
    real(R8), pointer          :: dataPtr2(:) => null()
    real(R8), pointer          :: ifrac(:) => null()
    real(R8), pointer          :: ofrac(:) => null()
    integer                    :: i, j, n, n1, ncnt
    character(len=*),parameter :: subname='(med_phases_prep_atm)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 3, mastertask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- map ocn and ice to atm
    !---------------------------------------
    if (is_local%wrap%med_coupling_active(compocn,compatm)) then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compocn,compocn), &
            FBDst=is_local%wrap%FBImp(compocn,compatm), &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            field_NormOne=is_local%wrap%field_normOne(compocn,compatm,:), &
            packed_data=is_local%wrap%packed_data(compocn,compatm,:), &
            routehandles=is_local%wrap%RH(compocn,compatm,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
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

    !---------------------------------------
    !--- map ocean albedos from ocn to atm grid if appropriate
    !---------------------------------------
    if (trim(coupling_mode) == 'cesm') then
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBMed_ocnalb_o, &
            FBDst=is_local%wrap%FBMed_ocnalb_a, &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            field_normOne=is_local%wrap%field_normOne(compocn,compatm,:), &
            packed_data=is_local%wrap%packed_data_ocnalb_o2a(:), &
            routehandles=is_local%wrap%RH(compocn,compatm,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    !--- map atm/ocn fluxes from ocn to atm grid if appropriate
    !---------------------------------------
    if (trim(coupling_mode) == 'cesm' .or. trim(coupling_mode) == 'hafs') then
       if (is_local%wrap%aoflux_grid == 'ogrid') then
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBMed_aoflux_o, &
               FBDst=is_local%wrap%FBMed_aoflux_a, &
               FBFracSrc=is_local%wrap%FBFrac(compocn), &
               field_normOne=is_local%wrap%field_normOne(compocn,compatm,:), &
               packed_data=is_local%wrap%packed_data_aoflux_o2a(:), &
               routehandles=is_local%wrap%RH(compocn,compatm,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (is_local%wrap%aoflux_grid == 'agrid') then
          ! do nothing - is_local%wrap%FBMed_aoflux_a has been computed in med_aofluxes_init_agrid
       else if (is_local%wrap%aoflux_grid == 'xgrid') then
          ! do nothing - is_local%wrap%FBMed_aoflux_a has been computed in med_aofluxes_init_agrid
       end if
    endif

    !---------------------------------------
    !--- merge all fields to atm
    !---------------------------------------
    if (trim(coupling_mode) == 'cesm' .or. trim(coupling_mode) == 'hafs') then
       call med_merge_auto(compatm, &
            is_local%wrap%med_coupling_active(:,compatm), &
            is_local%wrap%FBExp(compatm), &
            is_local%wrap%FBFrac(compatm), &
            is_local%wrap%FBImp(:,compatm), &
            fldListTo(compatm), &
            FBMed1=is_local%wrap%FBMed_ocnalb_a, &
            FBMed2=is_local%wrap%FBMed_aoflux_a, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (trim(coupling_mode) == 'nems_frac' .or. trim(coupling_mode) == 'nems_orig') then
       call med_merge_auto(compatm, &
            is_local%wrap%med_coupling_active(:,compatm), &
            is_local%wrap%FBExp(compatm), &
            is_local%wrap%FBFrac(compatm), &
            is_local%wrap%FBImp(:,compatm), &
            fldListTo(compatm), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBExp(compatm),string=trim(subname)//' FBexp(compatm) ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    ! set fractions to send back to atm
    if (FB_FldChk(is_local%wrap%FBExp(compatm), 'So_ofrac', rc=rc)) then
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compatm), fieldName='So_ofrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ofrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr1)
          dataptr1(n) = dataptr2(n)
       end do
    end if
    if (FB_FldChk(is_local%wrap%FBExp(compatm), 'Si_ifrac', rc=rc)) then
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compatm), fieldName='Si_ifrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ifrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr1)
          dataptr1(n) = dataptr2(n)
       end do
    end if
    if (FB_FldChk(is_local%wrap%FBExp(compatm), 'Sl_lfrac', rc=rc)) then
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compatm), fieldName='Sl_lfrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='lfrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr1)
          dataptr1(n) = dataptr2(n)
       end do
    end if

    ! Note - the following needs a custom merge since Faoo_fco2_ocn is scaled by (ifrac+ofrac)
    ! in the merge to the atm
    if ( FB_FldChk(is_local%wrap%FBExp(compatm)        , 'Faoo_fco2_ocn', rc=rc) .and. &
         FB_FldChk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fco2_ocn', rc=rc)) then
       call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ifrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=ifrac, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ofrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=ofrac, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compatm), fieldName='Faoo_fco2_ocn', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compatm), fieldName='Faoo_fco2_ocn', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr2)
          dataptr2(n) = (ifrac(n) + ofrac(n)) * dataptr1(n)
       end do
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_atm

end module med_phases_prep_atm_mod
