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
  use med_methods_mod       , only : FB_getfldptr=> med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_check_for_nans => med_methods_FB_check_for_nans
  use med_merge_mod         , only : med_merge_auto
  use med_map_mod           , only : med_map_field_packed
  use med_internalstate_mod , only : InternalState, maintask, logunit, samegrid_atmlnd
  use med_internalstate_mod , only : compatm, compocn, compice, compname, coupling_mode
  use esmFlds               , only : med_fldlist_GetfldListTo, med_fldlist_type
  use perf_mod              , only : t_startf, t_stopf
  use med_phases_aofluxes_mod, only : med_aofluxes_map_xgrid2agrid_output
  use med_phases_aofluxes_mod, only : med_aofluxes_map_ogrid2agrid_output

  implicit none
  private

  public :: med_phases_prep_atm
  public :: med_phases_prep_atm_enthalpy_correction

  real(r8), public :: global_htot_corr(1) = 0._r8  ! enthalpy correction from med_phases_prep_ocn

  character(len=13) :: fldnames_from_ocn(5) = (/'Faoo_fbrf_ocn','Faoo_fdms_ocn','Faoo_fco2_ocn',&
                                                'Faoo_fn2o_ocn','Faoo_fnh3_ocn'/)

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
    type(InternalState)        :: is_local
    real(R8), pointer          :: dataPtr1(:)
    real(R8), pointer          :: dataPtr2(:)
    real(R8), pointer          :: ifrac(:)
    real(R8), pointer          :: ofrac(:)
    integer                    :: n,nf
    type(med_fldlist_type), pointer :: fldList
    character(len=*),parameter :: subname='(med_phases_prep_atm)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 3, maintask)

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
            FBDat=is_local%wrap%FBData(compatm), &
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
    if (trim(coupling_mode) == 'cesm' .or. &
         trim(coupling_mode) == 'ufs.frac.aoflux') then
       if (is_local%wrap%aoflux_grid == 'ogrid') then
          call med_aofluxes_map_ogrid2agrid_output(gcomp, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else if (is_local%wrap%aoflux_grid == 'agrid') then
          ! Do nothing - fluxes are alread being computed on the agrid
       else if (is_local%wrap%aoflux_grid == 'xgrid') then
          call med_aofluxes_map_xgrid2agrid_output(gcomp, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    !---------------------------------------
    !--- merge all fields to atm
    !---------------------------------------
    fldList => med_fldList_GetfldListTo(compatm)
    call med_merge_auto(&
         is_local%wrap%med_coupling_active(:,compatm), &
         is_local%wrap%FBExp(compatm), &
         is_local%wrap%FBFrac(compatm), &
         is_local%wrap%FBImp(:,compatm), &
         fldList, &
         FBMed1=is_local%wrap%FBMed_ocnalb_a, &
         FBMed2=is_local%wrap%FBMed_aoflux_a, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
       if (samegrid_atmlnd) then
          call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='lfrac', field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='lfrin', field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr1)
          dataptr1(n) = dataptr2(n)
       end do
    end if

    ! Note - the following needs a custom merge since Faoo_fco2_ocn is scaled by (ifrac+ofrac)
    ! in the merge to the atm
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ifrac', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=ifrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ofrac', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=ofrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_from_ocn)
      if ( FB_FldChk(is_local%wrap%FBExp(compatm)        , trim(fldnames_from_ocn(nf)), rc=rc) .and. &
           FB_FldChk(is_local%wrap%FBImp(compocn,compocn), trim(fldnames_from_ocn(nf)), rc=rc)) then
        call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compatm), &
             fieldName=trim(fldnames_from_ocn(nf)), field=lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_FieldGet(lfield, farrayPtr=dataptr1, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_FieldBundleGet(is_local%wrap%FBExp(compatm), &
             fieldName=trim(fldnames_from_ocn(nf)), field=lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_FieldGet(lfield, farrayPtr=dataptr2, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        do n = 1,size(dataptr2)
          dataptr2(n) = (ifrac(n) + ofrac(n)) * dataptr1(n)
        end do
      end if
    end do

    ! Add enthalpy correction to sensible heat if appropriate
    if (FB_FldChk(is_local%wrap%FBExp(compatm), 'Faxx_sen', rc=rc)) then
       call FB_getfldptr(is_local%wrap%FBExp(compatm), 'Faxx_sen', dataptr1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(dataptr1)
          dataptr1(n) = dataptr1(n) + global_htot_corr(1)
       end do
    end if

    ! Check for nans in fields export to atm
    call FB_check_for_nans(is_local%wrap%FBExp(compatm), maintask, logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_atm

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_atm_enthalpy_correction (gcomp, hcorr, rc)

    ! Enthalpy correction term calculation called by med_phases_prep_ocn_accum in
    ! med_phases_prep_ocn_mod
    ! Note that this is only called if the following fields are in FBExp(compocn)
    ! 'Faxa_rain','Foxx_hrain','Faxa_snow' ,'Foxx_hsnow',
    ! 'Foxx_evap','Foxx_hevap','Foxx_hcond','Foxx_rofl',
    ! 'Foxx_hrofl','Foxx_rofi','Foxx_hrofi'

    use ESMF            , only : ESMF_VMAllreduce, ESMF_GridCompGet, ESMF_REDUCE_SUM
    use ESMF            , only : ESMF_VM

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    real(r8)            , intent(in)  :: hcorr(:)
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n
    real(r8)            :: local_htot_corr(1)
    type(ESMF_VM)       :: vm
    !---------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine sum of enthalpy correction for each hcorr index locally
    local_htot_corr(1) = 0._r8
    do n = 1,size(hcorr)
       local_htot_corr(1) = local_htot_corr(1) + hcorr(n)
    end do
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMAllreduce(vm, senddata=local_htot_corr, recvdata=global_htot_corr, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_prep_atm_enthalpy_correction

end module med_phases_prep_atm_mod
