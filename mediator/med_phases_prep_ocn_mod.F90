module med_phases_prep_ocn_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ocn export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod     , only : czero     =>med_constants_czero
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_merge_mod         , only : med_merge_auto, med_merge_field
  use med_map_mod           , only : med_map_field_packed
  use med_utils_mod         , only : memcheck      => med_memcheck
  use med_utils_mod         , only : chkerr        => med_utils_ChkErr
  use med_methods_mod       , only : FB_diagnose   => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_fldchk     => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_accum      => med_methods_FB_accum
  use med_methods_mod       , only : FB_average    => med_methods_FB_average
  use med_methods_mod       , only : FB_copy       => med_methods_FB_copy
  use med_methods_mod       , only : FB_reset      => med_methods_FB_reset
  use med_methods_mod       , only : FB_check_for_nans => med_methods_FB_check_for_nans
  use esmFlds               , only : med_fldList_GetfldListTo, med_fldlist_type
  use med_internalstate_mod , only : compocn, compatm, compice, coupling_mode
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_prep_ocn_init   ! called from med.F90
  public :: med_phases_prep_ocn_accum  ! called from run sequence
  public :: med_phases_prep_ocn_avg    ! called from run sequence

  private :: med_phases_prep_ocn_custom

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_init(gcomp, rc)

    use ESMF            , only : ESMF_GridComp, ESMF_SUCCESS
    use med_methods_mod , only : FB_Init  => med_methods_FB_init

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    character(len=*),parameter  :: subname=' (med_phases_prep_ocn_init) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    if (maintask) then
       write(logunit,'(a)') trim(subname)//' initializing ocean export accumulation FB for '
    end if
    call FB_init(is_local%wrap%FBExpAccumOcn, is_local%wrap%flds_scalar_name, &
         STgeom=is_local%wrap%NStateExp(compocn), STflds=is_local%wrap%NStateExp(compocn), &
         name='FBExpAccumOcn', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_reset(is_local%wrap%FBExpAccumOcn, value=czero, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_prep_ocn_init

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_accum(gcomp, rc)

    use ESMF                    , only : ESMF_GridComp, ESMF_FieldBundleGet
    use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_constants_mod       , only : shr_const_cpsw, shr_const_tkfrz, shr_const_pi
    use med_phases_prep_atm_mod , only : med_phases_prep_atm_enthalpy_correction

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n
    real(r8)            :: glob_area_inv
    real(r8), pointer   :: tocn(:)
    real(r8), pointer   :: rain(:), hrain(:)
    real(r8), pointer   :: snow(:), hsnow(:)
    real(r8), pointer   :: evap(:), hevap(:)
    real(r8), pointer   :: hcond(:)
    real(r8), pointer   :: rofl(:), hrofl(:)
    real(r8), pointer   :: rofi(:), hrofi(:)
    real(r8), pointer   :: rofl_glc(:), hrofl_glc(:)
    real(r8), pointer   :: rofi_glc(:), hrofi_glc(:)
    real(r8), pointer   :: areas(:)
    real(r8), allocatable :: hcorr(:)
    type(med_fldlist_type), pointer :: fldList
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, maintask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    fldList => med_fldList_GetfldListTo(compocn)

    !---------------------------------------
    ! --- map atm to ocn, only if data stream is available
    !---------------------------------------
    if (is_local%wrap%med_coupling_active(compatm,compocn) .and. &
        is_local%wrap%med_data_active(compatm,compocn) .and. &
        is_local%wrap%med_data_force_first(compocn)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            FBDat=is_local%wrap%FBData(compocn), &
            use_data=is_local%wrap%med_data_force_first(compocn), &
            field_normOne=is_local%wrap%field_normOne(compatm,compocn,:), &
            packed_data=is_local%wrap%packed_data(compatm,compocn,:), &
            routehandles=is_local%wrap%RH(compatm,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ocn')

       ! Reset flag to use data
       is_local%wrap%med_data_force_first(compocn) = .false.
    end if

    !---------------------------------------
    !--- merge all fields to ocn
    !---------------------------------------
    call med_merge_auto(&
         is_local%wrap%med_coupling_active(:,compocn), &
         is_local%wrap%FBExp(compocn), &
         is_local%wrap%FBFrac(compocn), &
         is_local%wrap%FBImp(:,compocn), &
         fldList, &
         FBMed1=is_local%wrap%FBMed_aoflux_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------
    ! compute enthalpy associated with rain, snow, condensation and liquid river & glc runoff
    ! the sea-ice model already accounts for the enthalpy flux (as part of melth), so
    ! enthalpy from meltw **is not** included below
    if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_rain'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrain'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_snow'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hsnow'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hevap'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hcond'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc'  , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl_glc' , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc'  , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi_glc' , rc=rc)) then

       call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', tocn, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_rain' , rain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrain', hrain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap' , evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hevap', hevap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hcond', hcond, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_snow' , snow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', hsnow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', hrofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , rofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', hrofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc' , rofl_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofl_glc', hrofl_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc' , rofi_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofi_glc', hrofi_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do n = 1,size(tocn)
          ! Need max to ensure that will not have an enthalpy contribution if the water is below 0C
          hrain(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * rain(n)  * shr_const_cpsw
          hsnow(n)  = min((tocn(n) - shr_const_tkfrz), 0._r8) * snow(n)  * shr_const_cpsw
          hevap(n)  = (tocn(n) - shr_const_tkfrz) * min(evap(n), 0._r8)   * shr_const_cpsw
          hcond(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * max(evap(n), 0._r8)  * shr_const_cpsw
          hrofl(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * rofl(n)  * shr_const_cpsw
          hrofi(n)  = min((tocn(n) - shr_const_tkfrz), 0._r8) * rofi(n)  * shr_const_cpsw
          hrofl_glc(n) = max((tocn(n) - shr_const_tkfrz), 0._r8) * rofl_glc(n)  * shr_const_cpsw
          hrofi_glc(n) = min((tocn(n) - shr_const_tkfrz), 0._r8) * rofi_glc(n)  * shr_const_cpsw
       end do

       ! Determine enthalpy correction factor that will be added to the sensible heat flux sent to the atm
       ! Areas here in radians**2 - this is an instantaneous snapshot that will be sent to the atm - only
       ! need to calculate this if data is sent back to the atm

       if (FB_fldchk(is_local%wrap%FBExp(compatm), 'Faxx_sen', rc=rc)) then
          allocate(hcorr(size(tocn)))
          glob_area_inv = 1._r8 / (4._r8 * shr_const_pi)
          areas => is_local%wrap%mesh_info(compocn)%areas
          do n = 1,size(tocn)
             hcorr(n) = (hrain(n) + hsnow(n) + hcond(n) + hevap(n) + hrofl(n) + hrofi(n) + hrofl_glc(n) + hrofi_glc(n)) * &
                        areas(n) * glob_area_inv
          end do
          call med_phases_prep_atm_enthalpy_correction(gcomp, hcorr, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          deallocate(hcorr)
       end if

    end if

    ! custom merges to ocean
    call med_phases_prep_ocn_custom(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ocean accumulator
    call FB_accum(is_local%wrap%FBExpAccumOcn, is_local%wrap%FBExp(compocn), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%ExpAccumOcnCnt = is_local%wrap%ExpAccumOcnCnt + 1

    ! diagnose output
    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBExpAccumOcn, string=trim(subname)//' FBExpAccumOcn accumulation ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_avg(gcomp, rc)

    ! Prepare the OCN import Fields.

    use ESMF , only : ESMF_GridComp, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    integer                    :: ncnt
    logical, save              :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_prep_ocn_avg)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Count the number of fields outside of scalar data, if zero, then return
    call ESMF_FieldBundleGet(is_local%wrap%FBExpAccumOcn, fieldCount=ncnt, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       ! average ocn accumulator
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccumOcn, &
               string=trim(subname)//' FBExpAccumOcn before avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call FB_average(is_local%wrap%FBExpAccumOcn, is_local%wrap%ExpAccumOcnCnt, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccumOcn, &
               string=trim(subname)//' FBExpAccumOcn after avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! copy to FBExp(compocn)
       call FB_copy(is_local%wrap%FBExp(compocn), is_local%wrap%FBExpAccumOcn, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Check for nans in fields export to ocn
       if(.not. first_call) then
          call FB_check_for_nans(is_local%wrap%FBExp(compocn), maintask, logunit, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       ! zero accumulator
       is_local%wrap%ExpAccumOcnCnt = 0
       call FB_reset(is_local%wrap%FBExpAccumOcn, value=czero, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)
    first_call = .false.

  end subroutine med_phases_prep_ocn_avg

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_custom(gcomp, rc)

    !---------------------------------------
    ! custom calculations for cesm
    !---------------------------------------

    use ESMF , only : ESMF_GridComp, ESMF_StateGet, ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_VMBroadCast
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: lfield
    real(R8), pointer   :: ifrac(:)
    real(R8), pointer   :: ofrac(:)
    real(R8), pointer   :: ifracr(:)
    real(R8), pointer   :: ofracr(:)
    real(R8), pointer   :: avsdr(:)
    real(R8), pointer   :: avsdf(:)
    real(R8), pointer   :: anidr(:)
    real(R8), pointer   :: anidf(:)
    real(R8), pointer   :: Faxa_swvdf(:)
    real(R8), pointer   :: Faxa_swndf(:)
    real(R8), pointer   :: Faxa_swvdr(:)
    real(R8), pointer   :: Faxa_swndr(:)
    real(R8), pointer   :: Foxx_swnet(:)
    real(R8), pointer   :: Foxx_swnet_afracr(:)
    real(R8), pointer   :: Foxx_swnet_vdr(:)
    real(R8), pointer   :: Foxx_swnet_vdf(:)
    real(R8), pointer   :: Foxx_swnet_idr(:)
    real(R8), pointer   :: Foxx_swnet_idf(:)
    real(R8), pointer   :: Fioi_swpen_vdr(:)
    real(R8), pointer   :: Fioi_swpen_vdf(:)
    real(R8), pointer   :: Fioi_swpen_idr(:)
    real(R8), pointer   :: Fioi_swpen_idf(:)
    real(R8), pointer   :: Fioi_swpen(:)
    real(R8), pointer   :: dataptr(:)
    real(R8), pointer   :: dataptr_scalar_ocn(:,:)
    real(R8)            :: frac_sum
    real(R8)            :: ifrac_scaled, ofrac_scaled
    real(R8)            :: ifracr_scaled, ofracr_scaled
    logical             :: export_swnet_by_bands
    logical             :: import_swpen_by_bands
    logical             :: export_swnet_afracr
    real(R8)            :: precip_fact(1)
    character(CS)       :: cvalue
    real(R8)            :: fswabsv, fswabsi
    integer             :: scalar_id
    integer             :: n
    integer             :: lsize
    real(R8)            :: c1,c2,c3,c4
    character(len=64), allocatable :: fldnames(:)
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_custom)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 5, maintask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check that the necessary export field is present
    if ( .not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc) .and. &
         .not. (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', rc=rc))) then
       return
    end if

    call t_startf('MED:'//subname)

    !---------------------------------------
    ! Compute netsw for ocean
    !---------------------------------------
    ! netsw_for_ocn = downsw_from_atm * (1-ocn_albedo) * (1-ice_fraction) + pensw_from_ice * (ice_fraction)

    ! Input from atm
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', Faxa_swvdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', Faxa_swndr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', Faxa_swvdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', Faxa_swndf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    lsize = size(Faxa_swvdr)

    ! Input from mediator, ocean albedos
    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Output to ocean swnet total
    if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  Foxx_swnet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       lsize = size(Faxa_swvdr)
       allocate(Foxx_swnet(lsize))
    end if

    ! Output to ocean swnet by radiation bands
    if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc)) then
       export_swnet_by_bands = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', Foxx_swnet_vdr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', Foxx_swnet_vdf, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', Foxx_swnet_idr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', Foxx_swnet_idf, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       export_swnet_by_bands = .false.
    end if

    ! -----------------------
    ! If cice IS NOT PRESENT
    ! -----------------------
    if (.not. is_local%wrap%comp_present(compice)) then
       ! Compute total swnet to ocean independent of swpen from sea-ice
       do n = 1,lsize
          fswabsv  = Faxa_swvdr(n) * (1.0_R8 - avsdr(n)) + Faxa_swvdf(n) * (1.0_R8 - avsdf(n))
          fswabsi  = Faxa_swndr(n) * (1.0_R8 - anidr(n)) + Faxa_swndf(n) * (1.0_R8 - anidf(n))
          Foxx_swnet(n) = fswabsv + fswabsi
       end do
       ! Compute sw export to ocean bands if required
       if (export_swnet_by_bands) then
          if (trim(coupling_mode) == 'cesm') then
             c1 = 0.285; c2 = 0.285; c3 = 0.215; c4 = 0.215
             Foxx_swnet_vdr(:) = c1 * Foxx_swnet(:)
             Foxx_swnet_vdf(:) = c2 * Foxx_swnet(:)
             Foxx_swnet_idr(:) = c3 * Foxx_swnet(:)
             Foxx_swnet_idf(:) = c4 * Foxx_swnet(:)
          else
             Foxx_swnet_vdr(:) = Faxa_swvdr(:) * (1.0_R8 - avsdr(:))
             Foxx_swnet_vdf(:) = Faxa_swvdf(:) * (1.0_R8 - avsdf(:))
             Foxx_swnet_idr(:) = Faxa_swndr(:) * (1.0_R8 - anidr(:))
             Foxx_swnet_idf(:) = Faxa_swndf(:) * (1.0_R8 - anidf(:))
          end if
       end if
    end if

    ! -----------------------
    ! If cice IS PRESENT
    ! -----------------------
    if (is_local%wrap%comp_present(compice)) then

       ! Input from mediator, ice-covered ocean and open ocean fractions
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', Fioi_swpen, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (FB_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc)) then
          import_swpen_by_bands = .true.
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          import_swpen_by_bands = .false.
       end if

       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr',rc=rc)) then
          ! Swnet without swpen from sea-ice
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr', Foxx_swnet_afracr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          export_swnet_afracr = .true.
       else
          export_swnet_afracr = .false.
       end if

       do n = 1,lsize
          ! Compute total swnet to ocean independent of swpen from sea-ice
          fswabsv  = Faxa_swvdr(n) * (1.0_R8 - avsdr(n)) + Faxa_swvdf(n) * (1.0_R8 - avsdf(n))
          fswabsi  = Faxa_swndr(n) * (1.0_R8 - anidr(n)) + Faxa_swndf(n) * (1.0_R8 - anidf(n))
          Foxx_swnet(n) = fswabsv + fswabsi

          ! Add swpen from sea ice
          ifrac_scaled = ifrac(n)
          ofrac_scaled = ofrac(n)
          frac_sum = ifrac(n) + ofrac(n)
          if (frac_sum /= 0._R8) then
             ifrac_scaled = ifrac(n) / (frac_sum)
             ofrac_scaled = ofrac(n) / (frac_sum)
          endif
          ifracr_scaled = ifracr(n)
          ofracr_scaled = ofracr(n)
          frac_sum = ifracr(n) + ofracr(n)
          if (frac_sum /= 0._R8) then
             ifracr_scaled = ifracr(n) / (frac_sum)
             ofracr_scaled = ofracr(n) / (frac_sum)
          endif
          Foxx_swnet(n) = ofracr_scaled*(fswabsv + fswabsi) + ifrac_scaled*Fioi_swpen(n)

          if (export_swnet_afracr) then
             Foxx_swnet_afracr(n) = ofracr_scaled*(fswabsv + fswabsi)
          end if

          ! Compute sw export to ocean bands if required
          if (export_swnet_by_bands) then
             if (import_swpen_by_bands) then
                ! use each individual band for swpen coming from the sea-ice
                Foxx_swnet_vdr(n) = Faxa_swvdr(n)*(1.0_R8-avsdr(n))*ofracr_scaled + Fioi_swpen_vdr(n)*ifrac_scaled
                Foxx_swnet_vdf(n) = Faxa_swvdf(n)*(1.0_R8-avsdf(n))*ofracr_scaled + Fioi_swpen_vdf(n)*ifrac_scaled
                Foxx_swnet_idr(n) = Faxa_swndr(n)*(1.0_R8-anidr(n))*ofracr_scaled + Fioi_swpen_idr(n)*ifrac_scaled
                Foxx_swnet_idf(n) = Faxa_swndf(n)*(1.0_R8-anidf(n))*ofracr_scaled + Fioi_swpen_idf(n)*ifrac_scaled
             else
                ! scale total Foxx_swnet to get contributions from each band
                c1 = 0.285; c2 = 0.285; c3 = 0.215; c4 = 0.215
                Foxx_swnet_vdr(n) = c1 * Foxx_swnet(n)
                Foxx_swnet_vdf(n) = c2 * Foxx_swnet(n)
                Foxx_swnet_idr(n) = c3 * Foxx_swnet(n)
                Foxx_swnet_idf(n) = c4 * Foxx_swnet(n)
             end if
          end if
       end do

       ! Output to ocean per ice thickness fraction and sw penetrating into ocean
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afrac', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afrac', fldptr1=dataptr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr(:) = ofrac(:)
       end if
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afracr', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afracr', fldptr1=dataptr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr(:) = ofracr(:)
       end if

    end if  ! if sea-ice is present

    ! Deallocate Foxx_swnet if it was allocated in this subroutine
    if (.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
       deallocate(Foxx_swnet)
    end if

    ! Apply precipitation factor from ocean (that scales atm rain and snow back to ocn ) if appropriate
    if (trim(coupling_mode) == 'cesm' .and. is_local%wrap%flds_scalar_index_precip_factor /= 0) then

       ! Note that in med_internal_mod.F90 all is_local%wrap%flds_scalar_index_precip_factor
       ! is initialized to 0.
       ! In addition, in med.F90, if this attribute is not present as a mediator component attribute,
       ! it is set to 0.
       if (maintask) then
          call ESMF_StateGet(is_local%wrap%NstateImp(compocn), &
               itemName=trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_ocn, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          scalar_id=is_local%wrap%flds_scalar_index_precip_factor
          precip_fact(1) = dataptr_scalar_ocn(scalar_id,1)
          if (precip_fact(1) /= 1._r8) then
             write(logunit,'(a,f21.13)')&
                  '(merge_to_ocn): Scaling rain, snow, liquid and ice runoff by non-unity precip_fact ',&
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

       ! Scale rain and snow to ocn from atm by the precipitation factor received from the ocean
       allocate(fldnames(4))
       fldnames = (/'Faxa_rain', 'Faxa_snow', 'Foxx_rofl', 'Foxx_rofi'/)
       do n = 1,size(fldnames)
          if (FB_fldchk(is_local%wrap%FBExp(compocn), trim(fldnames(n)), rc=rc)) then
             call FB_GetFldPtr(is_local%wrap%FBExp(compocn), trim(fldnames(n)) , dataptr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             dataptr(:) = dataptr(:) * is_local%wrap%flds_scalar_precip_factor
          end if
       end do
       deallocate(fldnames)
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_custom

end module med_phases_prep_ocn_mod
