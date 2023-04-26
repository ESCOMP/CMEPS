module med_phases_prep_wav_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing wav export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod     , only : czero     =>med_constants_czero
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_merge_mod         , only : med_merge_auto, med_merge_field
  use med_map_mod           , only : med_map_field_packed
  use med_utils_mod         , only : memcheck      => med_memcheck
  use med_utils_mod         , only : chkerr        => med_utils_ChkErr
  use med_methods_mod       , only : FB_diagnose   => med_methods_FB_diagnose
!PSH begin
  use med_methods_mod       , only : FB_fldchk     => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr
!PSH end
  use med_methods_mod       , only : FB_accum      => med_methods_FB_accum
  use med_methods_mod       , only : FB_average    => med_methods_FB_average
  use med_methods_mod       , only : FB_copy       => med_methods_FB_copy
  use med_methods_mod       , only : FB_reset      => med_methods_FB_reset
!PSH begin
  use esmFlds               , only : med_fldList_GetfldListTo
  use med_internalstate_mod , only : compwav
!  use esmFlds               , only : med_fldList_GetfldListTo, med_fldlist_type
!  use med_internalstate_mod , only : compwav, compocn, compatm, compice, coupling_mode
!PSH end
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_prep_wav_init   ! called from med.F90
  public :: med_phases_prep_wav_accum  ! called from run sequence
  public :: med_phases_prep_wav_avg    ! called from run sequence

!PSH begin
!  private :: med_phases_prep_wav_custom_cesm
!PSH end

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_wav_init(gcomp, rc)

    use ESMF            , only : ESMF_GridComp, ESMF_SUCCESS
    use med_methods_mod , only : FB_Init  => med_methods_FB_init

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    character(len=*),parameter  :: subname=' (med_phases_prep_wav_init) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       write(logunit,'(a)') trim(subname)//' initializing wave export accumulation FB for '
    end if
    call FB_Init(is_local%wrap%FBExpAccumWav, is_local%wrap%flds_scalar_name, &
         STgeom=is_local%wrap%NStateExp(compwav), STflds=is_local%wrap%NStateExp(compwav), &
         name='FBExpAccumWav', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_reset(is_local%wrap%FBExpAccumWav, value=czero, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_prep_wav_init

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_wav_accum(gcomp, rc)

    use ESMF , only : ESMF_GridComp, ESMF_FieldBundleGet
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n, ncnt
!PSH begin
!    type(med_fldlist_type), pointer :: fldList
!PSH end
    character(len=*), parameter    :: subname='(med_phases_prep_wav_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, mastertask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!PSH begin
!    fldList => med_fldList_GetfldListTo(compwav)
!PSH end
    ! auto merges to wav
!PSH begin
    call med_merge_auto(&
         is_local%wrap%med_coupling_active(:,compwav), &
         is_local%wrap%FBExp(compwav), &
         is_local%wrap%FBFrac(compwav), &
         is_local%wrap%FBImp(:,compwav), &
         med_fldList_GetfldListTo(compwav), rc=rc)
!       call med_merge_auto(&
!            is_local%wrap%med_coupling_active(:,compwav), &
!            is_local%wrap%FBExp(compwav), &
!            is_local%wrap%FBFrac(compwav), &
!            is_local%wrap%FBImp(:,compwav), &
!            fldList, &
!            FBMed1=is_local%wrap%FBMed_aoflux_o, rc=rc)
!PSH end
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!PSH begin
!    ! custom merges to ocean
!    if (trim(coupling_mode) == 'cesm') then
!       call med_phases_prep_wav_custom_cesm(gcomp, rc) 
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    end if
!PSH end

    ! wave accumulator
    call FB_accum(is_local%wrap%FBExpAccumWav, is_local%wrap%FBExp(compwav), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%ExpAccumWavCnt = is_local%wrap%ExpAccumWavCnt + 1

    ! diagnose output
    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBExpAccumWav, string=trim(subname)//' FBExpAccumWav accumulation ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_wav_accum

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_wav_avg(gcomp, rc)

    ! Prepare the wav import Fields.

    use ESMF , only : ESMF_GridComp, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    integer                    :: ncnt
    character(len=*),parameter :: subname='(med_phases_prep_wav_avg)'
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
    call ESMF_FieldBundleGet(is_local%wrap%FBExpAccumWav, fieldCount=ncnt, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       ! average wav accumulator
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccumWav, &
               string=trim(subname)//' FBExpAccumWav before avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call FB_average(is_local%wrap%FBExpAccumWav, is_local%wrap%ExpAccumWavCnt, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccumWav, &
               string=trim(subname)//' FBExpAccumWav after avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! copy to FBExp(compwav)
       call FB_copy(is_local%wrap%FBExp(compwav), is_local%wrap%FBExpAccumWav, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! zero accumulator
       is_local%wrap%ExpAccumWavCnt = 0
       call FB_reset(is_local%wrap%FBExpAccumWav, value=czero, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_wav_avg
  !-----------------------------------------------------------------------------
!  subroutine med_phases_prep_wav_custom_cesm(gcomp, rc)
!
!    !---------------------------------------
!    ! custom calculations for cesm
!    !---------------------------------------
!
!    use ESMF , only : ESMF_GridComp, ESMF_StateGet, ESMF_Field, ESMF_FieldGet
!    use ESMF , only : ESMF_VMBroadCast
!    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
!    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR
!
!    ! input/output variables
!    type(ESMF_GridComp)  :: gcomp
!    integer, intent(out) :: rc
!
!    ! local variables
!    type(InternalState) :: is_local
!    type(ESMF_Field)    :: lfield
!    real(R8), pointer   :: ifrac(:)
!    real(R8), pointer   :: ofrac(:)
!    real(R8), pointer   :: ifracr(:)
!    real(R8), pointer   :: ofracr(:)
!    real(R8), pointer   :: avsdr(:)
!    real(R8), pointer   :: avsdf(:)
!    real(R8), pointer   :: anidr(:)
!    real(R8), pointer   :: anidf(:)
!    real(R8), pointer   :: Faxa_swvdf(:)
!    real(R8), pointer   :: Faxa_swndf(:)
!    real(R8), pointer   :: Faxa_swvdr(:)
!    real(R8), pointer   :: Faxa_swndr(:)
!    real(R8), pointer   :: Foxx_swnet(:)
!    real(R8), pointer   :: Foxx_swnet_afracr(:)
!    real(R8), pointer   :: Foxx_swnet_vdr(:)
!    real(R8), pointer   :: Foxx_swnet_vdf(:)
!    real(R8), pointer   :: Foxx_swnet_idr(:)
!    real(R8), pointer   :: Foxx_swnet_idf(:)
!    real(R8), pointer   :: Fioi_swpen_vdr(:)
!    real(R8), pointer   :: Fioi_swpen_vdf(:)
!    real(R8), pointer   :: Fioi_swpen_idr(:)
!    real(R8), pointer   :: Fioi_swpen_idf(:)
!    real(R8), pointer   :: Fioi_swpen(:)
!    real(R8), pointer   :: dataptr(:)
!    real(R8), pointer   :: dataptr_scalar_ocn(:,:)
!    real(R8)            :: frac_sum
!    real(R8)            :: ifrac_scaled, ofrac_scaled
!    real(R8)            :: ifracr_scaled, ofracr_scaled
!    logical             :: export_swnet_by_bands
!    logical             :: import_swpen_by_bands
!    logical             :: export_swnet_afracr
!    real(R8)            :: precip_fact(1)
!    character(CS)       :: cvalue
!    real(R8)            :: fswabsv, fswabsi
!    integer             :: scalar_id
!    integer             :: n
!    integer             :: lsize
!    real(R8)            :: c1,c2,c3,c4
!    character(len=64), allocatable :: fldnames(:)
!    character(len=*), parameter    :: subname='(med_phases_prep_wav_custom_cesm)'
!    !---------------------------------------
!
!    rc = ESMF_SUCCESS
!
!    call t_startf('MED:'//subname)
!    if (dbug_flag > 20) then
!       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
!    end if
!    call memcheck(subname, 5, mastertask)
!
!    ! Get the internal state
!    nullify(is_local%wrap)
!    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!
!    !---------------------------------------
!    ! Compute netsw for ocean
!    !---------------------------------------
!    ! netsw_for_ocn = downsw_from_atm * (1-ocn_albedo) * (1-ice_fraction) + pensw_from_ice * (ice_fraction)
!
!    ! Input from atm
!    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', Faxa_swvdr, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', Faxa_swndr, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', Faxa_swvdf, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', Faxa_swndf, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    lsize = size(Faxa_swvdr)
!
!    ! Input from mediator, ocean albedos
!    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!
!    ! Output to ocean swnet total
!    if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
!       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  Foxx_swnet, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    else
!       lsize = size(Faxa_swvdr)
!       allocate(Foxx_swnet(lsize))
!    end if
!
!    ! Output to ocean swnet by radiation bands
!    if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc)) then
!       export_swnet_by_bands = .true.
!       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', Foxx_swnet_vdr, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', Foxx_swnet_vdf, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', Foxx_swnet_idr, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', Foxx_swnet_idf, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    else
!       export_swnet_by_bands = .false.
!    end if
!
!    ! -----------------------
!    ! If cice IS NOT PRESENT
!    ! -----------------------
!    if (.not. is_local%wrap%comp_present(compice)) then
!       ! Compute total swnet to ocean independent of swpen from sea-ice
!       do n = 1,lsize
!          fswabsv  = Faxa_swvdr(n) * (1.0_R8 - avsdr(n)) + Faxa_swvdf(n) * (1.0_R8 - avsdf(n))
!          fswabsi  = Faxa_swndr(n) * (1.0_R8 - anidr(n)) + Faxa_swndf(n) * (1.0_R8 - anidf(n))
!          Foxx_swnet(n) = fswabsv + fswabsi
!       end do
!       ! Compute sw export to ocean bands if required
!       if (export_swnet_by_bands) then
!          c1 = 0.285; c2 = 0.285; c3 = 0.215; c4 = 0.215
!          Foxx_swnet_vdr(:) = c1 * Foxx_swnet(:)
!          Foxx_swnet_vdf(:) = c2 * Foxx_swnet(:)
!          Foxx_swnet_idr(:) = c3 * Foxx_swnet(:)
!          Foxx_swnet_idf(:) = c4 * Foxx_swnet(:)
!       end if
!    end if
!
!    ! -----------------------
!    ! If cice IS PRESENT
!    ! -----------------------
!    if (is_local%wrap%comp_present(compice)) then
!
!       ! Input from mediator, ice-covered ocean and open ocean fractions
!       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!
!       call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', Fioi_swpen, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       if (FB_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc)) then
!          import_swpen_by_bands = .true.
!          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       else
!          import_swpen_by_bands = .false.
!       end if
!
!       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr',rc=rc)) then
!          ! Swnet without swpen from sea-ice
!          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr', Foxx_swnet_afracr, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          export_swnet_afracr = .true.
!       else
!          export_swnet_afracr = .false.
!       end if
!
!       do n = 1,lsize
!          ! Compute total swnet to ocean independent of swpen from sea-ice
!          fswabsv  = Faxa_swvdr(n) * (1.0_R8 - avsdr(n)) + Faxa_swvdf(n) * (1.0_R8 - avsdf(n))
!          fswabsi  = Faxa_swndr(n) * (1.0_R8 - anidr(n)) + Faxa_swndf(n) * (1.0_R8 - anidf(n))
!          Foxx_swnet(n) = fswabsv + fswabsi
!
!          ! Add swpen from sea ice
!          ifrac_scaled = ifrac(n)
!          ofrac_scaled = ofrac(n)
!          frac_sum = ifrac(n) + ofrac(n)
!          if (frac_sum /= 0._R8) then
!             ifrac_scaled = ifrac(n) / (frac_sum)
!             ofrac_scaled = ofrac(n) / (frac_sum)
!          endif
!          ifracr_scaled = ifracr(n)
!          ofracr_scaled = ofracr(n)
!          frac_sum = ifracr(n) + ofracr(n)
!          if (frac_sum /= 0._R8) then
!             ifracr_scaled = ifracr(n) / (frac_sum)
!             ofracr_scaled = ofracr(n) / (frac_sum)
!          endif
!          Foxx_swnet(n) = ofracr_scaled*(fswabsv + fswabsi) + ifrac_scaled*Fioi_swpen(n)
!
!          if (export_swnet_afracr) then
!             Foxx_swnet_afracr(n) = ofracr_scaled*(fswabsv + fswabsi)
!          end if
!
!          ! Compute sw export to ocean bands if required
!          if (export_swnet_by_bands) then
!             if (import_swpen_by_bands) then
!                ! use each individual band for swpen coming from the sea-ice
!                Foxx_swnet_vdr(n) = Faxa_swvdr(n)*(1.0_R8-avsdr(n))*ofracr_scaled + Fioi_swpen_vdr(n)*ifrac_scaled
!                Foxx_swnet_vdf(n) = Faxa_swvdf(n)*(1.0_R8-avsdf(n))*ofracr_scaled + Fioi_swpen_vdf(n)*ifrac_scaled
!                Foxx_swnet_idr(n) = Faxa_swndr(n)*(1.0_R8-anidr(n))*ofracr_scaled + Fioi_swpen_idr(n)*ifrac_scaled
!                Foxx_swnet_idf(n) = Faxa_swndf(n)*(1.0_R8-anidf(n))*ofracr_scaled + Fioi_swpen_idf(n)*ifrac_scaled
!             else
!                ! scale total Foxx_swnet to get contributions from each band
!                c1 = 0.285; c2 = 0.285; c3 = 0.215; c4 = 0.215
!                Foxx_swnet_vdr(n) = c1 * Foxx_swnet(n)
!                Foxx_swnet_vdf(n) = c2 * Foxx_swnet(n)
!                Foxx_swnet_idr(n) = c3 * Foxx_swnet(n)
!                Foxx_swnet_idf(n) = c4 * Foxx_swnet(n)
!             end if
!          end if
!       end do
!
!       ! Output to ocean per ice thickness fraction and sw penetrating into ocean
!       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afrac', rc=rc)) then
!          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afrac', fldptr1=dataptr, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          dataptr(:) = ofrac(:)
!       end if
!       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afracr', rc=rc)) then
!          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afracr', fldptr1=dataptr, rc=rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!          dataptr(:) = ofracr(:)
!       end if
!
!    end if  ! if sea-ice is present
!
!    ! Deallocate Foxx_swnet if it was allocated in this subroutine
!    if (.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
!       deallocate(Foxx_swnet)
!    end if
!
!    ! Apply precipitation factor from ocean (that scales atm rain and snow back to ocn ) if appropriate
!    if (trim(coupling_mode) == 'cesm' .and. is_local%wrap%flds_scalar_index_precip_factor /= 0) then
!
!       ! Note that in med_internal_mod.F90 all is_local%wrap%flds_scalar_index_precip_factor
!       ! is initialized to 0.
!       ! In addition, in med.F90, if this attribute is not present as a mediator component attribute,
!       ! it is set to 0.
!       if (mastertask) then
!          call ESMF_StateGet(is_local%wrap%NstateImp(compocn), &
!               itemName=trim(is_local%wrap%flds_scalar_name), field=lfield, rc=rc)
!          if (chkerr(rc,__LINE__,u_FILE_u)) return
!          call ESMF_FieldGet(lfield, farrayPtr=dataptr_scalar_ocn, rc=rc)
!          if (chkerr(rc,__LINE__,u_FILE_u)) return
!          scalar_id=is_local%wrap%flds_scalar_index_precip_factor
!          precip_fact(1) = dataptr_scalar_ocn(scalar_id,1)
!          if (precip_fact(1) /= 1._r8) then
!             write(logunit,'(a,f21.13)')&
!                  '(merge_to_ocn): Scaling rain, snow, liquid and ice runoff by non-unity precip_fact ',&
!                  precip_fact(1)
!          end if
!       end if
!       call ESMF_VMBroadCast(is_local%wrap%vm, precip_fact, 1, 0, rc=rc)
!       if (chkerr(rc,__LINE__,u_FILE_u)) return
!       is_local%wrap%flds_scalar_precip_factor = precip_fact(1)
!       if (dbug_flag > 5) then
!          write(cvalue,*) precip_fact(1)
!          call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO)
!       end if
!
!       ! Scale rain and snow to ocn from atm by the precipitation factor received from the ocean
!       allocate(fldnames(4))
!       fldnames = (/'Faxa_rain', 'Faxa_snow', 'Foxx_rofl', 'Foxx_rofi'/)
!       do n = 1,size(fldnames)
!          if (FB_fldchk(is_local%wrap%FBExp(compocn), trim(fldnames(n)), rc=rc)) then
!             call FB_GetFldPtr(is_local%wrap%FBExp(compocn), trim(fldnames(n)) , dataptr, rc=rc)
!             if (ChkErr(rc,__LINE__,u_FILE_u)) return
!             dataptr(:) = dataptr(:) * is_local%wrap%flds_scalar_precip_factor
!          end if
!       end do
!       deallocate(fldnames)
!    end if
!
!    if (dbug_flag > 20) then
!       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
!    end if
!    call t_stopf('MED:'//subname)
!
!  end subroutine med_phases_prep_wav_custom_cesm

end module med_phases_prep_wav_mod
