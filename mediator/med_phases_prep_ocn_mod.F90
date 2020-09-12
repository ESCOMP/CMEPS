module med_phases_prep_ocn_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ocn export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod     , only : czero=>med_constants_czero
  use med_constants_mod     , only : dbug_flag     => med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_merge_mod         , only : med_merge_auto, med_merge_field
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_utils_mod         , only : memcheck      => med_memcheck
  use med_utils_mod         , only : chkerr        => med_utils_ChkErr
  use med_methods_mod       , only : FB_diagnose   => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_getNumFlds => med_methods_FB_getNumFlds
  use med_methods_mod       , only : FB_fldchk     => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_accum      => med_methods_FB_accum
  use med_methods_mod       , only : FB_average    => med_methods_FB_average
  use med_methods_mod       , only : FB_copy       => med_methods_FB_copy
  use med_methods_mod       , only : FB_reset      => med_methods_FB_reset
  use esmFlds               , only : fldListFr, fldListTo
  use esmFlds               , only : compocn, compatm, compice, ncomps, compname
  use esmFlds               , only : coupling_mode
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_prep_ocn_map
  public :: med_phases_prep_ocn_merge
  public :: med_phases_prep_ocn_accum_fast
  public :: med_phases_prep_ocn_accum_avg

  private :: med_phases_prep_ocn_custom_cesm
  private :: med_phases_prep_ocn_custom_nems
  private :: med_phases_prep_ocn_custom_nemsdata

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_map(gcomp, rc)

    !---------------------------------------
    ! Map all fields in from relevant source components to the ocean grid
    !---------------------------------------

    use ESMF , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    integer                     :: n1, ncnt
    character(len=*), parameter :: subname='(med_phases_prep_ocn_map)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 5, mastertask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Count the number of fields outside of scalar data, if zero, then return
    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       ! map all fields in FBImp that have active ocean coupling to the ocean grid
       do n1 = 1,ncomps
          if (is_local%wrap%med_coupling_active(n1,compocn)) then
             call med_map_FB_Regrid_Norm( &
                  fldsSrc=fldListFr(n1)%flds,&
                  srccomp=n1, destcomp=compocn, &
                  FBSrc=is_local%wrap%FBImp(n1,n1), &
                  FBDst=is_local%wrap%FBImp(n1,compocn), &
                  FBFracSrc=is_local%wrap%FBFrac(n1), &
                  FBNormOne=is_local%wrap%FBNormOne(n1,compocn,:), &
                  RouteHandles=is_local%wrap%RH(n1,compocn,:), &
                  string=trim(compname(n1))//'2'//trim(compname(compocn)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo
    endif

    call t_stopf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if

  end subroutine med_phases_prep_ocn_map

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_merge(gcomp, rc)

    use ESMF , only : ESMF_GridComp, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n, ncnt
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_merge)'
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

    ! Count the number of fields outside of scalar data, if zero, then return
    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       ! merges to ocean
       !---------------------------------------

       ! auto merges to ocn
       if (trim(coupling_mode) == 'cesm' .or. &
           trim(coupling_mode) == 'nems_orig_data' .or. &
           trim(coupling_mode) == 'hafs') then
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), &
               FBMed1=is_local%wrap%FBMed_aoflux_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(coupling_mode) == 'nems_frac' .or. trim(coupling_mode) == 'nems_orig') then
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! custom merges to ocean
       if (trim(coupling_mode) == 'cesm') then
          call med_phases_prep_ocn_custom_cesm(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(coupling_mode) == 'nems_orig' .or. trim(coupling_mode) == 'nems_frac') then
          call med_phases_prep_ocn_custom_nems(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(coupling_mode) == 'nems_orig_data') then
          call med_phases_prep_ocn_custom_nemsdata(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if  ! end of nems_orig_data custom

       ! diagnose output
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compocn), string=trim(subname)//' FBexp(compocn) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    endif

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_merge

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_accum_fast(gcomp, rc)

    ! Carry out fast accumulation for the ocean

    use ESMF , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_Clock, ESMF_Time
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    type(InternalState)         :: is_local
    integer                     :: i,j,n,ncnt
    character(len=*), parameter :: subname='(med_phases_accum_fast)'
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
    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then
       ! ocean accumulator
       call FB_accum(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExp(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       is_local%wrap%FBExpAccumCnt(compocn) = is_local%wrap%FBExpAccumCnt(compocn) + 1

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum accumulation ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_fast

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_accum_avg(gcomp, rc)

    ! Prepare the OCN import Fields.

    use ESMF , only : ESMF_GridComp, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    integer                    :: ncnt
    character(len=*),parameter :: subname='(med_phases_prep_ocn_accum_avg)'
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
    call FB_getNumFlds(is_local%wrap%FBExpAccum(compocn), trim(subname)//"FBExpAccum(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       ! average ocn accumulator
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) before avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call FB_average(is_local%wrap%FBExpAccum(compocn), &
            is_local%wrap%FBExpAccumCnt(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) after avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! copy to FBExp(compocn)
       call FB_copy(is_local%wrap%FBExp(compocn), is_local%wrap%FBExpAccum(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! zero accumulator
       is_local%wrap%FBExpAccumFlag(compocn) = .true.
       is_local%wrap%FBExpAccumCnt(compocn) = 0
       call FB_reset(is_local%wrap%FBExpAccum(compocn), value=czero, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_avg

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_custom_cesm(gcomp, rc)

    !---------------------------------------
    ! custom calculations for cesm
    !---------------------------------------

    use ESMF , only : ESMF_GridComp
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    real(R8), pointer   :: ifrac(:), ofrac(:)
    real(R8), pointer   :: ifracr(:), ofracr(:)
    real(R8), pointer   :: avsdr(:), avsdf(:)
    real(R8), pointer   :: anidr(:), anidf(:)
    real(R8), pointer   :: Faxa_swvdf(:), Faxa_swndf(:)
    real(R8), pointer   :: Faxa_swvdr(:), Faxa_swndr(:)
    real(R8), pointer   :: Foxx_swnet(:)
    real(R8), pointer   :: Foxx_swnet_afracr(:)
    real(R8), pointer   :: Foxx_swnet_vdr(:), Foxx_swnet_vdf(:)
    real(R8), pointer   :: Foxx_swnet_idr(:), Foxx_swnet_idf(:)
    real(R8), pointer   :: Fioi_swpen_vdr(:), Fioi_swpen_vdf(:)
    real(R8), pointer   :: Fioi_swpen_idr(:), Fioi_swpen_idf(:)
    real(R8), pointer   :: Fioi_swpen(:)
    real(R8), pointer   :: dataptr(:)
    real(R8), pointer   :: dataptr_o(:)
    real(R8)            :: frac_sum
    real(R8)            :: ifrac_scaled, ofrac_scaled
    real(R8)            :: ifracr_scaled, ofracr_scaled
    logical             :: export_swnet_by_bands
    logical             :: import_swpen_by_bands
    logical             :: export_swnet_afracr
    logical             :: first_precip_fact_call = .true.
    real(R8)            :: precip_fact
    character(CS)       :: cvalue
    real(R8)            :: fswabsv, fswabsi
    integer             :: n
    integer             :: lsize
    real(R8)            :: c1,c2,c3,c4
    character(len=64), allocatable :: fldnames(:)
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_custom_cesm)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 5, mastertask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
          c1 = 0.285; c2 = 0.285; c3 = 0.215; c4 = 0.215
          Foxx_swnet_vdr(:) = c1 * Foxx_swnet(:)
          Foxx_swnet_vdf(:) = c2 * Foxx_swnet(:)
          Foxx_swnet_idr(:) = c3 * Foxx_swnet(:)
          Foxx_swnet_idf(:) = c4 * Foxx_swnet(:)
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

       ! Swnet without swpen from sea-ice
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr',rc=rc)) then
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
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afrac', fldptr1=dataptr_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr_o(:) = ofrac(:)
       end if
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afracr', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afracr', fldptr1=dataptr_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr_o(:) = ofracr(:)
       end if

    end if  ! if sea-ice is present

    ! Deallocate Foxx_swnet if it was allocated in this subroutine
    if (.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
       deallocate(Foxx_swnet)
    end if

    !---------------------------------------
    ! application of precipitation factor from ocean
    !---------------------------------------
    precip_fact = 1.0_R8
    if (precip_fact /= 1.0_R8) then
       if (first_precip_fact_call .and. mastertask) then
          write(logunit,'(a)')'(merge_to_ocn): Scaling rain, snow, liquid and ice runoff by precip_fact '
          first_precip_fact_call = .false.
       end if
       write(cvalue,*) precip_fact
       call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO)

       allocate(fldnames(4))
       fldnames = (/'Faxa_rain','Faxa_snow', 'Foxx_rofl', 'Foxx_rofi'/)
       do n = 1,size(fldnames)
          if (FB_fldchk(is_local%wrap%FBExp(compocn), trim(fldnames(n)), rc=rc)) then
             call FB_GetFldPtr(is_local%wrap%FBExp(compocn), trim(fldnames(n)) , dataptr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             dataptr(:) = dataptr(:) * precip_fact
          end if
       end do
       deallocate(fldnames)
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_custom_cesm

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_custom_nems(gcomp, rc)

    ! ----------------------------------------------
    ! Custom calculation for nems_orig or nems_frac
    ! ----------------------------------------------

    use ESMF , only : ESMF_GridComp
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    real(R8), pointer   :: ocnwgt1(:)
    real(R8), pointer   :: icewgt1(:)
    real(R8), pointer   :: wgtp01(:)
    real(R8), pointer   :: wgtm01(:)
    real(R8), pointer   :: customwgt(:)
    real(R8), pointer   :: ifrac(:)
    real(R8), pointer   :: ofrac(:)
    integer             :: lsize
    real(R8)        , parameter    :: const_lhvap = 2.501e6_R8  ! latent heat of evaporation ~ J/kg
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_custom_nems)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 5, mastertask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get ice and open ocean fractions on the ocn mesh
    call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(ofrac)
    allocate(customwgt(lsize))

    !call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_rain',  &
    !     FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_rain' , wgtA=ofrac, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_snow',  &
    !     FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_snow' , wgtA=ofrac, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_lwnet',  &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_lwnet', wgtA=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    customwgt(:) = -ofrac(:) / const_lhvap
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_evap', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_lat' , wgtA=customwgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    customwgt(:) = -ofrac(:)
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_sen',  &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_sen', wgtA=customwgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_taux',  &
         FBinA=is_local%wrap%FBImp(compice,compocn), fnameA='Fioi_taux' , wgtA=ifrac, &
         FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_taux' , wgtB=customwgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_tauy',  &
         FBinA=is_local%wrap%FBImp(compice,compocn), fnameA='Fioi_tauy' , wgtA=ifrac, &
         FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_tauy' , wgtB=customwgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! netsw_for_ocn = [downsw_from_atm*(1-ice_fraction)*(1-ocn_albedo)] + [pensw_from_ice*(ice_fraction)]
    customwgt(:) = ofrac(:) * (1.0 - 0.06)
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_vdr', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swvdr'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_vdr', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_vdf', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swvdf'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_vdf', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_idr', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swndr'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_idr', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_idf', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swndf'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_idf', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(customwgt)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_custom_nems

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ocn_custom_nemsdata(gcomp, rc)

    ! ----------------------------------------------
    ! Custom calculation for nems_orig_data
    ! ----------------------------------------------

    use ESMF , only : ESMF_GridComp
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    real(R8), pointer   :: ocnwgt1(:)   ! NEMS_orig_data
    real(R8), pointer   :: icewgt1(:)   ! NEMS_orig_data
    real(R8), pointer   :: wgtp01(:)    ! NEMS_orig_data
    real(R8), pointer   :: wgtm01(:)    ! NEMS_orig_data
    real(R8), pointer   :: customwgt(:) ! NEMS_orig_data
    real(R8), pointer   :: ifrac(:)
    real(R8), pointer   :: ofrac(:)
    integer             :: lsize
    integer             :: n
    real(R8)        , parameter    :: const_lhvap = 2.501e6_R8  ! latent heat of evaporation ~ J/kg
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_custom_nemsdata)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if
    call memcheck(subname, 5, mastertask)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get ice and open ocean fractions on the ocn mesh
    call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(ofrac)
    allocate(customwgt(lsize))

    ! open ocean (i.e. atm)  and ice fraction
    ! ocnwgt and icewgt are the "normal" fractions
    ! ocnwgt1, icewgt1, and wgtp01 are the fractions that switch between atm and mediator fluxes
    ! ocnwgt1+icewgt1+wgtp01 = 1.0 always
    ! wgtp01 and wgtm01 are the same just one is +1 and the other is -1 to change sign depending on the ice fraction.
    !   wgtp01 = 1 and wgtm01 = -1 when ice fraction = 0
    !   wgtp01 = 0 and wgtm01 =  0 when ice fraction > 0

    allocate(ocnwgt1(lsize))
    allocate(icewgt1(lsize))
    allocate(wgtp01(lsize))
    allocate(wgtm01(lsize))
    allocate(customwgt(lsize))

    do n = 1,lsize
       if (ifrac(n) <= 0._R8) then
          ! ice fraction is 0
          ocnwgt1(n) =  1.0_R8
          icewgt1(n) =  0.0_R8
          wgtp01(n)  =  0.0_R8
          wgtm01(n)  =  0.0_R8
       else
          ! ice fraction is > 0
          ocnwgt1(n) = ofrac(n)
          icewgt1(n) = ifrac(n)
          wgtp01(n)  = 0.0_R8
          wgtm01(n)  = 0.0_R8
       end if

       ! check wgts do add to 1 as expected
       ! TODO: check if this condition is still required
       if(ofrac(n)+ifrac(n) /= 0._R8)then
         if ( abs( ofrac(n) + ifrac(n) - 1.0_R8) > 1.0e-12 .or. &
              abs( ocnwgt1(n) + icewgt1(n) + wgtp01(n) - 1.0_R8) > 1.0e-12 .or. &
              abs( ocnwgt1(n) + icewgt1(n) - wgtm01(n) - 1.0_R8) > 1.0e-12) then

           write(6,100)trim(subname)//'ERROR: n, ofrac, ifrac, sum',&
                n,ofrac(n),ifrac(n),ofrac(n)+ifrac(n)
           write(6,101)trim(subname)//'ERROR: n, ocnwgt1, icewgt1, wgtp01, sum ', &
                n,ocnwgt1(n),icewgt1(n),wgtp01(n),ocnwgt1(n)+icewgt1(n)+wgtp01(n)
           write(6,101)trim(subname)//'ERROR: n, ocnwgt1, icewgt1, -wgtm01, sum ', &
                n,ocnwgt1(n),icewgt1(n),-wgtp01(n),ocnwgt1(n)+icewgt1(n)-wgtm01(n)
100        format(a,i8,2x,3(d20.13,2x))
101        format(a,i8,2x,4(d20.13,2x))

           call ESMF_LogWrite(trim(subname)//": ERROR atm + ice fracs inconsistent", &
                ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
           rc = ESMF_FAILURE
           return
         endif
       endif
    end do

    customwgt(:) = wgtm01(:) / const_lhvap
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_evap', &
         FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_evap', wgtA=ocnwgt1, &
         FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lat' , wgtB=customwgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_sen',    &
         FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_sen ', wgtA=ocnwgt1, &
         FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_sen' , wgtB=wgtm01, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_taux',  &
         FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_taux ', wgtA=ocnwgt1, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_taux' , wgtB=icewgt1, &
         FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_taux' , wgtC=wgtm01, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_tauy',  &
         FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_tauy ', wgtA=ocnwgt1, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_tauy' , wgtB=icewgt1, &
         FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_tauy' , wgtC=wgtm01, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! If there is no ice on the ocn gridcell (ocnwgt1=0) - sum Faxa_lwdn and Faxa_lwup
    ! If there is ice on the ocn gridcell -  merge Faox_lwup and Faxa_lwdn and ignore Faxa_lwup
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_lwnet', &
         FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_lwup ', wgtA=ocnwgt1, &
         FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lwdn' , wgtB=ocnwgt1, &
         FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_lwnet', wgtC=wgtp01, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_rain' , &
         FBInA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_rain' , wgtA=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_snow' , &
         FBInA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_snow' , wgtA=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! netsw_for_ocn = [downsw_from_atm*(1-ice_fraction)*(1-ocn_albedo)] + [pensw_from_ice*(ice_fraction)]
    customwgt(:) = ofrac(:) * (1.0 - 0.06)
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_vdr', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swvdr'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_vdr', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_vdf', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swvdf'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_vdf', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_idr', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swndr'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_idr', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_swnet_idf', &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_swndf'    , wgtA=customwgt, &
         FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_swpen_idf', wgtB=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(ocnwgt1)
    deallocate(icewgt1)
    deallocate(wgtp01)
    deallocate(wgtm01)
    deallocate(customwgt)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_custom_nemsdata

end module med_phases_prep_ocn_mod
