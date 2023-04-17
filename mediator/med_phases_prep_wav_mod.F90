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
!  use esmFlds               , only : med_fldList_GetfldListTo
!  use med_internalstate_mod , only : compwav
  use esmFlds               , only : med_fldList_GetfldListTo, med_fldlist_type
  use med_internalstate_mod , only : compwav, compocn, compatm, compice, coupling_mode
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
    type(med_fldlist_type), pointer :: fldList
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
    fldList => med_fldList_GetfldListTo(compwav)
!PSH end
    ! auto merges to wav
!PSH begin
!    call med_merge_auto(&
!         is_local%wrap%med_coupling_active(:,compwav), &
!         is_local%wrap%FBExp(compwav), &
!         is_local%wrap%FBFrac(compwav), &
!         is_local%wrap%FBImp(:,compwav), &
!         med_fldList_GetfldListTo(compwav), rc=rc)
       call med_merge_auto(&
            is_local%wrap%med_coupling_active(:,compwav), &
            is_local%wrap%FBExp(compwav), &
            is_local%wrap%FBFrac(compwav), &
            is_local%wrap%FBImp(:,compwav), &
            fldList, &
            FBMed1=is_local%wrap%FBMed_aoflux_o, rc=rc)
!PSH end
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
end module med_phases_prep_wav_mod
