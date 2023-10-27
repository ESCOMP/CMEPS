module med_phases_cdeps_mod

  use ESMF, only: ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
  use ESMF, only: ESMF_Mesh
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_LogWrite
  use ESMF, only: ESMF_Field, ESMF_FieldGet
  use ESMF, only: ESMF_FieldBundleGet
  use ESMF, only: ESMF_StateIsCreated
  use ESMF, only: ESMF_GridCompGetInternalState
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGMSG_INFO

  use med_internalstate_mod, only: InternalState
  use med_internalstate_mod, only: compname, compatm, compocn
  use perf_mod             , only: t_startf, t_stopf
  use med_kind_mod         , only: cl => shr_kind_cl
  use med_kind_mod         , only: r8 => shr_kind_r8
  use med_constants_mod    , only: dbug_flag => med_constants_dbug_flag
  use med_utils_mod        , only: chkerr => med_utils_ChkErr
  use med_methods_mod      , only: med_methods_FB_FldChk 
  use med_methods_mod      , only: med_methods_FB_getFieldN
  use med_methods_mod      , only: FB_init_pointer => med_methods_FB_Init_pointer

  use dshr_strdata_mod     , only: shr_strdata_type
  use dshr_strdata_mod     , only: shr_strdata_init_from_inline
  use dshr_strdata_mod     , only: shr_strdata_advance

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_phases_cdeps_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type config
     integer :: year_first
     integer :: year_last
     integer :: year_align
     integer :: offset
     real(r8) :: dtlimit
     character(len=cl) :: mesh_filename
     character(len=cl), allocatable :: data_filename(:)
     character(len=cl), allocatable :: fld_list(:)
     character(len=cl), allocatable :: fld_list_model(:)
     character(len=cl) :: mapalgo
     character(len=cl) :: taxmode
     character(len=cl) :: tintalgo
     character(len=cl) :: name
  end type config

  type(config), allocatable           :: stream(:) ! stream configuration
  type(shr_strdata_type), allocatable :: sdat(:)   ! input data stream

  character(*),parameter :: u_FILE_u = __FILE__

!============================================================================
contains
!============================================================================

  subroutine med_phases_cdeps_run(gcomp, rc)

    !------------------------------------------------------------------------
    ! Use CDEPS inline capability to read in data
    !------------------------------------------------------------------------

    use ESMF, only : ESMF_GridComp

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: currTime
    type(ESMF_Mesh)             :: meshdst
    type(ESMF_Field)            :: flddst
    integer                     :: n1, n2, item, localPet
    integer                     :: curr_ymd, sec
    integer                     :: year, month, day, hour, minute, second
    logical, save               :: first_time = .true.
    character(len=*), parameter :: subname = '(med_phases_cdeps_run)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Query component
    call ESMF_GridCompGet(gcomp, clock=clock, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialize cdeps inline
    if (first_time) then 
       ! Set components in both side
       ! TODO: This needs to be dynamic and read from hconfig file
       n1 = compocn
       n2 = compatm

       ! Allocate data structures
       ! TODO: The number of stream will come from config file
       if (.not. allocated(sdat)) allocate(sdat(1))
       if (.not. allocated(stream)) allocate(stream(1))

       ! Check coupling direction
       if (n1 /= n2) then
          if (is_local%wrap%med_coupling_active(n1,n2)) then
             ! Get destination field
             call med_methods_FB_getFieldN(is_local%wrap%FBImp(n2,n2), 1, flddst, rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Get destination field mesh
             call ESMF_FieldGet(flddst, mesh=meshdst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Initialize cdeps inline
             call shr_strdata_init_from_inline(sdat(1), my_task = localPet, logunit = 6, &
                compname = trim(compname(n2)), &
                model_clock = clock, model_mesh = meshdst, &
                stream_meshfile     = 'INPUT_DATA/ESMFmesh.nc', &
                stream_filenames    = (/ 'INPUT_DATA/tsfc_fv3grid_202318612_sub.nc' /), &
                stream_yearFirst    = 2023, &
                stream_yearLast     = 2023,               &
                stream_yearAlign    = 2023,              &
                stream_fldlistFile  = (/ 'twsfc' /),                 &
                stream_fldListModel = (/ 'twsfc' /),                 &
                stream_lev_dimname  = 'null',                              &
                stream_mapalgo      = 'bilinear',                          &
                stream_offset       = 0,                                   &
                stream_taxmode      = 'cycle',                             &
                stream_dtlimit      = 1.5d0,                           &
                stream_tintalgo     = 'linear',                            &
                stream_name         = 'fvcom great lakes',         &
                rc                  = rc) 

             ! Create FB to store data
             if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n2), rc=rc)) then
                call FB_init_pointer(is_local%wrap%NStateExp(n2), is_local%wrap%FBExpInline(n2), &
                   is_local%wrap%flds_scalar_name, name='FBExpInline'//trim(compname(n2)), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
       end if

       ! Set flag to false
       first_time = .false.
    end if

    ! Get current time
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Query current time
    call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, h=hour, m=minute, s=second, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    curr_ymd = abs(year)*10000+month*100+day
    sec = hour*3600+minute*60+second

    ! Run inline cdeps and read data
    n1 = compocn
    n2 = compatm

    if (n1 /= n2) then
       if (is_local%wrap%med_coupling_active(n1,n2)) then
          call shr_strdata_advance(sdat(1), ymd=curr_ymd, tod=sec, logunit=6, istr=trim(compname(n2)), rc=rc) 
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Loop over fields provided by CDEPS inline and add it to FB
          do item = 1, 1 !size(config%stream_fldListFile)
             ! Get field
             !call ESMF_FieldBundleGet(sdat(1)%pstrm(1)%fldbun_model, fieldName=trim(config%stream_fldListFile(item)), field=flddst, rc=rc)
             call ESMF_FieldBundleGet(sdat(1)%pstrm(1)%fldbun_model, fieldName='So_t', field=flddst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Check FB for field
             !if (med_methods_FB_FldChk(is_local%wrap%FBExpInline(n2), trim(config%stream_fldListFile(item)))) then
             !   
             !end if

          end do
       end if
    end if      

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_cdeps_run

  !==========================================================================

  subroutine read_config()

    !------------------------------------------------------------------------
    ! Read YAML based Hconfig file
    !------------------------------------------------------------------------



  end subroutine read_config

end module med_phases_cdeps_mod
