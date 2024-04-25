module med_phases_cdeps_mod

  use ESMF, only: ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
  use ESMF, only: ESMF_Mesh
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_LogWrite
  use ESMF, only: ESMF_Field, ESMF_FieldGet
  use ESMF, only: ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
  use ESMF, only: ESMF_FieldBundleCreate
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGMSG_INFO

  use med_internalstate_mod, only: InternalState
  use med_internalstate_mod, only: logunit, maintask
  use med_internalstate_mod, only: ncomps, compname, compatm, compocn
  use perf_mod             , only: t_startf, t_stopf
  use med_kind_mod         , only: cl => shr_kind_cl
  use med_kind_mod         , only: r8 => shr_kind_r8
  use med_constants_mod    , only: dbug_flag => med_constants_dbug_flag
  use med_utils_mod        , only: chkerr => med_utils_ChkErr
  use med_methods_mod      , only: FB_FldChk => med_methods_FB_FldChk 
  use med_methods_mod      , only: FB_getFieldN => med_methods_FB_getFieldN
  use med_methods_mod      , only: FB_getNumflds => med_methods_FB_getNumflds 
  use med_methods_mod      , only: FB_init => med_methods_FB_Init
  use med_methods_mod      , only: FB_diagnose => med_methods_FB_diagnose
  use med_methods_mod      , only: FB_write => med_methods_FB_write
  use med_methods_mod      , only: FB_GetFldPtr => med_methods_FB_GetFldPtr

  use dshr_mod             , only: dshr_pio_init
  use dshr_strdata_mod     , only: shr_strdata_type
  use dshr_strdata_mod     , only: shr_strdata_init_from_inline
  use dshr_strdata_mod     , only: shr_strdata_advance
  use dshr_stream_mod      , only: shr_stream_init_from_esmfconfig

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

  type(config) :: stream ! stream configuration
  type(shr_strdata_type), allocatable :: sdat(:,:) ! input data stream

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
    type(InternalState)            :: is_local
    type(ESMF_Clock)               :: clock
    type(ESMF_Time)                :: currTime
    type(ESMF_Mesh)                :: meshdst
    type(ESMF_Field)               :: flddst
    integer                        :: i, j, k, l, nflds, streamid
    integer                        :: n1, n2, item, nstreams, localPet
    integer                        :: curr_ymd, sec
    integer                        :: year, month, day, hour, minute, second
    logical                        :: found 
    logical, save                  :: first_time = .true.
    character(len=cl), allocatable :: fileList(:), varList(:,:)
    character(len=cl)              :: streamfilename, suffix, fldname
    type(shr_strdata_type)         :: sdat_config
    character(len=*), parameter    :: subname = '(med_phases_cdeps_run)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Query component
    call ESMF_GridCompGet(gcomp, clock=clock, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialize sdat streams
    if (.not. allocated(sdat)) allocate(sdat(ncomps,ncomps))
    sdat(:,:)%mainproc = (localPet == 0)

    ! Initialize cdeps inline
    if (first_time) then 
       ! Init PIO
       call dshr_pio_init(gcomp, sdat_config, logunit, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Read stream configuration file
       ! TODO: At this point it only suports ESMF config format (XML?)
       streamfilename = 'stream.config' 
       call shr_stream_init_from_esmfconfig(streamfilename, sdat_config%stream, logunit, &
          sdat_config%pio_subsystem, sdat_config%io_type, sdat_config%io_format, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get number of streams
       nstreams = size(sdat_config%stream)

       ! Loop over coupling directions and try to find field match in given streams
       do n1 = 1, ncomps
          do n2 = 1, ncomps
             ! Check for coupling direction and background fill
             if (n1 /= n2 .and. is_local%wrap%med_coupling_active(n1,n2) .and. is_local%wrap%med_data_active(n1,n2)) then
                ! Get number of fields
                call FB_getNumflds(is_local%wrap%FBImp(n1,n2), trim(subname), nflds, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                
                ! Loop over fields and try to find it in the given stream
                found = .false.
                do i = 1, nflds
                   ! Query destination field
                   call FB_getFieldN(is_local%wrap%FBImp(n1,n2), i, flddst, rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return

                   ! Query destination field name and its mesh
                   call ESMF_FieldGet(flddst, mesh=meshdst, name=fldname, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   if (maintask) write(logunit,'(a)') trim(subname)//": extracting destination mesh from "//trim(fldname)

                   ! Check if any field in FB in the given stream
                   ! NOTE: Single stream could provide multiple fields !!!
                   streamid = 0
                   do j = 1, nstreams
                      do k = 1, sdat_config%stream(j)%nvars
                         if (trim(sdat_config%stream(j)%varlist(k)%nameinmodel) == trim(fldname)) then
                            streamid = j
                         end if
                      end do
                   end do

                   ! If match is found and previously not initialized, then initialize cdeps inline for the stream
                   if (size(sdat(n1,n2)%stream) == 0 .and. streamid /= 0) then
                      ! Debug print
                      if (maintask) then
                         write(logunit,'(a,i3)') trim(subname)//": initialize stream ", streamid
                      end if

                      ! Allocate temporary variable to store file names in the stream
                      allocate(fileList(sdat_config%stream(streamid)%nfiles)) 
                      allocate(varList(sdat_config%stream(streamid)%nvars,2))
                      
                      ! Fill file abd variable lists with data
                      do l = 1, sdat_config%stream(streamid)%nfiles
                         fileList(l) = trim(sdat_config%stream(streamid)%file(l)%name) 
                         if (maintask) write(logunit,'(a,i2,2x,a)') trim(subname)//": file     ", l, trim(fileList(l))
                      end do
                      do l = 1, sdat_config%stream(streamid)%nvars
                         varList(l,1) = trim(sdat_config%stream(streamid)%varlist(l)%nameinfile)
                         varList(l,2) = trim(sdat_config%stream(streamid)%varlist(l)%nameinmodel)
                         if (maintask) write(logunit,'(a,i2,2x,a)') trim(subname)//": variable ", l, trim(varList(l,1))//" -> "//trim(varList(l,2))
                      end do

                      ! Set PIO related variables
                      sdat(n1,n2)%pio_subsystem => sdat_config%pio_subsystem
                      sdat(n1,n2)%io_type = sdat_config%io_type
                      sdat(n1,n2)%io_format = sdat_config%io_format

                      ! Init stream
                      call shr_strdata_init_from_inline(sdat(n1,n2), my_task=localPet, logunit=logunit, &
                         compname = 'cmeps', model_clock=clock, model_mesh=meshdst, &
                         stream_meshfile=trim(sdat_config%stream(streamid)%meshfile), &
                         stream_filenames=fileList, &
                         stream_yearFirst=sdat_config%stream(streamid)%yearFirst, &
                         stream_yearLast=sdat_config%stream(streamid)%yearLast, &
                         stream_yearAlign=sdat_config%stream(streamid)%yearAlign, &
                         stream_fldlistFile=varList(:,1), &
                         stream_fldListModel=varList(:,2), &
                         stream_lev_dimname=trim(sdat_config%stream(streamid)%lev_dimname), &
                         stream_mapalgo=trim(sdat_config%stream(streamid)%mapalgo), &
                         stream_offset=sdat_config%stream(streamid)%offset, &
                         stream_taxmode=trim(sdat_config%stream(streamid)%taxmode), &
                         stream_dtlimit=sdat_config%stream(streamid)%dtlimit, &
                         stream_tintalgo=trim(sdat_config%stream(streamid)%tInterpAlgo), &
                         stream_name=trim(compname(n1))//'_'//trim(compname(n2)), &
                         stream_src_mask=sdat_config%stream(streamid)%src_mask_val, &
                         stream_dst_mask=sdat_config%stream(streamid)%dst_mask_val, &
                         rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return

                      ! Print out source and destination mask used in the stream
                      if (maintask) write(logunit,'(a,2i2)') trim(subname)//": mask values src, dst ", &
                         sdat_config%stream(streamid)%src_mask_val, sdat_config%stream(streamid)%dst_mask_val

                      ! Remove temporary variables
                      deallocate(fileList)
                      deallocate(varList)

                      ! Set flag
                      found = .true.
                   end if
                end do ! nflds

                ! Create empty FB
                if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBData(n2), rc=rc) .and. found) then
                   is_local%wrap%FBData(n2) = ESMF_FieldBundleCreate(name="inline_"//trim(compname(n2)), rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
          end do ! n2
       end do ! n1

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

    ! Read data if stream initialized
    do n1 = 1, ncomps
       do n2 = 1, ncomps
          if (size(sdat(n1,n2)%stream) > 0) then
             ! Debug print
             if (maintask) then
                 write(logunit,'(a)') trim(subname)//": read stream "//trim(compname(n1))//" -> "//trim(compname(n2))  
             end if

             ! Read data
             call shr_strdata_advance(sdat(n1,n2), ymd=curr_ymd, tod=sec, logunit=logunit, &
                istr=trim(compname(n1))//'_'//trim(compname(n2)), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Check FB
             call FB_diagnose(sdat(n1,n2)%pstrm(1)%fldbun_model, &
                trim(subname)//':'//trim(compname(n1))//'_'//trim(compname(n2)), rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Point FB from internal one
             is_local%wrap%FBData(n2) = sdat(n1,n2)%pstrm(1)%fldbun_model
  
             ! Write FB for debugging
             if (dbug_flag > 10) then
                write(suffix, fmt='(i4,a1,i2.2,a1,i2.2,a1,i5.5)') year, '-', month, '-', day, '-', sec
                call FB_write(is_local%wrap%FBData(n2), suffix, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
       end do
    end do

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    call t_stopf('MED:'//subname)

  end subroutine med_phases_cdeps_run

end module med_phases_cdeps_mod
