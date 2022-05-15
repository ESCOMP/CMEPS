  module ufs_io_mod

  use ESMF,                  only : operator(-)
  use ESMF,                  only : ESMF_VM, ESMF_VMGet, ESMF_VMGetCurrent, ESMF_LogWrite
  use ESMF,                  only : ESMF_GridComp, ESMF_GridCompGet, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF,                  only : ESMF_Field, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
  use ESMF,                  only : ESMF_Grid, ESMF_Decomp_Flag, ESMF_DECOMP_SYMMEDGEMAX
  use ESMF,                  only : ESMF_GridCreateMosaic, ESMF_INDEX_GLOBAL, ESMF_TYPEKIND_R8
  use ESMF,                  only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
  use ESMF,                  only : ESMF_GridCompGetInternalState, ESMF_KIND_R8
  use ESMF,                  only : ESMF_ArraySpec, ESMF_ArraySpecSet, ESMF_MESHLOC_ELEMENT
  use ESMF,                  only : ESMF_FieldCreate, ESMF_FieldGet, ESMF_FieldDestroy
  use ESMF,                  only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated
  use ESMF,                  only : ESMF_MeshGet, ESMF_FieldRegridStore, ESMF_FieldRedist
  use ESMF,                  only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF,                  only : ESMF_FieldWriteVTK, ESMF_VMAllFullReduce, ESMF_REDUCE_SUM
  use ESMF,                  only : ESMF_Calendar, ESMF_Clock, ESMF_ClockGet
  use ESMF,                  only : ESMF_ClockGetNextTime, ESMF_TimeIntervalGet
  use ESMF,                  only : ESMF_Time, ESMF_TimeGet, ESMF_TimeInterval
  use ESMF,                  only : ESMF_FieldBundleIsCreated
  use NUOPC,                 only : NUOPC_CompAttributeGet
  use NUOPC_Mediator,        only : NUOPC_MediatorGet

  use fms_mod,               only : fms_init
  use fms2_io_mod,           only : open_file, FmsNetcdfFile_t
  use mosaic2_mod,           only : get_mosaic_ntiles, get_mosaic_grid_sizes
  use mosaic2_mod,           only : get_mosaic_contact, get_mosaic_ncontacts
  use mpp_mod,               only : mpp_pe, mpp_root_pe, mpp_error, FATAL
  use mpp_domains_mod,       only : mpp_get_compute_domain
  use mpp_domains_mod,       only : mpp_domains_init, mpp_define_mosaic, domain2d
  use mpp_io_mod,            only : MPP_RDONLY, MPP_NETCDF, MPP_SINGLE, MPP_MULTI
  use mpp_io_mod,            only : mpp_get_info, mpp_get_fields, mpp_get_atts
  use mpp_io_mod,            only : mpp_open, mpp_read, fieldtype

  use med_kind_mod,          only : r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use med_utils_mod,         only : chkerr => med_utils_chkerr
  use med_constants_mod,     only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod, only : InternalState, mastertask, logunit
  use med_io_mod,            only : med_io_write, med_io_wopen, med_io_enddef, med_io_read
  use med_io_mod,            only : med_io_close, med_io_write_time, med_io_define_time
  use med_io_mod,            only : med_io_date2yyyymmdd, med_io_sec2hms, med_io_ymd2date
  use ufs_const_mod,         only : shr_const_cday
  use med_methods_mod,       only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod,       only : fldbun_diagnose => med_methods_FB_diagnose
  use med_methods_mod,       only : FB_fldchk => med_methods_FB_FldChk
  use med_methods_mod,       only : FB_getfldptr => med_methods_FB_GetFldPtr

  use MED_data,              only : physics

  implicit none

  private ! default private

  public read_initial
  public read_restart
  public write_restart

  type domain_type
     type(ESMF_Grid)        :: grid          ! ESMF grid object from mosaic file
     type(ESMF_RouteHandle) :: rh            ! ESMF route handle object to transfer data from grid to mesh
     type(domain2d)         :: mosaic_domain ! domain object created by FMS
     integer                :: layout(2)     ! layout for domain decomposition
     integer, allocatable   :: nit(:)        ! size of tile in i direction
     integer, allocatable   :: njt(:)        ! size of tile in j direction
     integer                :: ntiles        ! number of tiles in case of having CS grid
     integer                :: ncontacts     ! number of contacts in case of having CS grid
     integer, allocatable   :: tile1(:)      ! list of tile numbers in tile 1 of each contact
     integer, allocatable   :: tile2(:)      ! list of tile numbers in tile 2 of each contact
     integer, allocatable   :: istart1(:)    ! list of starting i-index in tile 1 of each contact
     integer, allocatable   :: iend1(:)      ! list of ending i-index in tile 1 of each contact
     integer, allocatable   :: jstart1(:)    ! list of starting j-index in tile 1 of each contact
     integer, allocatable   :: jend1(:)      ! list of ending j-index in tile 1 of each contact
     integer, allocatable   :: istart2(:)    ! list of starting i-index in tile 2 of each contact
     integer, allocatable   :: iend2(:)      ! list of ending i-index in tile 2 of each contact
     integer, allocatable   :: jstart2(:)    ! list of starting j-index in tile 2 of each contact
     integer, allocatable   :: jend2(:)      ! list of ending j-index in tile 2 of each contact
  end type domain_type

  type(ESMF_FieldBundle), save :: FBrst
  character(cs) :: prefix = 'ccpp'
  integer       :: file_ind = 10
  character(cl) :: case_name = 'unset'  ! case name

  character(*), parameter :: modName = "(ufs_io)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine read_initial(gcomp, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer, intent(inout)           :: rc

    ! local variables
    type(domain_type)                :: domain  
    type(ESMF_Field)                 :: field
    real(ESMF_KIND_R8), pointer      :: ptr(:,:,:)
    character(len=cl)                :: filename
    character(len=*), parameter      :: subname = trim(modName)//': (read_initial) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Create domain
    ! ---------------------

    call create_fms_domain(gcomp, domain, rc)

    ! ---------------------
    ! Create grid 
    ! ---------------------

    call create_grid(domain, rc)

    !----------------------
    ! Set file name for initial conditions
    !----------------------

    ! TODO: make file name configurable
    filename = 'INPUT/sfc_data.tile'
    call ESMF_LogWrite(subname//' read initial conditions from '//trim(filename)//'*', ESMF_LOGMSG_INFO)

    !----------------------
    ! Read surface friction velocity 
    !----------------------

    call read_tiled_file(gcomp, filename, 'uustar', domain, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    physics%sfcprop%uustar(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read surface roughness length
    !----------------------

    call read_tiled_file(gcomp, filename, 'zorl', domain, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    physics%sfcprop%zorl(:) = ptr(:,1,1)
    physics%sfcprop%zorlw(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read sea surface temperature, composite
    !----------------------

    call read_tiled_file(gcomp, filename, 'tsea', domain, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    physics%sfcprop%tsfco(:) = ptr(:,1,1)
    physics%sfcprop%tsfc(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read precipitation
    !----------------------

    call read_tiled_file(gcomp, filename, 'tprcp', domain, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    physics%sfcprop%tprcp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine read_initial

  !===============================================================================
  subroutine read_restart(gcomp, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp       ! gridded component
    integer, intent(inout)          :: rc          ! return code

    ! local variables
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    type(ESMF_Clock)  :: mclock
    type(ESMF_Time)   :: currtime
    type(ESMF_TimeInterval) :: timeStep 
    type(InternalState) :: is_local
    integer           :: n, yr, mon, day, sec
    real(r8), pointer :: ptr(:)
    logical           :: isPresent, isSet
    character(len=cl) :: cvalue
    character(len=cl) :: rest_file
    character(len=cl) :: currtime_str
    character(len=cs), allocatable :: flds(:)
    character(len=*), parameter :: subname=trim(modName)//': (read_restart) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Query VM 
    !----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Set restart file name
    !----------------------

    if (trim(case_name) == 'unset') then
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ccpp_restart_file', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       rest_file = trim(cvalue)
    else
       call NUOPC_MediatorGet(gcomp, mediatorClock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet(mclock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
       call ESMF_TimeGet(currTime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(currtime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       rest_file = trim(case_name)//'.cpl.ccpp.'//trim(currtime_str)//'.nc'
    end if

    !----------------------
    ! Now read in the restart file
    !----------------------

    if (mastertask) then
       write(logunit,'(a)') 'Reading CCPP restart file: '//trim(rest_file)
    end if

    ! create FB
    FBrst = ESMF_FieldBundleCreate(rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add fields
    allocate(flds(12))
    flds = (/ 'zorl  ', &
              'uustar', &
              'qss   ' /)
    do n = 1,size(flds)
       field = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
               name=trim(flds(n)), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       call ESMF_FieldGet(field, farrayptr=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0.0_r8
       nullify(ptr)
       call ESMF_FieldBundleAdd(FBrst, (/field/), rc=rc)
    end do 

    ! read file to FB
    call med_io_read(rest_file, vm, FBrst,  pre=trim(prefix), rc=rc)     
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//' diagnose at '//trim(currtime_str), ESMF_LOGMSG_INFO)
       call fldbun_diagnose(FBrst, string=trim(subname)//' CCPP FBrst ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Fill internal data structures
    !----------------------

    do n = 1,size(flds)
       if (FB_FldChk(FBrst, trim(flds(n)), rc=rc)) then              
          call FB_getfldptr(FBrst, trim(flds(n)), ptr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (mastertask) write(logunit,'(a)') 'Reading: '//trim(flds(n))
          if (trim(flds(n)) == 'zorl'  ) physics%sfcprop%zorl(:)  = ptr(:)
          if (trim(flds(n)) == 'uustar') physics%sfcprop%uustar(:)= ptr(:)
          if (trim(flds(n)) == 'qss'   ) physics%sfcprop%qss(:)   = ptr(:)

          nullify(ptr)
       end if
    end do
    deallocate(flds)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_restart

  !===============================================================================
  subroutine create_fms_domain(gcomp, domain, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(domain_type), intent(inout) :: domain
    integer, intent(inout)           :: rc

    ! local variables
    type(ESMF_VM)                    :: vm
    type(FmsNetcdfFile_t)            :: mosaic_fileobj
    integer                          :: mpicomm
    integer                          :: n, ntiles
    integer                          :: halo = 0
    integer                          :: global_indices(4,6)
    integer                          :: layout2d(2,6)
    integer, allocatable             :: pe_start(:), pe_end(:)
    character(len=cl)                :: msg, mosaic_file
    character(len=*), parameter      :: subname = trim(modName)//': (create_mosaic) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Initialize FMS
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm=vm, mpiCommunicator=mpicomm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fms_init(mpicomm)

    ! ---------------------
    ! Open mosaic file and query some information 
    ! ---------------------

    ! TODO: make mosaic file name configurable
    mosaic_file = 'INPUT/C96_mosaic.nc'

    if (.not. open_file(mosaic_fileobj, trim(mosaic_file), 'read')) then
       call ESMF_LogWrite(trim(subname)//'error in opening file '//trim(mosaic_file), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! query number of tiles
    domain%ntiles = get_mosaic_ntiles(mosaic_fileobj)

    ! query domain sizes for each tile
    if (.not. allocated(domain%nit)) allocate(domain%nit(domain%ntiles))
    if (.not. allocated(domain%njt)) allocate(domain%njt(domain%ntiles))
    call get_mosaic_grid_sizes(mosaic_fileobj, domain%nit, domain%njt)

    ! query number of contacts
    domain%ncontacts = get_mosaic_ncontacts(mosaic_fileobj)

    ! allocate required arrays to create FMS domain from mosaic file
    if (.not. allocated(domain%tile1))   allocate(domain%tile1(domain%ncontacts))
    if (.not. allocated(domain%tile2))   allocate(domain%tile2(domain%ncontacts))
    if (.not. allocated(domain%istart1)) allocate(domain%istart1(domain%ncontacts))
    if (.not. allocated(domain%iend1))   allocate(domain%iend1(domain%ncontacts))
    if (.not. allocated(domain%jstart1)) allocate(domain%jstart1(domain%ncontacts))
    if (.not. allocated(domain%jend1))   allocate(domain%jend1(domain%ncontacts))
    if (.not. allocated(domain%istart2)) allocate(domain%istart2(domain%ncontacts))
    if (.not. allocated(domain%iend2))   allocate(domain%iend2(domain%ncontacts))
    if (.not. allocated(domain%jstart2)) allocate(domain%jstart2(domain%ncontacts))
    if (.not. allocated(domain%jend2))   allocate(domain%jend2(domain%ncontacts))

    ! query information about contacts
    call get_mosaic_contact(mosaic_fileobj, domain%tile1, domain%tile2, &
         domain%istart1, domain%iend1, domain%jstart1, domain%jend1, &
         domain%istart2, domain%iend2, domain%jstart2, domain%jend2)

    ! print out debug information
    if (dbug_flag > 5) then
       do n = 1, domain%ncontacts
          write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' : tile1, tile2 (', n ,') = ', domain%tile1(n), domain%tile2(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
          write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' : istart1, iend1, jstart1, jend1 (', n ,') = ', &
            domain%istart1(n), domain%iend1(n), domain%jstart1(n), domain%jend1(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
          write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' : istart2, iend2, jstart2, jend2 (', n ,') = ', &
            domain%istart2(n), domain%iend2(n), domain%jstart2(n), domain%jend2(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
       end do
    end if

    !----------------------
    ! Initialize domain 
    !----------------------

    call mpp_domains_init()

    !----------------------
    ! Set pe_start, pe_end 
    !----------------------

    ! TODO: make layout options configurable
    domain%layout(1) = 3
    domain%layout(2) = 8

    allocate(pe_start(domain%ntiles))
    allocate(pe_end(domain%ntiles))
    do n = 1, domain%ntiles
       pe_start(n) = mpp_root_pe()+(n-1)*domain%layout(1)*domain%layout(2)
       pe_end(n) = mpp_root_pe()+n*domain%layout(1)*domain%layout(2)-1
       if (dbug_flag > 5) then
          write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' pe_start, pe_end (', n ,') = ', pe_start(n), pe_end(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end if 
    enddo

    !----------------------
    ! Create FMS domain object
    !----------------------

    do n = 1, domain%ntiles
       layout2d(:,n) = domain%layout(:)
       global_indices(1,n) = 1
       global_indices(2,n) = domain%nit(n)
       global_indices(3,n) = 1
       global_indices(4,n) = domain%njt(n)
    enddo

    call mpp_define_mosaic(global_indices, layout2d, domain%mosaic_domain, &
         domain%ntiles, domain%ncontacts, domain%tile1, domain%tile2, &
         domain%istart1, domain%iend1, domain%jstart1, domain%jend1, &
         domain%istart2, domain%iend2, domain%jstart2, domain%jend2, &
         pe_start, pe_end, symmetry=.true., &
         whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
         name='atm domain')

    !----------------------
    ! Deallocate temporary arrays
    !----------------------

    deallocate(pe_start)
    deallocate(pe_end)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine create_fms_domain

  !===============================================================================
  subroutine create_grid(domain, rc)
    implicit none

    ! input/output variables
    type(domain_type), intent(inout) :: domain
    integer, intent(inout)           :: rc

    ! local variables
    type(ESMF_Decomp_Flag)           :: decompflagPTile(2,6)
    integer                          :: n
    integer                          :: decomptile(2,6)
    character(len=cl)                :: mosaic_file, input_dir
    character(len=*), parameter      :: subname = trim(modName)//': (create_grid) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! TODO: make mosaic file name and input folder configurable
    mosaic_file = 'INPUT/C96_mosaic.nc'
    input_dir = 'INPUT/'

    ! TODO: currently this is only tested with global application
    ! set decomposition
    do n = 1, domain%ntiles
       decomptile(1,n) = domain%layout(1)
       decomptile(2,n) = domain%layout(2)
       decompflagPTile(:,n) = (/ ESMF_DECOMP_SYMMEDGEMAX, ESMF_DECOMP_SYMMEDGEMAX /)
    end do

    ! create grid
    domain%grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file), &
       regDecompPTile=decomptile, tileFilePath=trim(input_dir), decompflagPTile=decompflagPTile, &
       staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
       indexflag=ESMF_INDEX_GLOBAL, name='input_grid', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine create_grid

  !===============================================================================
  subroutine read_tiled_file(gcomp, filename, varname, domain, field_dst, numrec, numlev, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)          :: gcomp
    character(len=*), intent(in)             :: filename
    character(len=*), intent(in)             :: varname 
    type(domain_type), intent(inout)         :: domain
    type(ESMF_Field), intent(inout)          :: field_dst
    integer, intent(in), optional            :: numrec
    integer, intent(in), optional            :: numlev
    integer, intent(inout), optional         :: rc

    ! local variables
    type(ESMF_Field)            :: field_src, field_tmp
    type(ESMF_ArraySpec)        :: arraySpec
    type(InternalState)         :: is_local
    type(fieldtype), allocatable:: vars(:) 
    integer                     :: funit, my_tile
    integer                     :: i, j, n, nt, nl
    integer                     :: isc, iec, jsc, jec
    integer                     :: ndim, nvar, natt, ntime
    logical                     :: not_found, is_root_pe
    real(ESMF_KIND_R8), pointer :: ptr(:), ptr3d(:,:,:)
    real(ESMF_KIND_R8), pointer :: ptr4d(:,:,:,:)
    real(r8), allocatable       :: rdata(:,:,:,:)
    character(len=cl)           :: cname, fname
    character(len=*), parameter :: subname=trim(modName)//': (read_tiled_file) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' reading '//trim(varname), ESMF_LOGMSG_INFO)

    !----------------------
    ! Get the internal state from the mediator component
    !----------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Define required variables
    !----------------------

    if (present(numrec)) then
       nt = numrec
    else
       nt = 1
    end if

    if (present(numlev)) then
       nl = numlev
    else
       nl = 1
    end if

    my_tile = int(mpp_pe()/(domain%layout(1)*domain%layout(2)))+1

    is_root_pe = .false.
    if (mpp_pe() == (my_tile-1)*(domain%layout(1)*domain%layout(2))) is_root_pe = .true.

    !----------------------
    ! Open file and query file attributes
    !----------------------
   
    write(cname, fmt='(A,I1,A)') trim(filename), my_tile, '.nc' 
    call mpp_open(funit, trim(cname), action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE, is_root_pe=is_root_pe)
    call mpp_get_info(funit, ndim, nvar, natt, ntime)
    allocate(vars(nvar))
    call mpp_get_fields(funit, vars(:))

    !----------------------
    ! Find and read requested variable
    !----------------------

    not_found = .true.
    do n = 1, nvar
       ! get variable name
       call mpp_get_atts(vars(n), name=cname)

       ! check variable name
       if (trim(cname) == trim(varname)) then
          ! get array bounds or domain
          call mpp_get_compute_domain(domain%mosaic_domain, isc, iec, jsc, jec)

          ! allocate data array and set initial value
          allocate(rdata(isc:iec,jsc:jec,nl,nt))
          rdata(:,:,:,:) = 0.0_r8

          ! read data
          do i = 1, nt
             call mpp_read(funit, vars(n), domain%mosaic_domain, rdata, 1)
          end do

          ! set missing values to zero
          where (rdata == 1.0e20)
             rdata(:,:,:,:) = 0.0_r8
          end where
       end if

       not_found = .false.
    end do

    if (not_found) then 
       call mpp_error(FATAL, 'File being read is not the expected one. '//trim(varname)//' is not found.')
    end if

    !----------------------
    ! Move data from grid to mesh
    !----------------------

    ! set type and rank for ESMF arrayspec
    call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create source field
    field_src = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
       indexflag=ESMF_INDEX_GLOBAL, ungriddedLBound=(/1,1/), ungriddedUBound=(/nl,nt/), &
       gridToFieldMap=(/1,2/), name=trim(varname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and fill it
    call ESMF_FieldGet(field_src, localDe=0, farrayPtr=ptr4d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ptr4d(:,:,:,:) = rdata(:,:,:,:)
    nullify(ptr4d)
    if (allocated(rdata)) deallocate(rdata)

    ! create destination field
    field_dst = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, & 
       name=trim(varname), meshloc=ESMF_MESHLOC_ELEMENT,  ungriddedLbound=(/1,1/), &
       ungriddedUbound=(/nl,nt/), gridToFieldMap=(/1/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create routehandle from grid to mesh
    if (.not. ESMF_RouteHandleIsCreated(domain%rh, rc=rc)) then
       call ESMF_FieldRegridStore(field_src, field_dst, routehandle=domain%rh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! redist field from ESMF Grid to Mesh
    call ESMF_FieldRedist(field_src, field_dst, domain%rh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! clean memory
    call ESMF_FieldDestroy(field_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return 

    !----------------------
    ! Output result field for debugging purpose
    !----------------------

    if (dbug_flag > 5) then
       ! TODO: ESMF_FieldWriteVTK() call does not support ungridded dimension
       ! The workaround is implemented in here but it would be nice to extend
       ! ESMF_FieldWriteVTK() call to handle it.
       field_tmp = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
          name=trim(varname), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(field_tmp, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(field_dst, localDe=0, farrayPtr=ptr3d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! write to different file along ungridded dimension
       do i = 1, nl
          do j = 1, nt
            ptr(:) = ptr3d(:,i,j)
            write(fname, fmt='(A,I2.2,A,I2.2)') trim(varname)//'_lev', i, '_time', j
            call ESMF_FieldWriteVTK(field_tmp, trim(fname), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end do
       end do

       ! clean memory
       nullify(ptr)
       nullify(ptr3d)
       call ESMF_FieldDestroy(field_tmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine read_tiled_file

  !===============================================================================
  subroutine write_restart(gcomp, restart_freq, rc)
    implicit none

    ! input/output variableswrite_restart
    type(ESMF_GridComp), intent(in) :: gcomp        ! gridded component
    integer, intent(in)             :: restart_freq ! restart interval in hours 
    integer, intent(inout)          :: rc           ! return code

    ! local variables
    type(ESMF_VM) :: vm
    type(ESMF_Field) :: field
    type(ESMF_Clock) :: mclock
    type(ESMF_Calendar) :: calendar
    type(ESMF_Time) :: currtime, starttime, nexttime
    type(ESMF_TimeInterval) :: timediff(2)
    type(InternalState) :: is_local
    integer           :: yr, mon, day, sec
    integer           :: m, ns, start_ymd
    character(cl)     :: time_units
    real(r8)          :: time_val
    real(r8)          :: time_bnds(2)
    real(r8), pointer :: ptr(:)
    logical :: whead(2) = (/.true. , .false./)
    logical :: wdata(2) = (/.false., .true. /)
    logical           :: isPresent, isSet
    character(len=cl) :: tmpstr
    character(len=cl) :: rest_file
    character(len=cl) :: nexttime_str
    integer, save     :: ns_total
    logical, save     :: first_call = .true.
    character(len=*), parameter :: subname=trim(modName)//': (write_restart) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Determine clock, starttime, currtime and nexttime
    !----------------------

    call NUOPC_MediatorGet(gcomp, mediatorClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return 
    call ESMF_ClockGet(mclock, currtime=currtime, starttime=starttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(mclock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Determine time units
    !----------------------

    call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_ymd2date(yr, mon, day, start_ymd)
    time_units = 'days since '//trim(med_io_date2yyyymmdd(start_ymd))//' '//med_io_sec2hms(sec, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Determine restart file name
    !----------------------

    if (trim(case_name) == 'unset') then
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    rest_file = trim(case_name)//'.cpl.ccpp.'//trim(nexttime_str)//'.nc'

    ! return if it is not time to write restart
    if (restart_freq < 0) return
    if (mod(sec, restart_freq) /= 0) return

    !----------------------
    ! Create restart file
    !----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_wopen(trim(rest_file), vm, clobber=.true., file_ind=file_ind)
    if (mastertask) then
       write(logunit,'(a)') 'CCPP restart file is created: '//trim(rest_file)
    end if

    !----------------------
    ! Define time dimension
    !----------------------

    timediff(1) = nexttime - starttime
    call ESMF_TimeIntervalGet(timediff(1), d=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    time_val = day + sec/real(shr_const_cday,r8)
    time_bnds(1) = time_val
    time_bnds(2) = time_val

    !----------------------
    ! Create FB and add fields to it
    !----------------------

    if (first_call) then
       ! create FB
       FBrst = ESMF_FieldBundleCreate(rc=rc)

       ! get total element count 
       call ESMF_MeshGet(is_local%wrap%aoflux_mesh, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMAllFullReduce(vm, (/ns/), ns_total, 1, ESMF_REDUCE_SUM, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! surface roughness length in cm
       field = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
               name='zorl', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, farrayptr=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%zorl(:)
       call ESMF_FieldBundleAdd(FBrst, (/field/), rc=rc)

       ! boundary layer parameter
       field = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
               name='uustar', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, farrayptr=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%uustar(:)
       nullify(ptr)
       call ESMF_FieldBundleAdd(FBrst, (/field/), rc=rc)

       ! surface air saturation specific humidity (kg/kg)
       field = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
               name='qss', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, farrayptr=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%qss(:)
       nullify(ptr)
       call ESMF_FieldBundleAdd(FBrst, (/field/), rc=rc)
    else
       call fldbun_getdata1d(FBrst, 'zorl', ptr, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%zorl(:)
       nullify(ptr)

       call fldbun_getdata1d(FBrst, 'uustar', ptr, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%uustar(:)
       nullify(ptr)

       call fldbun_getdata1d(FBrst, 'qss', ptr, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = physics%sfcprop%qss(:)
       nullify(ptr)
    end if

    ! diagnose
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//' diagnose at '//trim(nexttime_str), ESMF_LOGMSG_INFO)
       call fldbun_diagnose(FBrst, string=trim(subname)//' CCPP FBrst ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! debug
    

    !----------------------
    ! Write data
    !----------------------

    ! loop over whead/wdata phases
    do m = 1, 2
       if (m == 2) then
          call med_io_enddef(rest_file, file_ind=file_ind)
       end if

       ! write time values
       if (whead(m)) then
          call ESMF_ClockGet(mclock, calendar=calendar, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_define_time(time_units, calendar, file_ind=file_ind, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call med_io_write_time(time_val, time_bnds, nt=1, file_ind=file_ind, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! write data
       call med_io_write(rest_file, FBrst, whead(m), wdata(m), ns_total, 1, nt=1, pre=trim(prefix), file_ind=file_ind, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !----------------------
    ! Close file
    !----------------------

    call med_io_close(rest_file, vm, file_ind=file_ind, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       write(logunit,'(a)') 'CCPP restart file is closed: '//trim(rest_file)
    end if

  end subroutine write_restart

  end module ufs_io_mod
