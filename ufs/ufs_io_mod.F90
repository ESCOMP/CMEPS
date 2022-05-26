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
  use ESMF,                  only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated, ESMF_FieldRedist
  use ESMF,                  only : ESMF_MeshGet, ESMF_FieldRegrid, ESMF_FieldRegridStore
  use ESMF,                  only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF,                  only : ESMF_FieldWriteVTK, ESMF_VMAllFullReduce, ESMF_REDUCE_SUM
  use ESMF,                  only : ESMF_Mesh, ESMF_Calendar, ESMF_Clock, ESMF_ClockGet
  use ESMF,                  only : ESMF_ClockGetNextTime, ESMF_TimeIntervalGet
  use ESMF,                  only : ESMF_Time, ESMF_TimeGet, ESMF_TimeInterval
  use ESMF,                  only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
  use ESMF,                  only : ESMF_FieldBundleRemove, ESMF_FieldBundleDestroy
  use ESMF,                  only : ESMF_FieldWrite, ESMF_FieldBundleRead, ESMF_FieldBundleWrite
  use ESMF,                  only : ESMF_REGRIDMETHOD_CONSERVE_2ND, ESMF_MeshCreate
  use NUOPC,                 only : NUOPC_CompAttributeGet
  use NUOPC_Mediator,        only : NUOPC_MediatorGet

  use fms_mod,               only : fms_init
  use fms2_io_mod,           only : open_file, FmsNetcdfFile_t
  use mosaic2_mod,           only : get_mosaic_ntiles, get_mosaic_grid_sizes
  use mosaic2_mod,           only : get_mosaic_contact, get_mosaic_ncontacts
  use mpp_mod,               only : mpp_pe, mpp_root_pe, mpp_error, FATAL
  use mpp_domains_mod,       only : mpp_define_layout, mpp_get_compute_domain
  use mpp_domains_mod,       only : mpp_domains_init, mpp_define_mosaic, domain2d
  use mpp_io_mod,            only : MPP_RDONLY, MPP_NETCDF, MPP_SINGLE, MPP_MULTI
  use mpp_io_mod,            only : mpp_get_info, mpp_get_fields, mpp_get_atts
  use mpp_io_mod,            only : mpp_open, mpp_read, fieldtype

  use med_kind_mod,          only : r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use med_utils_mod,         only : chkerr => med_utils_chkerr
  use med_constants_mod,     only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod, only : InternalState, mastertask, logunit
  use med_internalstate_mod, only : compatm, compocn, mapconsf
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
     type(ESMF_Mesh)        :: mesh          ! ESMF mesh object from CS grid 
     type(ESMF_RouteHandle) :: rh            ! ESMF routehandle object to redist data from CS grid to mesh
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

  character(cl) :: case_name = 'unset'  ! case name

  character(*), parameter :: modName = "(ufs_io_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine read_initial(gcomp, ini_file, mosaic_file, input_dir, rc)

    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    character(len=cl), intent(in)    :: ini_file
    character(len=cl), intent(in)    :: mosaic_file
    character(len=cl), intent(in)    :: input_dir
    integer, intent(inout)           :: rc

    ! local variables
    type(domain_type)                :: domain
    type(InternalState)              :: is_local
    type(ESMF_RouteHandle)           :: rh
    type(ESMF_Field)                 :: lfield, field, field_dst
    real(ESMF_KIND_R8), pointer      :: ptr(:)
    integer                          :: n
    character(len=cs), allocatable   :: flds(:)
    character(len=*), parameter      :: subname = trim(modName)//': (read_initial) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create domain
    ! ---------------------

    call create_fms_domain(gcomp, domain, mosaic_file, rc)

    ! ---------------------
    ! Create grid 
    ! ---------------------

    call create_grid(gcomp, domain, mosaic_file, input_dir, rc)

    !----------------------
    ! Read data 
    !----------------------

    allocate(flds(2))
    flds = (/ 'zorl  ', &
              'uustar' /)
    do n = 1,size(flds)
       ! read from tiled file
       call read_tiled_file(gcomp, ini_file, trim(flds(n)), domain, field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create destination field
       field_dst = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
          name='uustar', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)        
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! map field
       if (is_local%wrap%aoflux_grid == 'ogrid' .or. is_local%wrap%aoflux_grid == 'xgrid') then
          ! create rh
          call ESMF_FieldRegridStore(field, field_dst, routehandle=rh, &
               regridmethod=ESMF_REGRIDMETHOD_CONSERVE_2ND, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! remap from atm to ocn/xgrid
          call ESMF_FieldRegrid(field, field_dst, rh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          ! do nothing, use source field
          field_dst = field
       end if

       ! debug
       if (dbug_flag > 5) then
          call ESMF_FieldWriteVTK(field_dst, 'ini_'//trim(flds(n))//'_'//trim(is_local%wrap%aoflux_grid), rc=rc)  
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! return pointer and fill variable
       call ESMF_FieldGet(field_dst, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit,'(a)') 'Reading: '//trim(flds(n))
       if (trim(flds(n)) == 'zorl'  ) physics%sfcprop%zorl(:)  = ptr(:)
       if (trim(flds(n)) == 'uustar') physics%sfcprop%uustar(:)= ptr(:)
       nullify(ptr)

       ! free memory
       call ESMF_FieldDestroy(field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! free memory
    if (allocated(flds)) deallocate(flds)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_initial

  !===============================================================================
  subroutine read_restart(gcomp, rst_file, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp    ! gridded component
    character(len=cl), intent(inout):: rst_file ! restart file
    integer, intent(inout)          :: rc       ! return code

    ! local variables
    type(ESMF_Field)  :: field, lfield
    type(ESMF_Clock)  :: mclock
    type(ESMF_Time)   :: currtime
    type(ESMF_TimeInterval) :: timeStep 
    type(ESMF_FieldBundle), save :: FBin
    type(InternalState) :: is_local
    integer           :: n, yr, mon, day, sec
    real(r8), pointer :: ptr(:)
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
    ! Set restart file name
    !----------------------

    if (trim(case_name) == 'unset') then
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (trim(rst_file) == 'unset') then
       call NUOPC_MediatorGet(gcomp, mediatorClock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet(mclock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
       call ESMF_TimeGet(currTime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(currtime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       rst_file = trim(case_name)//'.cpl.ccpp.'//trim(currtime_str)//'.nc'
    end if

    !----------------------
    ! Now read in the restart file
    !----------------------

    if (mastertask) then
       write(logunit,'(a)') 'Reading CCPP restart file: '//trim(rst_file)
    end if

    ! create FB
    FBin = ESMF_FieldBundleCreate(rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add fields
    allocate(flds(3))
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
       call ESMF_FieldBundleAdd(FBin, (/field/), rc=rc)
    end do 

    ! read file to FB
    call ESMF_FieldBundleRead(FBin, trim(rst_file), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! debug
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//' diagnose at '//trim(currtime_str), ESMF_LOGMSG_INFO)
       call fldbun_diagnose(FBin, string=trim(subname)//' CCPP FBin ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Fill internal data structures
    !----------------------

    do n = 1,size(flds)
       if (FB_FldChk(FBin, trim(flds(n)), rc=rc)) then              
          call FB_getfldptr(FBin, trim(flds(n)), ptr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (mastertask) write(logunit,'(a)') 'Reading: '//trim(flds(n))
          if (trim(flds(n)) == 'zorl'  ) physics%sfcprop%zorl(:)  = ptr(:)
          if (trim(flds(n)) == 'uustar') physics%sfcprop%uustar(:)= ptr(:)
          if (trim(flds(n)) == 'qss'   ) physics%sfcprop%qss(:)   = ptr(:)

          nullify(ptr)

          ! debug
          if (dbug_flag > 5) then
             call ESMF_FieldBundleGet(FBin, fieldName=trim(flds(n)), field=lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldWriteVTK(lfield, 'rst_'//trim(flds(n))//'_'//trim(is_local%wrap%aoflux_grid), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

    !----------------------
    ! Free memory
    !----------------------

    do n = 1,size(flds)
       if (FB_FldChk(FBin, trim(flds(n)), rc=rc)) then
          ! get field from FB
          call ESMF_FieldBundleGet(FBin, trim(flds(n)), field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          ! remove field from FB
          call ESMF_FieldBundleRemove(FBin, (/ trim(flds(n)) /), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          ! remove field
          call ESMF_FieldDestroy(field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do
    deallocate(flds)

    ! remove FB
    call ESMF_FieldBundleDestroy(FBin, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_restart

  !===============================================================================
  subroutine create_fms_domain(gcomp, domain, mosaic_file, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(domain_type), intent(inout) :: domain
    character(len=cl), intent(in)    :: mosaic_file
    integer, intent(inout)           :: rc

    ! local variables
    type(ESMF_VM)                    :: vm
    type(FmsNetcdfFile_t)            :: mosaic_fileobj
    integer                          :: mpicomm, npes_per_tile
    integer                          :: n, ntiles, npet
    integer                          :: halo = 0
    integer                          :: global_indices(4,6)
    integer                          :: layout2d(2,6)
    integer, allocatable             :: pe_start(:), pe_end(:)
    character(len=cl)                :: msg
    character(len=*), parameter      :: subname = trim(modName)//': (create_fms_domain) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Initialize FMS
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm=vm, mpiCommunicator=mpicomm, petCount=npet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fms_init(mpicomm)

    ! ---------------------
    ! Open mosaic file and query some information 
    ! ---------------------

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
    if (dbug_flag > 2) then
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
    ! Find out layout that will be used to read the data
    !----------------------

    ! setup global indices
    do n = 1, domain%ntiles
       global_indices(1,n) = 1
       global_indices(2,n) = domain%nit(n)
       global_indices(3,n) = 1
       global_indices(4,n) = domain%njt(n)
    end do

    ! check total number of PETs
    if (mod(npet, domain%ntiles) /= 0) then
       write(msg, fmt='(A,I5)') trim(subname)//' : nPet should be multiple of 6 to read initial conditions but it is ', npet 
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! calculate layout
    npes_per_tile = npet/domain%ntiles
    call mpp_define_layout(global_indices(:,1), npes_per_tile, domain%layout)

    ! set layout and print out debug information
    do n = 1, domain%ntiles
       layout2d(:,n) = domain%layout(:)
       if (dbug_flag > 2) then
          write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' layout (', n ,') = ', layout2d(1,n), layout2d(2,n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
          write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' global_indices (', n,') = ', &
            global_indices(1,n), global_indices(2,n), global_indices(3,n), global_indices(4,n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end if
    enddo

    !----------------------
    ! Set pe_start, pe_end 
    !----------------------

    allocate(pe_start(domain%ntiles))
    allocate(pe_end(domain%ntiles))
    do n = 1, domain%ntiles
       pe_start(n) = mpp_root_pe()+(n-1)*domain%layout(1)*domain%layout(2)
       pe_end(n) = mpp_root_pe()+n*domain%layout(1)*domain%layout(2)-1
       if (dbug_flag > 2) then
          write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' pe_start, pe_end (', n ,') = ', pe_start(n), pe_end(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end if 
    enddo

    !----------------------
    ! Create FMS domain object
    !----------------------

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
  subroutine create_grid(gcomp, domain, mosaic_file, input_dir, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(domain_type), intent(inout) :: domain
    character(len=cl), intent(in)    :: mosaic_file
    character(len=cl), intent(in)    :: input_dir
    integer, intent(inout)           :: rc

    ! local variables
    type(ESMF_Decomp_Flag)           :: decompflagPTile(2,6)
    integer                          :: n
    integer                          :: decomptile(2,6)
    character(len=*), parameter      :: subname = trim(modName)//': (create_grid) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

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

    ! create mesh
    domain%mesh = ESMF_MeshCreate(domain%grid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine create_grid

  !===============================================================================
  subroutine read_tiled_file(gcomp, filename, varname, domain, field_dst, rc)
    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)          :: gcomp
    character(len=*), intent(in)             :: filename
    character(len=*), intent(in)             :: varname 
    type(domain_type), intent(inout)         :: domain
    type(ESMF_Field), intent(inout)          :: field_dst
    integer, intent(inout), optional         :: rc

    ! local variables
    type(ESMF_Field)            :: field_src, field_tmp
    type(ESMF_ArraySpec)        :: arraySpec
    type(InternalState)         :: is_local
    type(fieldtype), allocatable:: vars(:) 
    integer                     :: funit, my_tile
    integer                     :: i, j, n
    integer                     :: isc, iec, jsc, jec
    integer                     :: ndim, nvar, natt, ntime
    logical                     :: not_found, is_root_pe
    real(ESMF_KIND_R8), pointer :: ptr2d(:,:)
    real(r8), allocatable       :: rdata(:,:)
    character(len=cl)           :: cname
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
    ! Set tile
    !----------------------
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
          allocate(rdata(isc:iec,jsc:jec))
          rdata(:,:) = 0.0_r8

          ! read data
          call mpp_read(funit, vars(n), domain%mosaic_domain, rdata, 1)

          ! set missing values to zero
          where (rdata == 1.0e20)
             rdata(:,:) = 0.0_r8
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
    call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create source field
    field_src = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
       indexflag=ESMF_INDEX_GLOBAL, name=trim(varname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and fill it
    call ESMF_FieldGet(field_src, localDe=0, farrayPtr=ptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ptr2d(:,:) = rdata(:,:)
    nullify(ptr2d)
    if (allocated(rdata)) deallocate(rdata)

    ! create destination field
    field_dst = ESMF_FieldCreate(domain%mesh, ESMF_TYPEKIND_R8, name=trim(varname), &
                meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create routehandle from grid to mesh
    if (.not. ESMF_RouteHandleIsCreated(domain%rh, rc=rc)) then
       call ESMF_FieldRegridStore(field_src, field_dst, routehandle=domain%rh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! redist field from ESMF Grid to Mesh
    call ESMF_FieldRedist(field_src, field_dst, domain%rh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Output result field for debugging purpose
    !----------------------

    if (dbug_flag > 2) then
       call ESMF_FieldWrite(field_dst, trim(varname)//'agrid', variableName=trim(varname), overwrite=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 5) then
       call ESMF_FieldWriteVTK(field_dst, trim(varname)//'agrid', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! clean memory
    call ESMF_FieldDestroy(field_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return 

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
    type(ESMF_FieldBundle), save :: FBout
    type(InternalState) :: is_local
    integer           :: yr, mon, day, sec
    integer           :: n, m, ns, start_ymd
    character(cl)     :: time_units
    real(r8)          :: time_val
    real(r8)          :: time_bnds(2)
    real(r8), pointer :: ptr(:)
    character(len=cl) :: tmpstr
    character(len=cl) :: rst_file
    character(len=cl) :: nexttime_str
    integer, save     :: ns_total
    logical, save     :: first_call = .true.
    character(len=cs), allocatable :: flds(:)
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
    rst_file = trim(case_name)//'.cpl.ccpp.'//trim(nexttime_str)//'.nc'

    ! return if it is not time to write restart
    if (restart_freq < 0) return
    if (mod(sec, restart_freq) /= 0) return

    !----------------------
    ! Create restart file
    !----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
       FBout = ESMF_FieldBundleCreate(rc=rc)

       ! get total element count 
       call ESMF_MeshGet(is_local%wrap%aoflux_mesh, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMAllFullReduce(vm, (/ns/), ns_total, 1, ESMF_REDUCE_SUM, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! add fields
       allocate(flds(3))
       flds = (/ 'zorl  ', &
                 'uustar', &
                 'qss   ' /)
       do n = 1,size(flds)
          ! create new field on aoflux mesh
          field = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
             name=trim(flds(n)), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
 
          ! get pointer out of field
          call ESMF_FieldGet(field, farrayptr=ptr, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! fill pointer
          if (trim(flds(n)) == 'zorl'  ) ptr(:) = physics%sfcprop%zorl(:)
          if (trim(flds(n)) == 'uustar') ptr(:) = physics%sfcprop%uustar(:)
          if (trim(flds(n)) == 'qss'   ) ptr(:) = physics%sfcprop%qss(:)
          nullify(ptr)

          ! add field to FB
          call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
    else
       do n = 1,size(flds)
          ! retrieve field pointer from FB
          call fldbun_getdata1d(FBout, trim(flds(n)), ptr, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! fill pointer
          if (trim(flds(n)) == 'zorl'  ) ptr(:) = physics%sfcprop%zorl(:)
          if (trim(flds(n)) == 'uustar') ptr(:) = physics%sfcprop%uustar(:)
          if (trim(flds(n)) == 'qss'   ) ptr(:) = physics%sfcprop%qss(:)
          nullify(ptr)
       end do
    end if

    ! debug
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//' diagnose at '//trim(nexttime_str), ESMF_LOGMSG_INFO)
       call fldbun_diagnose(FBout, string=trim(subname)//' CCPP FBout ', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! debug
    if (dbug_flag > 5) then
       do n = 1,size(flds)
          ! retrieve field from FB
          call ESMF_FieldBundleGet(FBout, fieldName=trim(flds(n)), field=field, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! write field in VTK format
          call ESMF_FieldWriteVTK(field, 'rst_'//trim(flds(n))//'_'//trim(is_local%wrap%aoflux_grid)//'_'//trim(nexttime_str), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
    end if

    !----------------------
    ! Write data
    !----------------------

    call ESMF_FieldBundleWrite(FBout, trim(rst_file), overwrite=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       write(logunit,'(a)') 'CCPP restart file is closed: '//trim(rst_file)
    end if

  end subroutine write_restart

  end module ufs_io_mod
