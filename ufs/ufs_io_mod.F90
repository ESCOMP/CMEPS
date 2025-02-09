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
  use ESMF,                  only : ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL 
  use ESMF,                  only : ESMF_RouteHandle, ESMF_FieldBundleRedistStore
  use ESMF,                  only : ESMF_FieldBundleRedist, ESMF_RouteHandleDestroy
  use ESMF,                  only : ESMF_TYPEKIND_I4, ESMF_TYPEKIND_R4
  use NUOPC,                 only : NUOPC_CompAttributeGet
  use NUOPC_Mediator,        only : NUOPC_MediatorGet

  use med_kind_mod,          only : r4=>SHR_KIND_R4, r8=>SHR_KIND_R8
  use med_kind_mod,          only : cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use med_utils_mod,         only : chkerr => med_utils_chkerr
  use med_constants_mod,     only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod, only : InternalState, maintask, logunit
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
     integer                :: layout(2)     ! layout for domain decomposition
     integer                :: ntiles        ! number of tiles in case of having CS grid
  end type domain_type

  type field_type
     real(r4), pointer  :: ptr1r4(:)         ! data pointer for 1d r4
     real(r8), pointer  :: ptr1r8(:)         ! data pointer for 1d r8
     integer , pointer  :: ptr1i4(:)         ! data pointer for 1d i4
     real(r4), pointer  :: ptr2r4(:,:)       ! data pointer for 2d r4
     real(r8), pointer  :: ptr2r8(:,:)       ! data pointer for 2d r8
     integer , pointer  :: ptr2i4(:,:)       ! data pointer for 2d i4
     character(len=128) :: short_name = ""   ! variable short name
     character(len=128) :: units = ""        ! variable unit
     character(len=128) :: long_name = ""    ! variable long name
     character(len=128) :: zaxis = ""        ! name of z-axis
     integer            :: nlev              ! number of layers in z-axis
     integer            :: nrec              ! number of record in file (time axis)
  end type field_type

  character(cl) :: case_name = 'unset'  ! case name

  character(*), parameter :: modName = "(ufs_io_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine read_initial(gcomp, ini_file, mosaic_file, input_dir, layout, rc)

    implicit none

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    character(len=cl), intent(in)    :: ini_file
    character(len=cl), intent(in)    :: mosaic_file
    character(len=cl), intent(in)    :: input_dir
    integer                          :: layout(2)
    integer, intent(inout)           :: rc

    ! local variables
    type(domain_type)                :: domain
    type(InternalState)              :: is_local
    type(ESMF_RouteHandle)           :: rh
    type(ESMF_Field)                 :: field_src, field_dst
    real(ESMF_KIND_R8), pointer      :: ptr_src(:), ptr_dst(:)
    integer                          :: n
    character(len=cs)                :: flds_name(2)
    type(field_type)                 :: flds(1)
    character(len=*), parameter      :: subname = trim(modName)//': (read_initial) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Set domain specific parameters
    ! TODO: This assumes global domain with six tile and needs to be 
    ! revisited to support regional apps with one tile
    ! ---------------------

    domain%ntiles = 6
    domain%layout(1) = layout(1)
    domain%layout(2) = layout(2)

    ! ---------------------
    ! Create grid 
    ! ---------------------

    call create_grid(gcomp, domain, mosaic_file, input_dir, rc)

    ! ---------------------
    ! Create field in source mesh
    ! ---------------------

    ! create field
    field_src = ESMF_FieldCreate(domain%mesh, ESMF_TYPEKIND_R8, name='field_src', &
      meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and init it
    call ESMF_FieldGet(field_src, localDe=0, farrayPtr=ptr_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ptr_src(:) = 0.0_r8

    ! ---------------------
    ! Create field in destination mesh
    ! ---------------------

    ! create destination field
    field_dst = ESMF_FieldCreate(is_local%wrap%aoflux_mesh, ESMF_TYPEKIND_R8, &
      name='field_dst', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and init it
    call ESMF_FieldGet(field_dst, localDe=0, farrayPtr=ptr_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ptr_src(:) = 0.0_r8

    ! ---------------------
    ! Create routehandle
    ! ---------------------

    call ESMF_FieldRegridStore(field_src, field_dst, routehandle=rh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read data 
    !----------------------

    ! list of fields that need to be read
    flds_name(1) = 'zorl'
    flds_name(2) = 'uustar'

    ! loop over fields and read them
    do n = 1,size(flds)
       ! read data
       flds(1)%short_name = trim(flds_name(n))
       flds(1)%ptr1r8 => ptr_src
       call read_tiled_file(domain, ini_file, flds, rc=rc) 

       ! map field
       if (is_local%wrap%aoflux_grid == 'agrid') then
          ! do nothing, just redist in case of having different decomp. in here and aoflux mesh
          call ESMF_FieldRedist(field_src, field_dst, rh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          ! remap from atm to ocn or exchange grid
          call ESMF_FieldRegrid(field_src, field_dst, rh, termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! debug
       if (dbug_flag > 5) then
          call ESMF_FieldWriteVTK(field_dst, 'ini_'//trim(flds_name(n))//'_'//trim(is_local%wrap%aoflux_grid), rc=rc)  
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! fill variables
       if (maintask) write(logunit,'(a)') 'Reading: '//trim(flds_name(n))
       if (trim(flds_name(n)) == 'zorl'  ) physics%sfcprop%zorl(:)  = ptr_dst(:)
       if (trim(flds_name(n)) == 'uustar') physics%sfcprop%uustar(:)= ptr_dst(:)
    end do

    !----------------------
    ! Free memory
    !----------------------

    call ESMF_FieldDestroy(field_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldDestroy(field_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    if (maintask) then
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

          if (maintask) write(logunit,'(a)') 'Reading: '//trim(flds(n))
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
  subroutine read_tiled_file(domain, filename, flds, rh, rc)

    ! input/output variables
    type(domain_type), intent(inout) :: domain
    character(len=*),  intent(in)    :: filename
    type(field_type),  intent(in)    :: flds(:)
    type(ESMF_RouteHandle), optional, intent(in) :: rh
    integer, optional, intent(inout) :: rc

    ! local variables
    integer                     :: i, j, k, rank, fieldCount
    integer, pointer            :: ptr_i4(:)
    real(r4), pointer           :: ptr_r4(:)
    real(r8), pointer           :: ptr_r8(:)
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpec
    type(ESMF_Field)            :: fgrid, fmesh, ftmp
    character(len=cl)           :: fname
    character(len=cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname = trim(modName)//': (read_tiled_file) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

    !----------------------
    ! Create field bundles
    !----------------------

    ! create empty field bundle on grid
    FBgrid = ESMF_FieldBundleCreate(name="fields_on_grid", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create empty field bundle on mesh
    FBmesh = ESMF_FieldBundleCreate(name="fields_on_mesh", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Loop over fields and add them to the field bundles
    !----------------------

    do i = 1, size(flds)
       ! 2d/r8 field (x,y)
       if (associated(flds(i)%ptr1r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/r4 field (x,y)
       else if (associated(flds(i)%ptr1r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr1r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/i4 field (x,y)
       else if (associated(flds(i)%ptr1i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_I4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr1i4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r8 field (x,y,rec)
       else if (associated(flds(i)%ptr2r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr2r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r4 field (x,y,rec)
       else if (associated(flds(i)%ptr2r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr2r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/i4 field (x,y,rec)
       else if (associated(flds(i)%ptr2i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_I4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(domain%mesh, flds(i)%ptr2i4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! debug print
       call ESMF_LogWrite(trim(subname)//' adding '//trim(flds(i)%short_name)//' to FB', ESMF_LOGMSG_INFO)

       ! add it to the field bundle on grid
       call ESMF_FieldBundleAdd(FBgrid, [fgrid], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! add it to the field bundle on mesh
       call ESMF_FieldBundleAdd(FBmesh, [fmesh], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !----------------------
    ! Read data
    !----------------------

    call ESMF_FieldBundleRead(FBgrid, fileName=trim(filename), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Create routehandle if it is not provided to transfer data from grid to mesh
    !----------------------

    if (present(rh)) then
       rh_local = rh
    else
       call ESMF_FieldBundleRedistStore(FBgrid, FBmesh, routehandle=rh_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Move data from ESMF grid to mesh
    !----------------------

    call ESMF_FieldBundleRedist(FBgrid, FBmesh, rh_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Debug output
    !----------------------

    if (dbug_flag > 5) then
       do i = 1, size(flds)
          ! get field from FB
          call ESMF_FieldBundleGet(FBmesh, fieldName=trim(flds(i)%short_name), field=fmesh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! check its rank
          call ESMF_FieldGet(fmesh, rank=rank, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! TODO: ESMF_FieldWriteVTK() call does not support ungridded dimension
          ! The workaround is implemented in here but it would be nice to extend
          ! ESMF_FieldWriteVTK() call to handle it.  
          if (rank > 1) then
             ! create temporary field
             if (associated(flds(i)%ptr2r4)) then
                ftmp = ESMF_FieldCreate(domain%mesh, typekind=ESMF_TYPEKIND_R4, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_r4, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else if (associated(flds(i)%ptr2r8)) then
                ftmp = ESMF_FieldCreate(domain%mesh, typekind=ESMF_TYPEKIND_R8, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_r8, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else if (associated(flds(i)%ptr2i4)) then
                ftmp = ESMF_FieldCreate(domain%mesh, typekind=ESMF_TYPEKIND_I4, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_i4, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! write all record to seperate VTK file
             do j = 1, flds(i)%nrec
                if (associated(flds(i)%ptr2i4)) ptr_i4(:) = flds(i)%ptr2i4(:,j)
                if (associated(flds(i)%ptr2r4)) ptr_r4(:) = flds(i)%ptr2r4(:,j)
                if (associated(flds(i)%ptr2r8)) ptr_r8(:) = flds(i)%ptr2r8(:,j)
                write(fname, fmt='(A,I2.2)') trim(flds(i)%short_name)//'_rec', j
                call ESMF_FieldWriteVTK(ftmp, trim(fname), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end do

             ! delete temporary field
             call ESMF_FieldDestroy(ftmp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             ! write field to VTK file
             call ESMF_FieldWriteVTK(fmesh, trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    !----------------------
    ! Empty FBs and destroy them 
    !----------------------

    ! FB grid
    call ESMF_FieldBundleGet(FBgrid, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBgrid, fieldNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do i = 1, fieldCount
       ! pull field from FB
       call ESMF_FieldBundleGet(FBgrid, fieldName=trim(fieldNameList(i)), field=ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy field
       call ESMF_FieldDestroy(ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! remove field from FB
       call ESMF_FieldBundleRemove(FBgrid, fieldNameList=[trim(fieldNameList(i))], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    deallocate(fieldNameList)

    ! destroy grid FB
    call ESMF_FieldBundleDestroy(FBgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! FB mesh 
    call ESMF_FieldBundleGet(FBmesh, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBmesh, fieldNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do i = 1, fieldCount
       ! pull field from FB
       call ESMF_FieldBundleGet(FBmesh, fieldName=trim(fieldNameList(i)), field=ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy field
       call ESMF_FieldDestroy(ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! remove field from FB
       call ESMF_FieldBundleRemove(FBmesh, fieldNameList=[trim(fieldNameList(i))], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    deallocate(fieldNameList)

    ! destroy grid FB
    call ESMF_FieldBundleDestroy(FBmesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Destroy route handle if it is created locally 
    !----------------------

    if (.not. present(rh)) then
       call ESMF_RouteHandleDestroy(rh_local, rc=rc)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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

    if (maintask) then
       write(logunit,'(a)') 'CCPP restart file is closed: '//trim(rst_file)
    end if

  end subroutine write_restart

  end module ufs_io_mod
