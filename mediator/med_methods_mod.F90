module med_methods_mod

  !-----------------------------------------------------------------------------
  ! Generic operation methods used by the Mediator Component.
  !-----------------------------------------------------------------------------

  use med_kind_mod       , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF               , only : operator(<), operator(/=), operator(+), operator(-), operator(*) , operator(>=)
  use ESMF               , only : operator(<=), operator(>), operator(==)
  use ESMF               , only : ESMF_FieldStatus_Flag
  use ESMF               , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF               , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError
  use ESMF               , only : ESMF_MAXSTR, ESMF_LOGMSG_WARNING
  use med_constants_mod  , only : dbug_flag => med_constants_dbug_flag
  use med_constants_mod  , only : czero => med_constants_czero
  use med_constants_mod  , only : spval_init => med_constants_spval_init
  use med_utils_mod      , only : ChkErr => med_utils_ChkErr
  use shr_sys_mod        , only : shr_sys_abort
  implicit none
  private

  interface med_methods_FieldPtr_compare ; module procedure &
    med_methods_FieldPtr_compare1, &
    med_methods_FieldPtr_compare2
  end interface

  interface med_methods_check_for_nans
     module procedure med_methods_check_for_nans_1d
     module procedure med_methods_check_for_nans_2d
  end interface med_methods_check_for_nans

  ! used/reused in module
  logical, public               :: mediator_checkfornans  ! set in med.F90 AdvertiseFields
  logical                       :: isPresent
  character(len=1024)           :: msgString
  type(ESMF_FieldStatus_Flag)   :: status
  character(*)      , parameter :: u_FILE_u = &
       __FILE__

  public med_methods_FB_copy
  public med_methods_FB_accum
  public med_methods_FB_average
  public med_methods_FB_init
  public med_methods_FB_init_pointer
  public med_methods_FB_reset
  public med_methods_FB_diagnose
  public med_methods_FB_write
  public med_methods_FB_FldChk
  public med_methods_FB_GetFldPtr
  public med_methods_FB_getNameN
  public med_methods_FB_getFieldN
  public med_methods_FB_getNumflds
  public med_methods_FB_Field_diagnose
  public med_methods_FB_GeomPrint
  public med_methods_FB_getdata2d
  public med_methods_FB_getdata1d
  public med_methods_FB_getmesh
  public med_methods_FB_check_for_nans

  public med_methods_State_reset
  public med_methods_State_diagnose
  public med_methods_State_GeomPrint
  public med_methods_State_SetScalar
  public med_methods_State_GetScalar
  public med_methods_State_GetNumFields

  public med_methods_Field_getdata1d
  public med_methods_Field_getdata2d
  public med_methods_Field_diagnose
  public med_methods_Field_GeomPrint
  public med_methods_FieldPtr_compare

  public med_methods_Clock_TimePrint

  private med_methods_Mesh_Print
  private med_methods_Grid_Print
  private med_methods_Field_GetFldPtr
#ifdef DIAGNOSE
  private med_methods_Array_diagnose
#endif
  private med_methods_check_for_nans

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_methods_FB_init_pointer(StateIn, FBout, flds_scalar_name, name, rc)

    ! ----------------------------------------------
    ! Create FBout from StateIn mesh and pointer
    ! ----------------------------------------------

    use ESMF , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleAdd, ESMF_FieldBundleCreate
    use ESMF , only : ESMF_State, ESMF_StateGet, ESMF_Mesh, ESMF_MeshLoc
    use ESMF , only : ESMF_AttributeGet, ESMF_INDEX_DELOCAL

    ! input/output variables
    type(ESMF_State)      , intent(in)           :: StateIn          ! input state
    type(ESMF_FieldBundle), intent(inout)        :: FBout            ! output field bundle
    character(len=*)      , intent(in)           :: flds_scalar_name ! name of scalar fields
    character(len=*)      , intent(in)           :: name
    integer               , intent(out)          :: rc

    ! local variables
    logical            :: isPresent
    integer            :: n,n1
    type(ESMF_Field)   :: lfield
    type(ESMF_Field)   :: newfield
    type(ESMF_MeshLoc) :: meshloc
    type(ESMF_Mesh)    :: lmesh
    integer            :: lrank
    integer            :: fieldCount
    integer            :: ungriddedCount
    integer            :: ungriddedLBound(1)
    integer            :: ungriddedUBound(1)
    real(R8), pointer  :: dataptr1d(:)
    real(R8), pointer  :: dataptr2d(:,:)
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(len=*), parameter :: subname='(med_methods_FB_init_pointer)'
    ! ----------------------------------------------

    ! Create empty FBout
    FBout = ESMF_FieldBundleCreate(name=trim(name), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get fields from StateIn
    call ESMF_StateGet(StateIn, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldNameList(fieldCount))
    call ESMF_StateGet(StateIn, itemNameList=lfieldNameList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Remove scalar field and blank fields from field bundle
    do n = 1, fieldCount
      if (trim(lfieldnamelist(n)) == trim(flds_scalar_name) .or. trim(lfieldnamelist(n)) == '') then
        do n1 = n, fieldCount-1
          lfieldnamelist(n1) = lfieldnamelist(n1+1)
        enddo
        fieldCount = fieldCount - 1
      endif
    enddo  ! n

    ! Only create the fieldbundle if the number of non-scalar fields is > 0
    if (fieldCount > 0) then

       ! Get mesh from first non-scalar field in StateIn (assumes all the fields have the same mesh)
       call ESMF_StateGet(StateIn, itemName=lfieldNameList(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh, meshloc=meshloc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Loop over fields in StateIn skipping the field with just scalar data
       do n = 1, fieldCount
          ! get field from StateIn
          call ESMF_StateGet(StateIn, itemName=lfieldNameList(n), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! determine rank of field
          call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (lrank == 2) then

             ! determine ungridded lower and upper bounds for lfield
             call ESMF_AttributeGet(lfield, name="UngriddedLBound", convention="NUOPC", &
                  purpose="Instance", itemCount=ungriddedCount,  isPresent=isPresent, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (ungriddedCount /= 1) then
                call shr_sys_abort(trim(subname)//": ERROR ungriddedCount for "// &
                     trim(lfieldnamelist(n))//" must be 1 if rank is 2 ")
             end if

             ! set ungridded dimensions for field
             call ESMF_AttributeGet(lfield, name="UngriddedLBound", convention="NUOPC", &
                  purpose="Instance", valueList=ungriddedLBound, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_AttributeGet(lfield, name="UngriddedUBound", convention="NUOPC", &
                  purpose="Instance", valueList=ungriddedUBound, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get 2d pointer for field
             call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! create new field with an ungridded dimension
             newfield = ESMF_FieldCreate(lmesh, dataptr2d, ESMF_INDEX_DELOCAL, &
                  meshloc=meshloc, name=lfieldNameList(n), &
                  ungriddedLbound=ungriddedLbound, ungriddedUbound=ungriddedUbound, gridToFieldMap=(/2/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          else if (lrank == 1) then

             ! get 1d pointer for field
             call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! create new field without an ungridded dimension
             newfield = ESMF_FieldCreate(lmesh, dataptr1d, ESMF_INDEX_DELOCAL, &
                  meshloc=meshloc, name=lfieldNameList(n), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          else
             call shr_sys_abort(trim(subname)//": ERROR only rank1 and rank2 are supported for rank of fields ")
          end if

          ! Add new field to FBout
          call ESMF_FieldBundleAdd(FBout, (/newfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

       end do ! end of loop over input state fields
    end if  ! end of fieldcount > 0

    deallocate(lfieldNameList)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": FBout from input State and field pointers", ESMF_LOGMSG_INFO)
    end if

  end subroutine med_methods_FB_init_pointer

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_init(FBout, flds_scalar_name, fieldNameList, FBgeom, STgeom, FBflds, STflds, name, rc)

    ! ----------------------------------------------
    ! Create FBout from fieldNameList, FBflds, STflds, FBgeom or STgeom in that order or priority
    ! Pass in FBgeom OR STgeom, get mesh from that object
    ! ----------------------------------------------

    use ESMF , only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleGet
    use ESMF , only : ESMF_State, ESMF_Mesh, ESMF_StaggerLoc, ESMF_MeshLoc
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_FieldBundleAdd, ESMF_FieldCreate
    use ESMF , only : ESMF_TYPEKIND_R8, ESMF_FIELDSTATUS_EMPTY, ESMF_AttributeGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)        :: FBout            ! output field bundle
    character(len=*)      , intent(in)           :: flds_scalar_name ! name of scalar fields
    character(len=*)      , intent(in), optional :: fieldNameList(:) ! names of fields to use in output field bundle
    type(ESMF_FieldBundle), intent(in), optional :: FBgeom           ! input field bundle geometry to use
    type(ESMF_State)      , intent(in), optional :: STgeom           ! input state geometry to use
    type(ESMF_FieldBundle), intent(in), optional :: FBflds           ! input field bundle fields
    type(ESMF_State)      , intent(in), optional :: STflds           ! input state fields
    character(len=*)      , intent(in), optional :: name             ! name to use for output field bundle
    integer               , intent(out)          :: rc

    ! local variables
    integer                :: n,n1
    integer                :: fieldCount,fieldCountgeom
    character(ESMF_MAXSTR) :: lname
    type(ESMF_Field)       :: field,lfield
    type(ESMF_Mesh)        :: lmesh
    type(ESMF_MeshLoc)     :: meshloc
    integer                :: ungriddedCount
    integer                :: ungriddedCount_in
    integer, allocatable   :: ungriddedLBound(:)
    integer, allocatable   :: ungriddedUBound(:)
    logical                :: isPresent
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(len=*), parameter :: subname='(med_methods_FB_init)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lname = 'undefined'
    if (present(name)) then
      lname = trim(name)
    endif
    lname = 'FB '//trim(lname)

    !---------------------------------
    ! check argument consistency and
    ! verify that geom argument has a field
    !---------------------------------

    if (present(fieldNameList) .and. present(FBflds) .and. present(STflds)) then
      call shr_sys_abort(trim(subname)//": ERROR only fieldNameList, FBflds, or STflds can be an argument")
    endif

    if (present(FBgeom) .and. present(STgeom)) then
       call shr_sys_abort(trim(subname)//": ERROR FBgeom and STgeom cannot both be arguments")
    endif

    if (.not.present(FBgeom) .and. .not.present(STgeom)) then
       call shr_sys_abort(trim(subname)//": ERROR FBgeom or STgeom must be an argument")
    endif

    if (present(FBgeom)) then
      call ESMF_FieldBundleGet(FBgeom, fieldCount=fieldCountGeom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (present(STgeom)) then
      call ESMF_StateGet(STgeom, itemCount=fieldCountGeom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call shr_sys_abort(trim(subname)//": ERROR FBgeom or STgeom must be passed")
    endif

    !---------------------------------
    ! determine the names of fields that will be in FBout
    !---------------------------------

    if (present(fieldNameList)) then
      fieldcount = size(fieldNameList)
      allocate(lfieldNameList(fieldcount))
      lfieldNameList = fieldNameList
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from argument", ESMF_LOGMSG_INFO)
      end if
    elseif (present(FBflds)) then
      call ESMF_FieldBundleGet(FBflds, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_FieldBundleGet(FBflds, fieldNameList=lfieldNameList, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from FBflds", ESMF_LOGMSG_INFO)
      end if
    elseif (present(STflds)) then
      call ESMF_StateGet(STflds, itemCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_StateGet(STflds, itemNameList=lfieldNameList, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from STflds", ESMF_LOGMSG_INFO)
      end if
    elseif (present(FBgeom)) then
      call ESMF_FieldBundleGet(FBgeom, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_FieldBundleGet(FBgeom, fieldNameList=lfieldNameList, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from FBgeom", ESMF_LOGMSG_INFO)
      end if
    elseif (present(STgeom)) then
      call ESMF_StateGet(STgeom, itemCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_StateGet(STgeom, itemNameList=lfieldNameList, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from STgeom", ESMF_LOGMSG_INFO)
      end if
    else
       call shr_sys_abort(trim(subname)//": ERROR fieldNameList, FBflds, STflds, FBgeom, or STgeom must be passed")
    endif

    !---------------------------------
    ! remove scalar field and blank fields from field bundle
    !---------------------------------

    do n = 1, fieldCount
      if (trim(lfieldnamelist(n)) == trim(flds_scalar_name) .or. &
          trim(lfieldnamelist(n)) == '') then
        do n1 = n, fieldCount-1
           lfieldnamelist(n1) = lfieldnamelist(n1+1)
        enddo
        fieldCount = fieldCount - 1
      endif
    enddo  ! n

    !---------------------------------
    ! create the mesh(lmesh) that will be used for FBout fields
    !---------------------------------

    if (fieldcount > 0 .and. fieldcountgeom > 0) then

      ! Look at only the first field in either the FBgeom and STgeom to get the mesh
      if (present(FBgeom)) then
        call med_methods_FB_getFieldN(FBgeom, 1, lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" mesh from FBgeom", ESMF_LOGMSG_INFO)
        end if
      elseif (present(STgeom)) then
        call med_methods_State_getNameN(STgeom, 1, lname, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_StateGet(STgeom, itemName=lname, field=lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" mesh from STgeom", ESMF_LOGMSG_INFO)
        end if
      else
        call shr_sys_abort(trim(subname)//": ERROR FBgeom or STgeom must be passed")
      endif

      ! Make sure the field is not empty - if it is abort with an error
      call ESMF_FieldGet(lfield, status=status, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (status == ESMF_FIELDSTATUS_EMPTY) then
         call shr_sys_abort(trim(subname)//":"//trim(lname)//": ERROR field does not have a geom yet ", &
              line=__LINE__, file=u_FILE_u)
      endif

      ! Assume field is on mesh
      call ESMF_FieldGet(lfield, mesh=lmesh, meshloc=meshloc, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" use mesh", ESMF_LOGMSG_INFO)
      end if

    endif  ! fieldcount > 0

    !---------------------------------
    ! create FBout
    !---------------------------------

    FBout = ESMF_FieldBundleCreate(name=trim(lname), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fieldcountgeom > 0) then

       ! Now loop over all the fields in the field name list
       do n = 1, fieldCount

          ! Note that input fields come from ONE of FBFlds, STflds, or fieldNamelist input argument
          if (present(FBFlds) .or. present(STflds)) then

             ! ungridded dimensions might be present in the input states or field bundles
             if (present(FBflds)) then
                call ESMF_FieldBundleGet(FBflds, fieldName=lfieldnamelist(n), field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             elseif (present(STflds)) then
                call med_methods_State_getNameN(STflds, n, lname, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_StateGet(STflds, itemName=lname, field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             end if

             ! Determine ungridded lower and upper bounds for lfield
             call ESMF_AttributeGet(lfield, name="UngriddedUBound", convention="NUOPC", &
                  purpose="Instance", itemCount=ungriddedCount_in,  isPresent=isPresent, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (isPresent) then
                ungriddedCount = ungriddedCount_in
             else
                ungriddedCount=0  ! initialize in case it was not set
             end if

             ! Create the field on a lmesh
             if (ungriddedCount > 0) then
                ! ungridded dimensions in field
                allocate(ungriddedLBound(ungriddedCount), ungriddedUBound(ungriddedCount))
                call ESMF_AttributeGet(lfield, name="UngriddedLBound", convention="NUOPC", &
                     purpose="Instance", valueList=ungriddedLBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_AttributeGet(lfield, name="UngriddedUBound", convention="NUOPC", &
                     purpose="Instance", valueList=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=lfieldNameList(n), &
                     ungriddedLbound=ungriddedLbound, ungriddedUbound=ungriddedUbound, gridToFieldMap=(/2/))
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                deallocate( ungriddedLbound, ungriddedUbound)
             else
                ! No ungridded dimensions in field
                field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=lfieldNameList(n), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             end if

          else if (present(fieldNameList)) then

             ! Assume no ungridded dimensions if just the field name list is give
             field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=lfieldNameList(n), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          end if

          ! Add the created field bundle FBout
          if (dbug_flag > 1) then
             call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" adding field "//trim(lfieldNameList(n)), &
                  ESMF_LOGMSG_INFO)
          end if
          call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

       enddo  ! fieldCount
    endif  ! fieldcountgeom

    deallocate(lfieldNameList)

    call med_methods_FB_reset(FBout, value=spval_init, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_init

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_getNameN(FB, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get name of field number fieldnum in input field bundle FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    character(len=*)      , intent(out)   :: fieldname
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(med_methods_FB_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldnum > fieldCount) then
      call shr_sys_abort(trim(subname)//": ERROR fieldnum > fieldCount ")
    endif
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    fieldname = lfieldnamelist(fieldnum)
    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_getNameN

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_getFieldN(FB, fieldnum, field, rc)

    ! ----------------------------------------------
    ! Get field with number fieldnum in input field bundle FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    character(len=*),parameter :: subname='(med_methods_FB_getFieldN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call med_methods_FB_getNameN(FB, fieldnum, name, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_getFieldN

  !-----------------------------------------------------------------------------

  subroutine med_methods_State_getNameN(State, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get field number fieldnum name out of State
    ! ----------------------------------------------

    use ESMF, only : ESMF_State, ESMF_StateGet

    type(ESMF_State), intent(in)    :: State
    integer         , intent(in)    :: fieldnum
    character(len=*), intent(out)   :: fieldname
    integer         , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(med_methods_State_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '
    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldnum > fieldCount) then
      call shr_sys_abort(trim(subname)//": ERROR fieldnum > fieldCount ")
    endif
    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    fieldname = lfieldnamelist(fieldnum)
    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_State_getNameN

  !-----------------------------------------------------------------------------

  subroutine med_methods_State_getNumFields(State, fieldnum, rc)

    ! ----------------------------------------------
    ! Get field number fieldnum name out of State
    ! ----------------------------------------------

    use NUOPC , only : NUOPC_GetStateMemberLists
    use ESMF  , only : ESMF_State, ESMF_Field, ESMF_StateGet, ESMF_STATEITEM_FIELD
    use ESMF  , only : ESMF_StateItem_Flag

    type(ESMF_State), intent(in)    :: State
    integer         , intent(inout) :: fieldnum
    integer         , intent(out)   :: rc

    ! local variables
    type(ESMF_Field), pointer          :: fieldList(:)
    character(len=*),parameter         :: subname='(med_methods_State_getNumFields)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    nullify(fieldList)
    call NUOPC_GetStateMemberLists(state, fieldList=fieldList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    fieldnum = 0
    if (associated(fieldList)) then
       fieldnum = size(fieldList)
       deallocate(fieldList)
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_State_getNumFields

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_reset(FB, value, rc)
    ! ----------------------------------------------
    ! Set all fields to value in FB
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field
    use ESMF, only : ESMF_FieldBundleIsCreated

    ! intput/output variables
    type(ESMF_FieldBundle) , intent(inout)        :: FB
    real(R8)               , intent(in), optional :: value
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8)                        :: lvalue
    type(ESMF_Field)                :: lfield
    integer                         :: lrank
    real(R8), pointer               :: fldptr1(:)
    real(R8), pointer               :: fldptr2(:,:)
    character(len=*),parameter      :: subname='(med_methods_FB_reset)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lvalue = czero
    if (present(value)) then
      lvalue = value
    endif

    if (ESMF_FieldBundleIsCreated(fieldbundle=FB)) then
       call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lfieldnamelist(fieldCount))
       call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       fieldCount=0
    endif

    do n = 1, fieldCount
       call ESMF_FieldBundleGet(FB, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_methods_Field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          fldptr1 = lvalue
       elseif (lrank == 2) then
          fldptr2 = lvalue
       else
          call shr_sys_abort(trim(subname)//": ERROR in rank "//trim(lfieldnamelist(n)), &
                line=__LINE__, file=u_FILE_u)
       endif
    enddo

    if (ESMF_FieldBundleIsCreated(fieldbundle=FB)) then
      deallocate(lfieldnamelist)
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_reset

  !-----------------------------------------------------------------------------

  subroutine med_methods_State_reset(State, value, rc)

    ! ----------------------------------------------
    ! Set all fields to value in State
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------

    use ESMF, only : ESMF_State, ESMF_StateGet, ESMF_Field

    ! intput/output variables
    type(ESMF_State) , intent(inout)        :: State
    real(R8)         , intent(in), optional :: value
    integer          , intent(out)          :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8)                        :: lvalue
    type(ESMF_Field)                :: lfield
    integer                         :: lrank
    real(R8), pointer               :: fldptr1(:)
    real(R8), pointer               :: fldptr2(:,:)
    character(len=*),parameter      :: subname='(med_methods_State_reset)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lvalue = czero
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldCount
       call ESMF_StateGet(State, itemName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_methods_Field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          fldptr1 = lvalue
       elseif (lrank == 2) then
          fldptr2 = lvalue
       else
          call shr_sys_abort(trim(subname)//": ERROR in rank "//trim(lfieldnamelist(n)), &
               line=__LINE__, file=u_FILE_u)
       endif
    enddo
    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_State_reset

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_average(FB, count, rc)

    ! ----------------------------------------------
    ! Set all fields to zero in FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout) :: FB
    integer               , intent(in)    :: count
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8), pointer               :: dataPtr1(:)
    real(R8), pointer               :: dataPtr2(:,:)
    type(ESMF_Field)                :: lfield
    character(len=*),parameter      :: subname='(med_methods_FB_average)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    if (count == 0) then

       if (dbug_flag > 10) then
          call ESMF_LogWrite(trim(subname)//": WARNING count is 0", ESMF_LOGMSG_INFO)
       end if
       !call ESMF_LogWrite(trim(subname)//": WARNING count is 0 set avg to spval", ESMF_LOGMSG_INFO)
       !call med_methods_FB_reset(FB, value=spval, rc=rc)
       !if (chkerr(rc,__LINE__,u_FILE_u)) return

    else

      call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldnamelist(fieldCount))
      call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      do n = 1, fieldCount
         call ESMF_FieldBundleGet(FB, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call med_methods_Field_GetFldPtr(lfield, fldptr1=dataptr1, fldptr2=dataptr2, rank=lrank, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

        if (lrank == 0) then
          ! no local data
        elseif (lrank == 1) then
          do i=lbound(dataptr1,1),ubound(dataptr1,1)
            dataptr1(i) = dataptr1(i) / real(count, R8)
          enddo
        elseif (lrank == 2) then
          do j=lbound(dataptr2,2),ubound(dataptr2,2)
          do i=lbound(dataptr2,1),ubound(dataptr2,1)
            dataptr2(i,j) = dataptr2(i,j) / real(count, R8)
          enddo
          enddo
        else
          call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
        endif
      enddo
      deallocate(lfieldnamelist)

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_average

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_diagnose(FB, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    type(ESMF_Field)                :: lfield
    character(len=*), parameter     :: subname='(med_methods_FB_diagnose)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string) // ' '
    endif

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call ESMF_FieldBundleGet(FB, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_methods_Field_GetFldPtr(lfield, fldptr1=dataptr1d, fldptr2=dataptr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine med_methods_FB_diagnose

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_write(FB, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF, only : ESMF_Field, ESMF_FieldGet
    use ESMF, only : ESMF_FieldWriteVTK

    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    type(ESMF_Field)                :: lfield
    character(len=*), parameter     :: subname='(med_methods_FB_write)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call ESMF_FieldBundleGet(FB, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 1) then
          call ESMF_FieldWriteVTK(lfield, trim(lfieldnamelist(n))//'_'//trim(lstring), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
       endif
    end do

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine med_methods_FB_write

  !-----------------------------------------------------------------------------
#ifdef DIAGNOSE
  subroutine med_methods_Array_diagnose(array, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of Array
    ! ----------------------------------------------

    use ESMF, only : ESMF_Array, ESMF_ArrayGet

    ! input/output variables
    type(ESMF_Array), intent(inout)        :: array
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    character(len=CS) :: lstring
    real(R8), pointer :: dataPtr3d(:,:,:)
    character(len=*),parameter  :: subname='(med_methods_Array_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! this is not working yet, not sure about dataPtr dim/type
    return

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_ArrayGet(Array, farrayPtr=dataPtr3d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    write(msgString,'(A,3g14.7)') trim(subname)//' '//trim(lstring), &
        minval(dataPtr3d), maxval(dataPtr3d), sum(dataPtr3d)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    end if

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Array_diagnose
#endif
  !-----------------------------------------------------------------------------

  subroutine med_methods_State_diagnose(State, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    use ESMF, only : ESMF_State, ESMF_StateGet, ESMF_Field

    ! input/output variables
    type(ESMF_State), intent(in)           :: State
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=CS)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    type(ESMF_Field)                :: lfield
    character(len=*),parameter      :: subname='(med_methods_State_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    endif

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldCount
       call ESMF_StateGet(State, itemName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_methods_Field_GetFldPtr(lfield, fldptr1=dataptr1d, fldptr2=dataptr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif
       else
          call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
       endif

       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
   endif

  end subroutine med_methods_State_diagnose

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_Field_diagnose(FB, fieldname, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)  :: FB
    character(len=*), intent(in)           :: fieldname
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    character(len=CS) :: lstring
    real(R8), pointer :: dataPtr1d(:)
    real(R8), pointer :: dataPtr2d(:,:)
    type(ESMF_Field)  :: lfield
    integer           :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    character(len=*),parameter :: subname='(med_methods_FB_Field_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_FieldBundleGet(FB, fieldName=fieldname, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (ungriddedUBound(1) > 0) then
       call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (size(dataptr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    else
       call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    end if
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_Field_diagnose

  !-----------------------------------------------------------------------------

  subroutine med_methods_Field_diagnose(field, fieldname, string, rc)

    ! ----------------------------------------------
    ! Diagnose Field
    ! ----------------------------------------------

    use ESMF, only : ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_Field) , intent(inout)        :: field
    character(len=*) , intent(in)           :: fieldname
    character(len=*) , intent(in), optional :: string
    integer          , intent(out)          :: rc

    ! local variables
    integer                    :: lrank
    character(len=CS)          :: lstring
    real(R8), pointer          :: dataPtr1d(:)
    real(R8), pointer          :: dataPtr2d(:,:)
    character(len=*),parameter :: subname='(med_methods_Field_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_FieldGet(field, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
       ! no local data
    elseif (lrank == 1) then
       call ESMF_FieldGet(field, farrayPtr=dataPtr1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    elseif (lrank == 2) then
       call ESMF_FieldGet(field, farrayPtr=dataPtr2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (size(dataPtr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    else
       call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
    endif
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Field_diagnose

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_copy(FBout, FBin, rc)

    ! ----------------------------------------------
    ! Copy common field names from FBin to FBout
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    integer               , intent(out)   :: rc
    character(len=*), parameter :: subname='(med_methods_FB_copy)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    call med_methods_FB_accum(FBout, FBin, copy=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_copy

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_accum(FBout, FBin, copy, rc)

    ! ----------------------------------------------
    ! Accumulate common field names from FBin to FBout
    ! If copy is passed in and true, the this is a copy
    ! ----------------------------------------------

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    logical, optional     , intent(in)    :: copy
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lranki, lranko
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    logical                         :: exists
    logical                         :: lcopy
    real(R8), pointer               :: dataPtri1(:)
    real(R8), pointer               :: dataPtro1(:)
    real(R8), pointer               :: dataPtri2(:,:)
    real(R8), pointer               :: dataPtro2(:,:)
    type(ESMF_Field)                :: lfield
    character(len=*), parameter     :: subname='(med_methods_FB_accum)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lcopy = .false.  ! accumulate by default
    if (present(copy)) then
      lcopy = copy
    endif

    call ESMF_FieldBundleGet(FBout, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FBout, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call ESMF_FieldBundleGet(FBin, fieldName=lfieldnamelist(n), isPresent=exists, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (exists) then
        call ESMF_FieldBundleGet(FBin, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call med_methods_Field_GetFldPtr(lfield, fldptr1=dataptri1, fldptr2=dataptri2, rank=lranki, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        call ESMF_FieldBundleGet(FBout, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call med_methods_Field_GetFldPtr(lfield, fldptr1=dataptro1, fldptr2=dataptro2, rank=lranko, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        if (lranki == 0 .and. lranko == 0) then
           ! do nothing
          call ESMF_LogWrite(trim(subname)//": Both ranki and ranko are 0", ESMF_LOGMSG_INFO)
        elseif (lranki == 1 .and. lranko == 1) then

          if (.not.med_methods_FieldPtr_Compare(dataPtro1, dataPtri1, subname, rc)) then
            call shr_sys_abort(trim(subname)//": ERROR in dataPtr1 size ")
          endif

          if (lcopy) then
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtri1(i)
            enddo
          else
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtro1(i) + dataPtri1(i)
            enddo
          endif

        elseif (lranki == 2 .and. lranko == 2) then

          if (.not.med_methods_FieldPtr_Compare(dataPtro2, dataPtri2, subname, rc)) then
            call shr_sys_abort(trim(subname)//": ERROR in dataPtr2 size ")
          endif

          if (lcopy) then
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtri2(i,j)
            enddo
            enddo
          else
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtro2(i,j) + dataPtri2(i,j)
            enddo
            enddo
          endif

        else
          write(msgString,'(a,2i8)') trim(subname)//": ranki, ranko = ",lranki,lranko
          call shr_sys_abort(trim(subname)//": ERROR ranki ranko not supported "//trim(msgstring)//"\n"//trim(lfieldnamelist(n)))
        endif

      endif
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_accum

  !-----------------------------------------------------------------------------

  logical function med_methods_FB_FldChk(FB, fldname, rc)

    ! ----------------------------------------------
    ! Determine if field with fldname is in input field bundle
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF, only : ESMF_FieldBundleIsCreated

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    integer               , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_methods_FB_FldChk)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! If field bundle is not created then set return to .false.
    if (.not. ESMF_FieldBundleIsCreated(FB)) then
       med_methods_FB_FldChk = .false.
       return
    end if

    ! If field bundle is created determine if fldname is present in field bundle
    med_methods_FB_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) then
       call shr_sys_abort(string=trim(subname)//" Error checking field: "//trim(fldname), line=__LINE__,file=u_FILE_u)
    endif
    if (isPresent) then
       med_methods_FB_FldChk = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function med_methods_FB_FldChk

  !-----------------------------------------------------------------------------

  subroutine med_methods_Field_GetFldPtr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    use ESMF , only : ESMF_Field,ESMF_Mesh, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
    use ESMF , only : ESMF_GeomType_Flag

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(R8), pointer , intent(inout), optional :: fldptr1(:)
    real(R8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)  , optional :: rc

    ! local variables
    type(ESMF_Mesh)          :: lmesh
    integer                  :: lrank, nnodes, nelements
    logical                  :: labort
    type(ESMF_GeomType_Flag) :: geomtype
    character(len=*), parameter :: subname='(med_methods_Field_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call shr_sys_abort(trim(subname)//": ERROR rc not present ", &
            line=__LINE__, file=u_FILE_u)
    endif

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
      labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
      lrank = 0
      if (labort) then
         call shr_sys_abort(trim(subname)//": ERROR data not allocated ")
      else
        call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO)
      endif
    else

      call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_FieldGet(field, rank=lrank, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (geomtype == ESMF_GEOMTYPE_MESH) then
        call ESMF_FieldGet(field, rank=lrank, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (nnodes == 0 .and. nelements == 0) lrank = 0
      else
         call shr_sys_abort(trim(subname)//": ERROR geomtype not supported ")
      endif ! geomtype

      if (lrank == 0) then
         call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
              ESMF_LOGMSG_INFO)

      elseif (lrank == 1) then
        if (.not.present(fldptr1)) then
           call shr_sys_abort(trim(subname)//": ERROR missing rank=1 array ", &
                line=__LINE__, file=u_FILE_u)
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

      elseif (lrank == 2) then
        if (.not.present(fldptr2)) then
           call shr_sys_abort(trim(subname)//": ERROR missing rank=2 array ", &
                line=__LINE__, file=u_FILE_u)
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

      else
         call shr_sys_abort(trim(subname)//": ERROR in rank ", &
              line=__LINE__, file=u_FILE_u)
      endif

    endif  ! status

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Field_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, rank, field, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    ! ----------------------------------------------
    ! Get pointer to a field bundle field
    ! ----------------------------------------------

    type(ESMF_FieldBundle) , intent(in)              :: FB
    character(len=*)       , intent(in)              :: fldname
    real(R8), pointer      , intent(inout), optional :: fldptr1(:)
    real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
    integer                , intent(out),   optional :: rank
    integer                , intent(out),   optional :: rc
    type(ESMF_Field)       , intent(out),   optional :: field

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    character(len=*), parameter :: subname='(med_methods_FB_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call shr_sys_abort(trim(subname)//": ERROR rc not present "//trim(fldname), &
            line=__LINE__, file=u_FILE_u)
    endif

    rc = ESMF_SUCCESS

    if (.not. med_methods_FB_FldChk(FB, trim(fldname), rc=rc)) then

       call shr_sys_abort(trim(subname)//": ERROR field "//trim(fldname)//" not in FB ", &
            line=__LINE__, file=u_FILE_u)
    endif

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call med_methods_Field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) then
      rank = lrank
    endif
    if (present(field)) then
       field = lfield
    endif
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_GetFldPtr

  !-----------------------------------------------------------------------------

  logical function med_methods_FieldPtr_Compare1(fldptr1, fldptr2, cstring, rc)

    ! input/output variables
    real(R8)         , pointer, intent(in)  :: fldptr1(:)
    real(R8)         , pointer, intent(in)  :: fldptr2(:)
    character(len=*) ,          intent(in)  :: cstring
    integer          ,          intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_methods_FieldPtr_Compare1)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    med_methods_FieldPtr_Compare1 = .false.
    if (lbound(fldptr2,1) /= lbound(fldptr1,1) .or. ubound(fldptr2,1) /= ubound(fldptr1,1)) then
       write(msgString,*) trim(subname)//': fldptr1 ',lbound(fldptr1),ubound(fldptr1),"\n",&
            trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2),": ERROR in data size "//trim(cstring)
      call shr_sys_abort(msgstring,rc=rc)
    else
      med_methods_FieldPtr_Compare1 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function med_methods_FieldPtr_Compare1

  !-----------------------------------------------------------------------------

  logical function med_methods_FieldPtr_Compare2(fldptr1, fldptr2, cstring, rc)

    ! input/otuput variables
    real(R8), pointer , intent(in)  :: fldptr1(:,:)
    real(R8), pointer , intent(in)  :: fldptr2(:,:)
    character(len=*)  , intent(in)  :: cstring
    integer           , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_methods_FieldPtr_Compare2)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    med_methods_FieldPtr_Compare2 = .false.
    if (lbound(fldptr2,2) /= lbound(fldptr1,2) .or. lbound(fldptr2,1) /= lbound(fldptr1,1) .or. &
        ubound(fldptr2,2) /= ubound(fldptr1,2) .or. ubound(fldptr2,1) /= ubound(fldptr1,1)) then
      write(msgString,*) trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2),': fldptr1 ',lbound(fldptr1),ubound(fldptr1),&
           ": ERROR in data size "//trim(cstring)
      call shr_sys_abort(trim(msgString))
    else
      med_methods_FieldPtr_Compare2 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function med_methods_FieldPtr_Compare2

  !-----------------------------------------------------------------------------

  subroutine med_methods_State_GeomPrint(state, string, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Field)                :: lfield
    integer                         :: fieldcount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter  :: subname='(med_methods_State_GeomPrint)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldCount > 0) then
       call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lfieldnamelist(fieldCount))
       call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_StateGet(State, itemName=lfieldnamelist(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       deallocate(lfieldnamelist)
       call med_methods_Field_GeomPrint(lfield, string, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_State_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine med_methods_FB_GeomPrint(FB, string, rc)

    use ESMF, only : ESMF_FieldBundle, ESMF_Field, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Field)  :: lfield
    integer           :: fieldcount
    character(len=*),parameter  :: subname='(med_methods_FB_GeomPrint)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldCount > 0) then
      call med_methods_Field_GeomPrint(lfield, string, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_FB_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine med_methods_Field_GeomPrint(field, string, rc)

    use ESMF, only : ESMF_Field, ESMF_Grid, ESMF_Mesh
    use ESMF, only : ESMF_FieldGet, ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_EMPTY
    use ESMF, only : ESMF_GeomType_Flag

    ! input/output variables
    type(ESMF_Field), intent(in)  :: field
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Grid)          :: lgrid
    type(ESMF_Mesh)          :: lmesh
    integer                  :: lrank
    real(R8), pointer        :: dataPtr1(:)
    real(R8), pointer        :: dataPtr2(:,:)
    type(ESMF_GeomType_Flag) :: geomtype
    character(len=*),parameter  :: subname='(med_methods_Field_GeomPrint)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (status == ESMF_FIELDSTATUS_EMPTY) then
       call shr_sys_abort(trim(subname)//":"//trim(string)//": ERROR field does not have a geom yet ", &
            line=__LINE__, file=u_FILE_u)
    endif

    call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, grid=lgrid, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call med_methods_Grid_Print(lgrid, string, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (geomtype == ESMF_GEOMTYPE_MESH) then
      call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call med_methods_Mesh_Print(lmesh, string, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    call med_methods_Field_GetFldPtr(field, &
         fldptr1=dataPtr1, fldptr2=dataPtr2, rank=lrank, abort=.false., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=1 ",lbound(dataptr1,1),ubound(dataptr1,1)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    elseif (lrank == 2) then
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=1 ",lbound(dataptr2,1),ubound(dataptr2,1)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=2 ",lbound(dataptr2,2),ubound(dataptr2,2)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    elseif (lrank == 0) then
      ! means data allocation does not exist yet
      continue
    else
       call shr_sys_abort(trim(subname)//": ERROR rank not supported ", &
            line=__LINE__, file=u_FILE_u)
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Field_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine med_methods_Mesh_Print(mesh, string, rc)

    use ESMF, only: ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF, only: ESMF_DELayoutGet, ESMF_DELayout
    use ESMF, only: ESMF_MeshStatus_Flag, ESMF_MeshStatus_Complete

    type(ESMF_Mesh) , intent(in)  :: mesh
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Distgrid)         :: distgrid
    type(ESMF_DELayout)         :: delayout
    integer                     :: pdim, sdim, nnodes, nelements
    integer                     :: localDeCount
    integer                     :: DeCount
    integer                     :: dimCount, tileCount
    integer, allocatable        :: minIndexPTile(:,:), maxIndexPTile(:,:)
    type(ESMF_MeshStatus_Flag)  :: meshStatus
    logical                     :: elemDGPresent, nodeDGPresent
    character(len=*),parameter  :: subname='(med_methods_Mesh_Print)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh, elementDistGridIsPresent=elemDGPresent, &
         nodalDistgridIsPresent=nodeDGPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(mesh, status=meshStatus, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! first get the distgrid, which should be available
    if (elemDGPresent) then
       call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=element"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (nodeDGPresent) then
       call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=nodal"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (.not. elemDGPresent .and. .not. nodeDGPresent) then
       call ESMF_LogWrite(trim(subname)//": cannot print distgrid from mesh", &
            ESMF_LOGMSG_WARNING, rc=rc)
       return
    endif

    ! if mesh is complete, also get additional parameters
    if (meshStatus==ESMF_MESHSTATUS_COMPLETE) then
       ! access localDeCount to show this is a real Grid
       call ESMF_MeshGet(mesh, parametricDim=pdim, spatialDim=sdim, &
            numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": parametricDim=", pdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
       write (msgString,*) trim(subname)//":"//trim(string)//": spatialDim=", sdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedNodes=", nnodes
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedElements=", nelements
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Mesh_Print

  !-----------------------------------------------------------------------------
  subroutine med_methods_Grid_Print(grid, string, rc)

    use ESMF, only : ESMF_Grid, ESMF_DistGrid, ESMF_StaggerLoc
    use ESMF, only : ESMF_GridGet, ESMF_DistGridGet, ESMF_GridGetCoord
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
    use ESMF, only : ESMF_TypeKind_Flag, ESMF_TYPEKIND_R4, ESMF_TYPEKIND_R8

    ! input/output variabes
    type(ESMF_Grid) , intent(in)  :: grid
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Distgrid)        :: distgrid
    integer                    :: localDeCount
    integer                    :: DeCount
    integer                    :: dimCount, tileCount
    integer                    :: rank
    type(ESMF_StaggerLoc)      :: staggerloc
    type(ESMF_TypeKind_Flag)   :: coordTypeKind
    character(len=32)          :: staggerstr
    integer, allocatable       :: minIndexPTile(:,:), maxIndexPTile(:,:)
    real, pointer              :: fldptrR41D(:)
    real, pointer              :: fldptrR42D(:,:)
    real(R8), pointer          :: fldptrR81D(:)
    real(R8), pointer          :: fldptrR82D(:,:)
    integer                    :: n1,n2,n3
    character(len=*),parameter :: subname='(med_methods_Grid_Print)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! access grid coordinate type
    call ESMF_GridGet(grid, coordTypeKind=coordTypeKind, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! access localDeCount to show this is a real Grid
    call ESMF_GridGet(grid, localDeCount=localDeCount, distgrid=distgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": localDeCount=", localDeCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    ! get dimCount and tileCount
    call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, deCount=deCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write (msgString,*) trim(subname)//":"//trim(string)//": dimCount=", dimCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    write (msgString,*) trim(subname)//":"//trim(string)//": tileCount=", tileCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    write (msgString,*) trim(subname)//":"//trim(string)//": deCount=", deCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
    allocate(minIndexPTile(dimCount, tileCount), &
             maxIndexPTile(dimCount, tileCount))

    ! get minIndex and maxIndex arrays
    call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
       maxIndexPTile=maxIndexPTile, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write (msgString,*) trim(subname)//":"//trim(string)//": minIndexPTile=", minIndexPTile
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    write (msgString,*) trim(subname)//":"//trim(string)//": maxIndexPTile=", maxIndexPTile
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    deallocate(minIndexPTile, maxIndexPTile)

    ! get rank
    call ESMF_GridGet(grid, rank=rank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": rank=", rank
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    do n1 = 1,2
      if (n1 == 1) then
        staggerloc = ESMF_STAGGERLOC_CENTER
        staggerstr = 'ESMF_STAGGERLOC_CENTER'
      elseif (n1 == 2) then
        staggerloc = ESMF_STAGGERLOC_CORNER
        staggerstr = 'ESMF_STAGGERLOC_CORNER'
      else
         call shr_sys_abort(trim(subname)//":staggerloc failure")
      endif
      call ESMF_GridGetCoord(grid, staggerloc=staggerloc, isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      write (msgString,*) trim(subname)//":"//trim(staggerstr)//" present=",isPresent
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
      if (isPresent) then
        do n3 = 0,localDECount-1
        do n2 = 1,dimCount
          if (rank == 1) then
            if (coordTypeKind == ESMF_TYPEKIND_R4) then
              call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptrR41D,rc=rc)
              if (chkerr(rc,__LINE__,u_FILE_u)) return
              write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",&
                   n2,n3,minval(fldptrR41D),maxval(fldptrR41D)
            else if (coordTypeKind == ESMF_TYPEKIND_R8) then
              call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptrR81D,rc=rc)
              if (chkerr(rc,__LINE__,u_FILE_u)) return
              write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",&
                   n2,n3,minval(fldptrR81D),maxval(fldptrR81D)
            else
              write(msgString,*) trim(subname)//":"//" only R4 and R8 types are supported for grid coordinates (1D)!"
            end if
          endif
          if (rank == 2) then
            if (coordTypeKind == ESMF_TYPEKIND_R4) then
              call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptrR42D,rc=rc)
              if (chkerr(rc,__LINE__,u_FILE_u)) return
              write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",&
                   n2,n3,minval(fldptrR42D),maxval(fldptrR42D)
            else if (coordTypeKind == ESMF_TYPEKIND_R8) then
              call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptrR82D,rc=rc)
              if (chkerr(rc,__LINE__,u_FILE_u)) return
              write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",&
                   n2,n3,minval(fldptrR82D),maxval(fldptrR82D)
            else
              write(msgString,*) trim(subname)//":"//" only R4 and R8 types are supported for grid coordinates (2D)!"
            end if

          endif
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
        enddo
        enddo
      endif
    enddo

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Grid_Print

  !-----------------------------------------------------------------------------
  subroutine med_methods_Clock_TimePrint(clock,string,rc)

    use ESMF , only : ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF , only : ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet

    ! input/output variables
    type(ESMF_Clock) , intent(in)          :: clock
    character(len=*) , intent(in),optional :: string
    integer          , intent(out)         :: rc

    ! local variables
    type(ESMF_Time)         :: time
    type(ESMF_TimeInterval) :: timeStep
    character(len=CS)       :: timestr
    character(len=CL)       :: lstring
    character(len=*), parameter :: subname='(med_methods_Clock_TimePrint)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (present(string)) then
      lstring = trim(subname)//":"//trim(string)
    else
      lstring = trim(subname)
    endif

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": currtime = "//trim(timestr), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock,starttime=time,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": startime = "//trim(timestr), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock,stoptime=time,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": stoptime = "//trim(timestr), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock,timestep=timestep,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(timestep,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": timestep = "//trim(timestr), ESMF_LOGMSG_INFO)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_methods_Clock_TimePrint

!================================================================================

  subroutine med_methods_State_GetScalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    use ESMF , only : ESMF_SUCCESS, ESMF_State, ESMF_StateGet, ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LogWrite
    use ESMF , only : ESMF_LOGMSG_INFO, ESMF_VM, ESMF_VMBroadCast, ESMF_VMGetCurrent
    use ESMF , only : ESMF_VMGet

    ! input/output variables
    type(ESMF_State), intent(in)     :: state
    integer,          intent(in)     :: scalar_id
    real(R8),         intent(out)    :: scalar_value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask, icount
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(R8), pointer :: farrayptr(:,:)
    real(r8)          :: tmp(1)
    character(len=*), parameter :: subname='(med_methods_State_GetScalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! check item exist or not?
    call ESMF_StateGet(State, itemSearch=trim(flds_scalar_name), itemCount=icount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (icount > 0) then
      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call shr_sys_abort(trim(subname)//": ERROR in scalar_id", line=__LINE__, file=u_FILE_u)
        endif
        tmp(:) = farrayptr(scalar_id,:)
      endif
      call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      scalar_value = tmp(1)
    else
      scalar_value = 0.0_R8
      call ESMF_LogWrite(trim(subname)//": no ESMF_Field found named: "//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
    end if

  end subroutine med_methods_State_GetScalar

!================================================================================

  subroutine med_methods_State_SetScalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    use ESMF , only : ESMF_Field, ESMF_State, ESMF_StateGet, ESMF_FieldGet
    use ESMF , only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet

    ! input/output arguments
    real(R8),         intent(in)     :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: field
    type(ESMF_VM)     :: vm
    real(R8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(med_methods_State_SetScalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call shr_sys_abort(trim(subname)//": ERROR in scalar_id")
      endif
      farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine med_methods_State_SetScalar

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_getNumFlds(FB, string, nflds, rc)

    ! ----------------------------------------------
    ! Determine if fieldbundle is created and if so, the number of non-scalar
    ! fields in the field bundle
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: string
    integer                , intent(out)   :: nflds
    integer                , intent(inout) :: rc
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    if (.not. ESMF_FieldBundleIsCreated(FB)) then
       call ESMF_LogWrite(trim(string)//": has not been created, returning", ESMF_LOGMSG_INFO)
       nflds = 0
    else
       ! Note - the scalar field has been removed from all mediator
       ! field bundles - so this is why we check if the fieldCount is 0 and not 1 here

       call ESMF_FieldBundleGet(FB, fieldCount=nflds, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (nflds == 0) then
          call ESMF_LogWrite(trim(string)//": only has scalar data, returning", ESMF_LOGMSG_INFO)
       end if
    end if

  end subroutine med_methods_FB_getNumFlds

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_getdata2d(FB, fieldname, dataptr2d, rc)

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fieldname
    real(r8)               , pointer       :: dataptr2d(:,:)
    integer                , intent(inout) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldname=trim(fieldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_methods_FB_getdata2d

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_getdata1d(FB, fieldname, dataptr1d, rc)

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fieldname
    real(r8)               , pointer       :: dataptr1d(:)
    integer                , intent(inout) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldname=trim(fieldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_methods_FB_getdata1d

  !-----------------------------------------------------------------------------
  subroutine med_methods_Field_getdata2d(field, dataptr2d, rc)

    use ESMF, only : ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_Field) , intent(in)    :: field
    real(r8)         , pointer       :: dataptr2d(:,:)
    integer          , intent(inout) :: rc
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, farrayptr=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_methods_Field_getdata2d

  !-----------------------------------------------------------------------------
  subroutine med_methods_Field_getdata1d(field, dataptr1d, rc)

    use ESMF, only : ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_Field) , intent(in)    :: field
    real(r8)         , pointer       :: dataptr1d(:)
    integer          , intent(inout) :: rc
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_methods_Field_getdata1d

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_getmesh(FB, mesh, rc)

    use ESMF, only : ESMF_FieldBundle, ESMF_Field, ESMF_Mesh, ESMF_FieldBundleGet, ESMF_FieldGet

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    type(ESMF_Mesh)        , intent(out)   :: mesh
    integer                , intent(inout) :: rc

    ! local variables
    integer                   :: fieldCount
    type(ESMF_Field), pointer :: fieldlist(:)
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(FB, fieldlist=fieldlist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(fieldlist(1), mesh=mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(fieldlist)

  end subroutine med_methods_FB_getmesh

  !-----------------------------------------------------------------------------
  subroutine med_methods_FB_check_for_nans(FB, maintask, logunit, rc)
    use ESMF, only  : ESMF_FieldBundle, ESMF_Field, ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_LOGMSG_ERROR
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    logical                , intent(in)    :: maintask
    integer                , intent(in)    :: logunit
    integer                , intent(inout) :: rc

    ! local variables
    type(ESMF_Field)            :: field
    integer                     :: index
    integer                     :: fieldcount
    integer                     :: fieldrank
    character(len=CL)           :: fieldname
    real(r8) , pointer          :: dataptr1d(:)
    real(r8) , pointer          :: dataptr2d(:,:)
    integer                     :: nancount
    character(len=CS)           :: nancount_char
    character(len=CL)           :: msg_error
    logical                     :: nanfound
    character(len=*), parameter :: subname='(med_methods_FB_check_for_nans)'
    ! ----------------------------------------------
    rc = ESMF_SUCCESS

    if(.not. mediator_checkfornans) return
    
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    nanfound = .false.
    do index=1,fieldCount
       call med_methods_FB_getNameN(FB, index, fieldname, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleGet(FB, fieldName=fieldname, field=field, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, rank=fieldrank, name=fieldname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (fieldrank == 1) then
          call ESMF_FieldGet(field, farrayPtr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call med_methods_check_for_nans(dataptr1d, nancount)
       else
          call ESMF_FieldGet(field, farrayPtr=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call med_methods_check_for_nans(dataptr2d, nancount)
       end if
       if (nancount > 0) then
          write(nancount_char, '(i0)') nancount
          msg_error = "ERROR: " // trim(nancount_char) //" nans found in "//trim(fieldname)
          nanfound = .true.
          call ESMF_LogWrite(trim(msg_error), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       end if
    end do
    if (nanfound) then
       call shr_sys_abort('ABORTING JOB, see PET file for details', line=__LINE__, file=u_FILE_u)
    end if
    
  end subroutine med_methods_FB_check_for_nans

  !-----------------------------------------------------------------------------
  subroutine med_methods_check_for_nans_1d(dataptr, nancount)
    use shr_infnan_mod, only: shr_infnan_isnan
    ! input/output variables
    real(r8) , intent(in)  :: dataptr(:)
    integer  , intent(out) :: nancount
    ! local variables
    integer :: n

    nancount = 0
    do n = 1,size(dataptr)
       if (shr_infnan_isnan(dataptr(n))) then
          nancount = nancount + 1
       end if
    end do
  end subroutine med_methods_check_for_nans_1d

  subroutine med_methods_check_for_nans_2d(dataptr, nancount)
    use shr_infnan_mod, only: shr_infnan_isnan
    ! input/output variables
    real(r8) , intent(in)  :: dataptr(:,:)
    integer  , intent(out) :: nancount
    ! local variables
    integer :: n,k

    nancount = 0
    do k = 1,size(dataptr, dim=1)
       do n = 1,size(dataptr, dim=2)
          if (shr_infnan_isnan(dataptr(k,n))) then
             nancount = nancount + 1
          end if
       end do
    end do
  end subroutine med_methods_check_for_nans_2d

end module med_methods_mod
