module med_map_packed_mod

  !------------------------------------------
  ! Create packed field bundles needed for mapping
  ! and map the packed field bundles 
  !------------------------------------------

  use med_kind_mod      , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, i8=>shr_kind_i8, r8=>shr_kind_r8
  use med_constants_mod , only : dbug_flag => med_constants_dbug_flag
  use med_utils_mod     , only : chkerr    => med_utils_ChkErr
  use perf_mod          , only : t_startf, t_stopf

  implicit none
  private

  public :: med_map_packed_fieldbundles_create
  public :: med_map_packed_fieldbundles

  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains  
!===============================================================================

  subroutine med_map_packed_fieldbundles_create(fldsSrc, FBSrc, FBDst, fieldsrc_packed, fielddst_packed, rc)
    
    use ESMF
    use esmFlds, only : mapfcopy, med_fldList_entry_type

    ! input/output variables
    type(med_fldList_entry_type) , pointer       :: fldsSrc(:)
    type(ESMF_FieldBundle)       , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)       , intent(inout) :: FBDst
    type(ESMF_Field)             , intent(inout) :: fieldsrc_packed(:) ! array over mapping types
    type(ESMF_Field)             , intent(inout) :: fielddst_packed(:) ! array over mapping types   
    integer                      , intent(out)   :: rc

    ! local variables
    integer                    :: nf, nu
    integer                    :: npacked  
    integer                    :: fieldcount
    type(ESMF_Field)           :: lfield
    integer                    :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer          :: ptrsrc_packed(:,:)
    real(r8), pointer          :: ptrdst_packed(:,:)
    integer                    :: lsize_src
    integer                    :: lsize_dst
    type(ESMF_Mesh)            :: lmesh_src
    type(ESMF_Mesh)            :: lmesh_dst
    integer                    :: mapindex
    type(ESMF_Field), pointer  :: fieldlist_src(:)     
    type(ESMF_Field), pointer  :: fieldlist_dst(:)     
    character(len=*), parameter  :: subname=' (med_packed_fieldbundles_create) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: Just for now assume redistribution, need to leverage fldssrc which is simply not being used
    mapindex = mapfcopy
    
    ! Get field count for both FBsrc and FBdst
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_src(fieldcount))
    allocate(fieldlist_dst(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldlist=fieldlist_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBdst, fieldlist=fieldlist_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    ! Calculated size of packed field based on the fact that some fields have ungridded dimensions and
    ! need to unwrap them into separate fields for the purposes of packing
    npacked = 0
    do nf = 1, fieldCount
       call ESMF_FieldGet(fieldlist_src(nf), ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (ungriddedUBound(1) > 0) then
          do nu = 1,ungriddedUBound(1)
             npacked = npacked + 1
          end do
       else
          npacked = npacked + 1
       end if
    end do

    ! Determine local size and mesh of source fields
    ! Allocate a source fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call ESMF_FieldGet(fieldlist_src(1), mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, numOwnedElements=lsize_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ptrsrc_packed(npacked, lsize_src))
    fieldsrc_packed(mapindex) = ESMF_FieldCreate(lmesh_src, &
         ptrsrc_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

    ! Determine local size of destination fields
    ! Allocate a destination fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call ESMF_FieldGet(fieldlist_dst(1), mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, numOwnedElements=lsize_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ptrdst_packed(npacked, lsize_dst))
    fielddst_packed(mapindex) = ESMF_FieldCreate(lmesh_dst, &
         ptrdst_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

  end subroutine med_map_packed_fieldbundles_create

  !================================================================================
  subroutine med_map_packed_fieldbundles(FBSrc, FBDst, fieldsrc_packed, fielddst_packed, &
       routehandles, rc)
    
    ! -----------------------------------------------
    ! Do the redistribution on the packed field bundles
    ! -----------------------------------------------

    use ESMF
    use esmFlds         , only : mapfcopy
    use med_methods_mod , only : FB_getFieldN => med_methods_FB_getFieldN

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FBSrc
    type(ESMF_FieldBundle) , intent(inout) :: FBDst
    type(ESMF_Field)       , intent(inout) :: fieldsrc_packed(:) ! array over mapping types
    type(ESMF_Field)       , intent(inout) :: fielddst_packed(:) ! array over mapping types
    type(ESMF_RouteHandle) , intent(inout) :: routehandles(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer                   :: nf, nu
    integer                   :: npacked  
    integer                   :: fieldcount
    type(ESMF_Field)          :: lfield
    integer                   :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer         :: dataptr1d(:)
    real(r8), pointer         :: dataptr2d(:,:)
    real(r8), pointer         :: dataptr2d_packed(:,:)
    integer                   :: mapindex
    type(ESMF_Field), pointer :: fieldlist_src(:)     
    type(ESMF_Field), pointer :: fieldlist_dst(:)     
    character(len=*), parameter  :: subname=' (med_map_packed_fieldbundles) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: Just for now assume redistribution, need to leverage fldssrc which is simply not being used
    mapindex = mapfcopy

    ! Get field count for both FBsrc and FBdst
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_src(fieldcount))
    allocate(fieldlist_dst(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldlist=fieldlist_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBdst, fieldlist=fieldlist_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the pointer for the packed source data
    call ESMF_FieldGet(fieldsrc_packed(mapindex), farrayptr=dataptr2d_packed, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Copy the src fields into the packed field bundle
    call t_startf('MED:'//trim(subname)//' copy from src')
    npacked = 0
    do nf = 1,fieldcount
       call ESMF_FieldGet(fieldlist_src(nf), ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (ungriddedUBound(1) > 0) then
          call ESMF_FieldGet(fieldlist_src(nf), farrayptr=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nu = 1,ungriddedUBound(1)
             npacked = npacked + 1
             dataptr2d_packed(npacked,:) = dataptr2d(nu,:)
          end do
       else
          call ESMF_FieldGet(fieldlist_src(nf), farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          npacked = npacked + 1
          dataptr2d_packed(npacked,:) = dataptr1d(:)
       end if
    end do
    call t_stopf('MED:'//trim(subname)//' copy from src')

    ! Do the mapping
    call t_startf('MED:'//trim(subname)//' map')
    call ESMF_FieldRedist(fieldsrc_packed(mapindex), fielddst_packed(mapindex), routehandles(mapindex), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//trim(subname)//' map')

    ! Get the pointer for the packed destination data
    call ESMF_FieldGet(fielddst_packed(mapindex), farrayptr=dataptr2d_packed, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Copy the destination packed field bundle into the destination unpacked field bundle
    call t_startf('MED:'//trim(subname)//' copy to dest')
    npacked = 0
    do nf = 1,fieldcount
       call ESMF_FieldGet(fieldlist_dst(nf), ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (ungriddedUBound(1) > 0) then
          call ESMF_FieldGet(fieldlist_dst(nf), farrayptr=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nu = 1,ungriddedUBound(1)
             npacked = npacked + 1
             dataptr2d(nu,:) = dataptr2d_packed(npacked,:)
          end do
       else
          call ESMF_FieldGet(fieldlist_dst(nf), farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          npacked = npacked + 1
          dataptr1d(:) = dataptr2d_packed(npacked,:)
       end if
    end do
    call t_stopf('MED:'//trim(subname)//' copy to dest')

    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

  end subroutine med_map_packed_fieldbundles

end module med_map_packed_mod
