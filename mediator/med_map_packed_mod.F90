module med_map_packed_mod

  !------------------------------------------
  ! Create packed field bundles needed for mapping
  ! and map the packed field bundles 
  !------------------------------------------

  use med_kind_mod      , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, i8=>shr_kind_i8, r8=>shr_kind_r8
  use med_constants_mod , only : dbug_flag => med_constants_dbug_flag
  use med_utils_mod     , only : chkerr    => med_utils_ChkErr

  implicit none
  private

  public :: med_map_packed_fieldbundles_create
  public :: med_map_packed_fieldbundles

  integer    , parameter :: number_strlen = 2
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains  
!===============================================================================

  subroutine med_map_packed_fieldbundles_create(flds_scalar_name, fldsSrc, &
       FBSrc, FBSrc_packed, FBDst_packed, FBDst, rc)
    
    use ESMF
    use esmFlds         , only : mapfcopy, med_fldList_entry_type
    use med_methods_mod , only : FB_getFieldN => med_methods_FB_getFieldN

    ! input/output variables
    character(len=*)             , intent(in)    :: flds_scalar_name
    type(med_fldList_entry_type) , pointer       :: fldsSrc(:)
    type(ESMF_FieldBundle)       , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)       , intent(inout) :: FBSrc_packed(:) ! array over mapping types
    type(ESMF_FieldBundle)       , intent(inout) :: FBDst_packed(:) ! array over mapping types   
    type(ESMF_FieldBundle)       , intent(in)    :: FBDst
    integer                      , intent(out)   :: rc

    ! local variables
    integer                    :: n, n1, nf, nu
    integer                    :: npacked  
    integer                    :: fieldcount
    integer                    :: lrank
    type(ESMF_Field)           :: lfield
    character(CL)              :: fldname
    integer                    :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    integer                    :: gridToFieldMap(1)      ! currently the size must equal 1 for rank 2 fields
    character(CL), allocatable :: lfieldnamelist(:)
    character(CL), allocatable :: lfieldnamelist_packed(:)
    real(r8), pointer          :: ptrsrc_packed(:,:)
    real(r8), pointer          :: ptrdst_packed(:,:)
    integer                    :: lsize_src
    integer                    :: lsize_dst
    type(ESMF_Mesh)            :: lmesh_src
    type(ESMF_Mesh)            :: lmesh_dst
    integer                    :: mapindex
    character(len=number_strlen) :: cnumber
    character(len=*), parameter  :: subname=' (med_packed_fieldbundles_create) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: Just for now assume redistribution, need to leverage fldssrc which is simply not being used
    mapindex = mapfcopy
    
    ! determine field list in source field bundle
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldCount
       if (trim(lfieldnamelist(n)) == trim(flds_scalar_name) .or. trim(lfieldnamelist(n)) == '') then
          do n1 = n, fieldCount-1
             lfieldnamelist(n1) = lfieldnamelist(n1+1)
          enddo
          fieldCount = fieldCount - 1
       endif
    enddo  ! n

    ! Reset fieldcount based on the fact that some fields have ungridded dimensions and
    ! need to unwrap them into separate fields for the purposes of packing
    npacked = 0
    do nf = 1, fieldCount
       fldname = trim(lfieldnamelist(nf))
       call ESMF_FieldBundleGet(FBSrc, fieldName=trim(fldname), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (lrank == 2) then
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          npacked = npacked + ungriddedUBound(1)
       else
          npacked = npacked + 1
       end if
    end do

    ! Get list of fields in including new field names for each ungridded dimension
    npacked = 0
    allocate(lfieldnamelist_packed(npacked))
    do nf = 1, fieldCount
       fldname = trim(lfieldnamelist(nf))
       call ESMF_FieldBundleGet(FBSrc, fieldName=trim(fldname), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (lrank == 2) then
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nu = 1,ungriddedUBound(1)
             npacked = npacked + 1
             write(cnumber,'(i0)') nu
             lfieldnamelist_packed(npacked) = trim(fldname)//trim(cnumber)
          end do
       else
          npacked = npacked + 1
          lfieldnamelist_packed(npacked) = trim(fldname)
       end if
    end do

    ! Determine local size and mesh of source fields
    ! Allocate a source fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call FB_getFieldN(FBSrc, 1, lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, numOwnedElements=lsize_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ptrsrc_packed(npacked, lsize_src))
    FBSrc_packed(mapindex) = ESMF_FieldBundleCreate(lfieldnamelist_packed, &
         ptrsrc_packed, lmesh_src, lsize_src, gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

    ! Determine local size of destination fields
    ! Allocate a destination fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call FB_getFieldN(FBDst, 1, lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, numOwnedElements=lsize_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ptrdst_packed(npacked, lsize_dst))
    FBDst_packed(mapindex) = ESMF_FieldBundleCreate(lfieldnamelist_packed, &
         ptrdst_packed, lmesh_dst, lsize_dst, gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

    deallocate(lfieldnamelist)
    deallocate(lfieldnamelist_packed)

  end subroutine med_map_packed_fieldbundles_create

  !================================================================================
  subroutine med_map_packed_fieldbundles(flds_scalar_name, FBSrc, FBSrc_packed, FBDst_packed, FBDst, &
       RouteHandles, rc)
    
    ! -----------------------------------------------
    ! Do the redistribution on the packed field bundles
    ! -----------------------------------------------

    use ESMF
    use esmFlds         , only : mapfcopy
    use med_methods_mod , only : FB_getFieldByName => med_methods_FB_getFieldByName

    ! input/output variables
    character(len=*)       , intent(in)    :: flds_scalar_name
    type(ESMF_FieldBundle) , intent(in)    :: FBSrc
    type(ESMF_FieldBundle) , intent(inout) :: FBSrc_packed(:) ! array over mapping types
    type(ESMF_FieldBundle) , intent(inout) :: FBDst_packed(:) ! array over mapping types
    type(ESMF_FieldBundle) , intent(inout) :: FBDst
    type(ESMF_RouteHandle) , intent(inout) :: RouteHandles(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer                    :: n, n1, nf, nu
    integer                    :: npacked  
    integer                    :: fieldcount
    character(CL), allocatable :: lfieldnamelist(:)
    integer                    :: lrank
    type(ESMF_Field)           :: lfield
    type(ESMF_Field)           :: lfield_packed
    character(CL)              :: fldname
    integer                    :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer          :: dataptr1d(:)
    real(r8), pointer          :: dataptr2d(:,:)
    real(r8), pointer          :: dataptr1d_packed(:)
    integer                    :: lsize_src
    integer                    :: lsize_dst
    integer                    :: mapindex
    character(len=number_strlen) :: cnumber
    character(len=*), parameter  :: subname=' (med_map_packed_fieldbundles) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: Just for now assume redistribution, need to leverage fldssrc which is simply not being used
    mapindex = mapfcopy

    ! determine field list in source field bundle
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldCount
       if (trim(lfieldnamelist(n)) == trim(flds_scalar_name) .or. trim(lfieldnamelist(n)) == '') then
          do n1 = n, fieldCount-1
             lfieldnamelist(n1) = lfieldnamelist(n1+1)
          enddo
          fieldCount = fieldCount - 1
       endif
    enddo  ! n

    ! copy the src fields into the packed field bundle
    do n = 1,fieldcount
       fldname = lfieldnamelist(n)
       call FB_getFieldByName(FBSrc, trim(fldname), lfield, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (lrank == 2) then
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nu = 1,ungriddedUBound(1)
             ! Name of packed field is same as name of src field with
             ! the ungridded dimension appended
             write(cnumber,'(i0)') nu
             call FB_getFieldByName(FbSrc_packed(mapindex), trim(fldname)//trim(cnumber), lfield_packed, rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_packed, farrayptr=dataptr1d_packed, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr1d_packed(:) = dataptr2d(nu,:)
          end do
       else
          ! Name of packed field is same as name of src field
          call FB_getFieldByName(FBSrc_packed(mapindex), trim(fldname), lfield_packed, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_packed, farrayptr=dataptr1d_packed, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_packed(:) = dataptr1d(:)
       end if
    end do

    ! Do the mapping
    call ESMF_FieldBundleRedist(FBSrc_packed(mapindex), FBDst_packed(mapindex), Routehandles(mapindex), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Copy the destination packed field bundle into the destination unpacked field bundle
    do n = 1,fieldcount
       fldname = lfieldnamelist(n)
       call FB_getFieldByName(FBDst, trim(fldname), lfield, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (lrank == 2) then
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nu = 1,ungriddedUBound(1)
             ! Name of packed field is same as name of field with
             ! the ungridded dimension appended
             write(cnumber,'(i0)') nu
             call FB_getFieldByName(FBDst_packed(mapindex), &
                  trim(fldname)//trim(cnumber), lfield_packed, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_packed, farrayptr=dataPtr1d_packed, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayptr=dataPtr2d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr2d(nu,:) = dataptr1d_packed(:)
          end do
       else
          ! Name of packed field is same as name of src field
          call FB_getFieldByName(FBDst_packed(mapindex), trim(fldname), lfield_packed, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_packed, farrayptr=dataptr1d_packed, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d(:) = dataptr1d_packed(:)
       end if
    end do

    deallocate(lfieldnamelist)

  end subroutine med_map_packed_fieldbundles

end module med_map_packed_mod
