module med_map_packed_mod

  !------------------------------------------
  ! Create packed field bundles needed for mapping
  ! and map the packed field bundles
  !------------------------------------------

  use med_kind_mod      , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, i8=>shr_kind_i8, r8=>shr_kind_r8
  use med_constants_mod , only : dbug_flag => med_constants_dbug_flag
  use med_utils_mod     , only : chkerr    => med_utils_ChkErr
  use perf_mod          , only : t_startf, t_stopf
  use med_internalstate_mod, only : mastertask, logunit

  use shr_sys_mod, only : shr_sys_flush

  implicit none
  private

  public :: med_map_packed_field_create
  public :: med_map_packed_field_map

  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_map_packed_field_create(destcomp, flds_scalar_name, &
       fldsSrc, FBSrc, FBDst, packed_data, rc)

    use ESMF
    use esmFlds               , only : med_fldList_entry_type, nmappers
    use esmFlds               , only : ncomps, compatm, compice, compocn, compname, mapnames
    use med_internalstate_mod , only : packed_data_type

    ! input/output variables
    integer                      , intent(in)    :: destcomp
    character(len=*)             , intent(in)    :: flds_scalar_name 
    type(med_fldList_entry_type) , pointer       :: fldsSrc(:) ! array over mapping types 
    type(ESMF_FieldBundle)       , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)       , intent(inout) :: FBDst
    type(packed_data_type)       , intent(inout) :: packed_data(:) ! array over mapping types
    integer                      , intent(out)   :: rc

    ! local variables
    integer                    :: nf, nu, ns
    integer, allocatable       :: npacked(:)
    integer                    :: fieldcount
    type(ESMF_Field)           :: lfield
    integer                    :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer          :: ptrsrc_packed(:,:) => null()
    real(r8), pointer          :: ptrdst_packed(:,:) => null()
    integer                    :: lsize_src
    integer                    :: lsize_dst
    type(ESMF_Mesh)            :: lmesh_src
    type(ESMF_Mesh)            :: lmesh_dst
    integer                    :: mapindex
    type(ESMF_Field), pointer  :: fieldlist_src(:) => null()
    type(ESMF_Field), pointer  :: fieldlist_dst(:) => null()
    character(CL), allocatable :: fieldNameList(:)
    character(len=*), parameter  :: subname=' (med_packed_fieldbundles_create) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get field count for both FBsrc and FBdst
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get fields in source and destination field bundles
    allocate(fieldlist_src(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldlist=fieldlist_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_dst(fieldcount))
    call ESMF_FieldBundleGet(FBdst, fieldlist=fieldlist_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! field names are the same for the source and destination field bundles
    allocate(fieldnamelist(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldnamelist=fieldnamelist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine local size and mesh of source fields
    ! Allocate a source fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call ESMF_FieldGet(fieldlist_src(1), mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, numOwnedElements=lsize_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine local size of destination fields
    ! Allocate a destination fortran pointer for the new packed field bundle
    ! Create the packed source field bundle
    call ESMF_FieldGet(fieldlist_dst(1), mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, numOwnedElements=lsize_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Gather all fields that will be mapped with a target map index into a packed field
    ! Calculated size of packed field based on the fact that some fields have 
    ! ungridded dimensions and need to unwrap them into separate fields for the 
    ! purposes of packing

    ! Allocate memory to keep tracked of packing index for each mapping type
    allocate(npacked(nmappers))
    npacked(:) = 0

    if (mastertask) write(logunit,*)
    ! Loop over mapping types
    do mapindex = 1,nmappers

       ! Allocate the fldindex attribute of packed_indices if needed
       if (.not. allocated(packed_data(mapindex)%fldindex)) then
          allocate(packed_data(mapindex)%fldindex(fieldcount))
          packed_data(mapindex)%fldindex(:) = -999
       end if

       ! Loop over the fields in FBSrc
       do nf = 1, fieldCount

          ! Loop over the fldsSrc types
          do ns = 1,size(fldsSrc)

             if ( fldsSrc(ns)%mapindex(destcomp) == mapindex .and. &
                  fldsSrc(ns)%shortname /= flds_scalar_name  .and. &
                  trim(fldsSrc(ns)%shortname) == trim(fieldnamelist(nf))) then

                ! Determine mapnorm - the assumption is that there is only one
                ! mapnorm type for a source packed field bundle
                if (npacked(mapindex) == 0) then
                   packed_data(mapindex)%mapnorm = fldsSrc(ns)%mapnorm(destcomp)
                else
                   if (fldsSrc(ns)%mapnorm(destcomp) /= packed_data(mapindex)%mapnorm) then
                      call ESMF_LogWrite(trim(subname)//&
                           ": ERROR mapnorm "//trim(fldsSrc(ns)%mapnorm(destcomp))//" is not equal to "// &
                           trim(packed_data(mapindex)%mapnorm)//" for packed field bundle ", &  
                           ESMF_LOGMSG_ERROR, rc=rc)
                      rc = ESMF_FAILURE
                      return
                   end if
                end if

                ! Determine mapping of indices into packed field bundle
                ! Get source field
                call ESMF_FieldGet(fieldlist_src(nf), ungriddedUBound=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (ungriddedUBound(1) > 0) then
                   do nu = 1,ungriddedUBound(1)
                      npacked(mapindex) = npacked(mapindex) + 1
                      if (nu == 1) then
                         packed_data(mapindex)%fldindex(nf) = npacked(mapindex)
                      end if
                   end do
                else
                   npacked(mapindex) = npacked(mapindex) + 1
                   packed_data(mapindex)%fldindex(nf) = npacked(mapindex)
                end if

                if (mastertask) then
                   write(logunit,'(5(a,2x),2x,i4)') trim(subname)//&
                        'Packed field: destcomp,mapping,mapnorm,fldname,index: ', &
                        trim(compname(destcomp)), &
                        trim(mapnames(mapindex)), &
                        trim(fldsSrc(ns)%mapnorm(destcomp)), &
                        trim(fieldnamelist(nf)), &
                        packed_data(mapindex)%fldindex(nf)
                end if

             end if! end if source field is mapped to destination field with mapindex
          end do ! end loop over FBSrc fields
       end do ! end loop over fldSrc elements

       ! Create the source and destination packed fields for mapindex
       if (npacked(mapindex) > 0) then
          allocate(ptrsrc_packed(npacked(mapindex), lsize_src))
          packed_data(mapindex)%field_src = ESMF_FieldCreate(lmesh_src, &
               ptrsrc_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          allocate(ptrdst_packed(npacked(mapindex), lsize_dst))
          packed_data(mapindex)%field_dst = ESMF_FieldCreate(lmesh_dst, &
               ptrdst_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          packed_data(mapindex)%field_fracsrc = ESMF_FieldCreate(lmesh_src, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          packed_data(mapindex)%field_fracdst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       end if
    end do ! end loop over mapindex

    deallocate(npacked) 
    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

  end subroutine med_map_packed_field_create

  !================================================================================
  subroutine med_map_packed_field_map(FBSrc, FBDst, FBFracSrc, FBNormOne, &
       packed_data, routehandles, rc)

    ! -----------------------------------------------
    ! Do the redistribution on the packed field bundles
    ! -----------------------------------------------

    use ESMF
    use esmFlds               , only : mapunset, mapnames, nmappers
    use esmFlds               , only : mapnstod, mapnstod_consd, mapnstod_consf, mapnstod_consd
    use esmFlds               , only : ncomps, compatm, compice, compocn, compname
    use esmFlds               , only : mapfcopy, mapconsd, mapconsf, mapnstod
    use med_map_mod           , only : med_map_field_regrid
    use med_internalstate_mod , only : packed_data_type

    ! input/output variables
    type(ESMF_FieldBundle)    , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)    , intent(inout) :: FBDst
    type(ESMF_FieldBundle)    , intent(in)    :: FBNormOne(:)   ! array over mapping types 
    type(ESMF_FieldBundle)    , intent(in)    :: FBFracSrc      ! fraction field bundle for source 
    type(packed_data_type)    , intent(inout) :: packed_data(:) ! array over mapping types
    type(ESMF_RouteHandle)    , intent(inout) :: routehandles(:)
    integer                   , intent(out)   :: rc

    ! local variables
    integer                    :: nf, nu, np, n
    integer                    :: fieldcount
    integer                    :: mapindex
    integer                    :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer          :: dataptr1d(:) => null()
    real(r8), pointer          :: dataptr2d(:,:) => null()
    real(r8), pointer          :: dataptr2d_packed(:,:) => null()
    type(ESMF_Field)           :: lfield
    type(ESMF_Field)           :: frac_field_src
    type(ESMF_Field), pointer  :: fieldlist_src(:) => null()
    type(ESMF_Field), pointer  :: fieldlist_dst(:) => null()
    character(CL), allocatable :: fieldNameList(:)
    real(r8), pointer          :: data_src(:,:) => null()
    real(r8), pointer          :: data_srctmp(:,:) => null()
    real(r8), pointer          :: data_dst(:,:) => null()
    real(r8), pointer          :: data_fracsrc(:) => null()
    real(r8), pointer          :: data_fracdst(:) => null()
    real(r8), pointer          :: data_norm(:) => null()
    character(len=*), parameter  :: subname=' (med_map_packed_fieldbundles) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get field count for both FBsrc and FBdst
    call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldlist_src(fieldcount))
    call ESMF_FieldBundleGet(FBsrc, fieldlist=fieldlist_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldlist_dst(fieldcount))
    call ESMF_FieldBundleGet(FBdst, fieldlist=fieldlist_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Loop over mapping types
    do mapindex = 1,nmappers

       ! If packed field is created
       if (ESMF_FieldIsCreated(packed_data(mapindex)%field_src)) then

          ! -----------------------------------
          ! Copy the src fields into the packed field bundle
          ! -----------------------------------

          call t_startf('MED:'//trim(subname)//' copy from src')

          ! First get the pointer for the packed source data
          call ESMF_FieldGet(packed_data(mapindex)%field_src, farrayptr=dataptr2d_packed, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Now do the copy
          do nf = 1,fieldcount
             ! Get the indices into the packed data structure
             np = packed_data(mapindex)%fldindex(nf)
             if (np > 0) then
                call ESMF_FieldGet(fieldlist_src(nf), ungriddedUBound=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (ungriddedUBound(1) > 0) then
                   call ESMF_FieldGet(fieldlist_src(nf), farrayptr=dataptr2d, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   do nu = 1,ungriddedUBound(1)
                      dataptr2d_packed(np+nu-1,:) = dataptr2d(nu,:)
                   end do
                else
                   call ESMF_FieldGet(fieldlist_src(nf), farrayptr=dataptr1d, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   dataptr2d_packed(np,:) = dataptr1d(:)
                end if
             end if
          end do
          call t_stopf('MED:'//trim(subname)//' copy from src')

          ! -----------------------------------
          ! Do the mapping
          ! -----------------------------------

          call t_startf('MED:'//trim(subname)//' map')
          if (mapindex == mapfcopy) then

             ! Mapping is redistribution
             call ESMF_FieldRedist(packed_data(mapindex)%field_src, packed_data(mapindex)%field_dst, &
                  routehandles(mapindex), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          else if ( trim(packed_data(mapindex)%mapnorm) /= 'unset' .and. &
                    trim(packed_data(mapindex)%mapnorm) /= 'one'   .and. &
                    trim(packed_data(mapindex)%mapnorm) /= 'none') then
             
             ! Mapping is not redistribution
             ! ASSUME that  each packed field has only one normalization type
             ! normalize packed source data 
             ! - get a pointer (data_fracsrc) to the normalization array
             ! - get a pointer (data_src) to source field data in FBSrc
             ! - copy data_src to data_srctmp
             ! - normalize data_src by data_fracsrc
             call ESMF_FieldBundleGet(FBFracSrc, fieldName=packed_data(mapindex)%mapnorm, field=lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayPtr=data_fracsrc, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(packed_data(mapindex)%field_src, farrayPtr=data_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(data_srctmp(size(data_src,dim=1), size(data_src,dim=2)))
             data_srctmp(:,:) = data_src(:,:)
             do n = 1,size(data_fracsrc)
                data_src(:,n) = data_src(:,n) * data_fracsrc(n)
             end do

             ! regrid normalized packed source field
             call med_map_field_regrid (packed_data(mapindex)%field_src, packed_data(mapindex)%field_dst, &
                  routehandles, mapindex, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! restore original value to packed source field
             data_src(:,:) = data_srctmp(:,:)
             deallocate(data_srctmp)

             call t_startf('MED:'//trim(subname)//' map_nofrac')
             ! regrid fraction field from source to destination
             call ESMF_FieldBundleGet(FBFracSrc, fieldname=trim(packed_data(mapindex)%mapnorm), &
                  field=frac_field_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call med_map_field_regrid(frac_field_src, packed_data(mapindex)%field_fracdst, &
                  routehandles, mapindex, rc=rc)
             call t_stopf('MED:'//trim(subname)//' map_nofrac')

             call t_startf('MED:'//trim(subname)//' norm one')
             ! get pointer to mapped fraction and normalize
             ! destination mapped values by the reciprocal of the mapped fraction
             call ESMF_FieldGet(packed_data(mapindex)%field_fracdst, farrayPtr=data_fracdst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(packed_data(mapindex)%field_dst, farrayPtr=data_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_dst,dim=2)
                if (data_fracdst(n) == 0.0_r8) then
                   data_dst(:,n) = 0.0_r8
                else
                   data_dst(:,n) = data_dst(:,n)/data_fracdst(n)
                end if
             end do
             call t_stopf('MED:'//trim(subname)//' norm one')

          else if ( trim(packed_data(mapindex)%mapnorm) == 'one' .or. trim(packed_data(mapindex)%mapnorm) == 'none') then

             ! Mapping is not redistribution
             call med_map_field_regrid (packed_data(mapindex)%field_src, packed_data(mapindex)%field_dst, &
                  routehandles, mapindex, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! obtain unity normalization factor and multiply
             ! interpolated field by reciprocal of normalization factor
             if (trim(packed_data(mapindex)%mapnorm) == 'one') then
                call ESMF_FieldBundleGet(FBNormOne(mapindex), fieldName='one', field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(lfield, farrayPtr=data_norm, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(packed_data(mapindex)%field_dst, farrayPtr=data_dst, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                do n = 1,size(data_dst,dim=2)
                   if (data_norm(n) == 0.0_r8) then
                      data_dst(:,n) = 0.0_r8
                   else
                      data_dst(:,n) = data_dst(:,n)/data_norm(n)
                   end if
                end do
             end if

          end if 
          call t_stopf('MED:'//trim(subname)//' map')

          ! -----------------------------------
          ! Copy the destination packed field bundle into the destination unpacked field bundle
          ! -----------------------------------

          call t_startf('MED:'//trim(subname)//' copy to dest')

          ! First get the pointer for the packed destination data
          call ESMF_FieldGet(packed_data(mapindex)%field_dst, farrayptr=dataptr2d_packed, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Now do the copy back to FBDst
          do nf = 1,fieldcount
             ! Get the indices into the packed data structure
             np = packed_data(mapindex)%fldindex(nf)
             if (np > 0) then
                call ESMF_FieldGet(fieldlist_dst(nf), ungriddedUBound=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (ungriddedUBound(1) > 0) then
                   call ESMF_FieldGet(fieldlist_dst(nf), farrayptr=dataptr2d, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   do nu = 1,ungriddedUBound(1)
                      dataptr2d(nu,:) = dataptr2d_packed(np+nu-1,:)
                   end do
                else
                   call ESMF_FieldGet(fieldlist_dst(nf), farrayptr=dataptr1d, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   dataptr1d(:) = dataptr2d_packed(np,:)
                end if
             end if
          end do
          call t_stopf('MED:'//trim(subname)//' copy to dest')

       end if
    end do ! end of loop over mapindex

    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

  end subroutine med_map_packed_field_map

end module med_map_packed_mod
