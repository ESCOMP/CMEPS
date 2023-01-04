module med_map_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_kind_mod          , only : I4=>SHR_KIND_I4
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_Field
  use med_internalstate_mod , only : InternalState, logunit, mastertask
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr    => med_utils_ChkErr
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  ! public routines
  public :: med_map_routehandles_init
  public :: med_map_rh_is_created
  public :: med_map_packed_field_create
  public :: med_map_field_packed
  public :: med_map_field_normalized
  public :: med_map_field

  interface med_map_routehandles_init
     module procedure med_map_routehandles_initfrom_esmflds  ! called from med.F90
     module procedure med_map_routehandles_initfrom_fieldbundle
     module procedure med_map_routehandles_initfrom_field
  end interface

  interface med_map_RH_is_created
     module procedure med_map_RH_is_created_RH3d
     module procedure med_map_RH_is_created_RH1d
  end interface

  type(ESMF_Field) :: uv3d_src, uv3d_dst ! needed for 3d mapping of u,v vector pairs

  ! private module variables

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_map_RouteHandles_initfrom_esmflds(gcomp, flds_scalar_name, llogunit, rc)

    !---------------------------------------------
    ! Initialize route handles in the mediator and also
    ! initialize unity normalization fields and do the mapping for
    ! unity normalization up front
    !
    ! Assumptions:
    !   -  Route handles are created per target field bundles NOT
    !      per individual fields in the bundle
    !   -  ALL fields in the bundle are on identical grids
    !   -  MULTIPLE route handles are going to be generated for
    !      given field bundle source and destination grids
    !    - Route handles will ONLY be created if coupling_active is true between n1 and n2
    ! Algorithm
    !     n1=source component index
    !     n2=destination component index
    !     nf=field index
    !     fldListFr(n)%flds(nf) is queried to determine the mapindex and mapfile
    !     and the appropriate route handle is created
    !
    ! Regridding is done on a per-field basis AND only for those fields that have a
    ! valid mapping index for the destination component
    !     n = source field index field index
    !     destcomp = destination component index
    !     The fldsSrc input argument is queried for the mapping type of the field
    !     for the desination component
    !        mapindex = fldsSrc(n)%mapindex(destcomp)
    !     If the mapindex is 0 (there is no valid mapping) then NO mapping is done
    !        for the field
    !---------------------------------------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_Field
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleCreate
    use ESMF                  , only : ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldDestroy
    use ESMF                  , only : ESMF_Mesh, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT
    use med_methods_mod       , only : med_methods_FB_getFieldN, med_methods_FB_getNameN
    use med_constants_mod     , only : czero => med_constants_czero
    use esmFlds               , only : med_fldList_GetfldListFr, med_fldlist_type
    use esmFlds               , only : med_fld_GetFldInfo, med_fldList_entry_type
    use med_internalstate_mod , only : mapunset, compname, compocn, compatm
    use med_internalstate_mod , only : ncomps, nmappers, compname, mapnames, mapfcopy

    ! input/output variables
    type(ESMF_GridComp)          :: gcomp
    character(len=*), intent(in) :: flds_scalar_name
    integer, intent(in)          :: llogunit
    integer, intent(out)         :: rc

    ! local variables
    type(InternalState)       :: is_local
    type(ESMF_Field)          :: fldsrc
    type(ESMF_Field)          :: flddst
    integer                   :: n1,n2
    integer                   :: nf
    integer                   :: fieldCount
    type(ESMF_Field), pointer :: fieldlist(:)
    type(ESMF_Field)          :: field_src
    character(len=CX)         :: mapfile
    integer                   :: mapindex
    logical                   :: mapexists = .false.
    real(R8), pointer         :: dataptr(:)
    type(ESMF_Mesh)           :: mesh_src
    type(ESMF_Mesh)           :: mesh_dst
    type(med_fldlist_type), pointer :: FldListFr
    type(med_fldlist_entry_type), pointer :: fldptr
    character(len=*), parameter :: subname=' (module_med_map: RouteHandles_init) '
    !-----------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": start", ESMF_LOGMSG_INFO)
    endif

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! --------------------------------------------------------------
    ! Create the necessary route handles
    ! --------------------------------------------------------------

    ! First loop over source and destination components components
    if (mastertask) write(logunit,*) ' '
    do n1 = 1, ncomps
       do n2 = 1, ncomps
          if (n1 /= n2) then
             if (is_local%wrap%med_coupling_active(n1,n2)) then ! If coupling is active between n1 and n2
                ! Get source field
                call med_methods_FB_getFieldN(is_local%wrap%FBImp(n1,n1), 1, fldsrc, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Check number of fields in source FB on destination mesh and get destination field
                if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(n1,n2))) then
                   call ESMF_LogWrite(trim(subname)//'FBImp('//trim(compname(n1))//','//trim(compname(n2))//')'// &
                        ' has not been created', ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
                   rc = ESMF_FAILURE
                   return
                end if
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n2), fieldCount=fieldCount, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (fieldCount == 0) then
                  call med_methods_FB_getFieldN(is_local%wrap%FBExp(n2), 1, flddst, rc)
                  if (chkerr(rc,__LINE__,u_FILE_u)) return
                else
                  call med_methods_FB_getFieldN(is_local%wrap%FBImp(n1,n2), 1, flddst, rc)
                  if (chkerr(rc,__LINE__,u_FILE_u)) return
                end if

                ! Loop over fields
                fldListFr => med_fldList_getFldListFr(n1)
                fldptr => fldListFr%fields
                nf = 0
                do while(associated(fldptr))
                   ! Determine the mapping type for mapping field nf from n1 to n2
                   call med_fld_GetFldInfo(fldptr, compsrc=n2, mapindex=mapindex)
                   if (mapindex /= mapunset) then

                      ! determine if route handle has already been created
                      mapexists = med_map_RH_is_created(is_local%wrap%RH,n1,n2,mapindex,rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return

                      ! Create route handle for target mapindex if route handle is required
                      ! (i.e. mapindex /= mapunset) and route handle has not already been created
                      if (.not. mapexists) then
                         call med_fld_GetFldInfo(fldptr, compsrc=n2, mapfile=mapfile)
                         call med_map_routehandles_initfrom_field(n1, n2, fldsrc, flddst, &
                              mapindex, is_local%wrap%rh(n1,n2,:), mapfile=trim(mapfile), rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                      end if

                   end if ! end if mapindex is mapunset
                   fldptr => fldptr%next
                end do ! loop over fields

                
             end if ! if coupling active
          end if ! if n1 not equal to n2
       end do ! loop over n2
    end do ! loop over n1

    ! --------------------------------------------------------------
    ! Initialize unity normalization fields and do the mapping for
    ! unity normalization up front
    ! --------------------------------------------------------------

    if (mastertask) then
       write(logunit,*)
       write(logunit,'(a)') trim(subname)//"Initializing unity map normalizations"
    endif

    ! Create the destination normalization field
    do n1 = 1,ncomps

       ! Since coupling could be uni-directional, the import FB could be
       ! available but number of fields could be zero, so it is better to
       ! check export FB if this is the case
       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(n1,n1)) .or. &
           ESMF_FieldBundleIsCreated(is_local%wrap%FBExp(n1))) then

          ! Get source mesh
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (fieldCount == 0) then
            if (mastertask) then
              write(logunit,*) trim(subname)//' '//trim(compname(n1))//' import FB field count is = ', fieldCount
              write(logunit,*) trim(subname)//' '//trim(compname(n1))//' trying to use export FB'
            end if
            call ESMF_FieldBundleGet(is_local%wrap%FBExp(n1), fieldCount=fieldCount, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            allocate(fieldlist(fieldcount))
            call ESMF_FieldBundleGet(is_local%wrap%FBExp(n1), fieldlist=fieldlist, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
            allocate(fieldlist(fieldcount))
            call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldlist=fieldlist, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_FieldGet(fieldlist(1), mesh=mesh_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_src = ESMF_FieldCreate(mesh_src, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_src, farrayptr=dataPtr, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr(:) = 1.0_R8

          ! Loop over destination components
          do n2 = 1,ncomps
             if ( n1 /= n2 .and. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(n1,n2)) .and. &
                  is_local%wrap%med_coupling_active(n1,n2)) then

                ! Get destination mesh
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n2), fieldlist=fieldlist, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(fieldlist(1), mesh=mesh_dst, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Create is_local%wrap%field_NormOne(n1,n2,mapindex) if appropriate (don't create if mapping is redist)
                do mapindex = 1,nmappers
                   if (mapindex /= mapfcopy .and. med_map_RH_is_created(is_local%wrap%RH,n1,n2,mapindex,rc=rc)) then
                      is_local%wrap%field_NormOne(n1,n2,mapindex) = ESMF_FieldCreate(mesh_dst, &
                           ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      call ESMF_FieldGet(is_local%wrap%field_NormOne(n1,n2,mapindex), farrayptr=dataptr, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      dataptr(:) = czero
                      call med_map_field(field_src=field_src, field_dst=is_local%wrap%field_NormOne(n1,n2,mapindex), &
                           routehandles=is_local%wrap%RH(n1,n2,:), maptype=mapindex, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      if (mastertask) then
                         write(logunit,'(a)') trim(subname)//' created field_NormOne for '&
                              //compname(n1)//'->'//compname(n2)//' with mapping '//trim(mapnames(mapindex))
                      end if
                   end if
                end do ! end of loop over map_indiex mappers
             end if ! end of if block for creating destination field
          end do ! end of loop over n2

          ! Deallocate memory
          deallocate(fieldlist)
          call ESMF_FieldDestroy(field_src, rc=rc, noGarbage=.true.)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

       end if ! end of if-block for existence of field bundle
    end do ! end of loop over n1

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_RouteHandles_initfrom_esmflds

  !================================================================================
  subroutine med_map_routehandles_initfrom_fieldbundle(n1, n2, FBsrc, FBdst, mapindex, RouteHandle, rc)

    use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF            , only : ESMf_Field, ESMF_FieldBundle, ESMF_RouteHandle
    use med_methods_mod , only : med_methods_FB_getFieldN

    !---------------------------------------------
    ! Initialize initialize additional route handles for mapping fractions
    !---------------------------------------------

    ! input/output variables
    integer                , intent(in)    :: n1
    integer                , intent(in)    :: n2
    type(ESMF_FieldBundle) , intent(in)    :: FBsrc
    type(ESMF_FieldBundle) , intent(in)    :: fBdst
    integer                , intent(in)    :: mapindex
    type(ESMF_RouteHandle) , intent(inout) :: RouteHandle(:,:,:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)   :: fldsrc
    type(ESMF_Field)   :: flddst
    character(len=*), parameter :: subname=' (module_MED_map:med_map_routehandles_initfrom_fieldbundle) '
    !---------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": start", ESMF_LOGMSG_INFO)
    endif

    call med_methods_FB_getFieldN(FBsrc, 1, fldsrc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call med_methods_FB_getFieldN(FBDst, 1, flddst, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call med_map_routehandles_initfrom_field(n1, n2, fldsrc, flddst, mapindex, routehandle(n1,n2,:), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_routehandles_initfrom_fieldbundle

  !================================================================================
  subroutine med_map_routehandles_initfrom_field(n1, n2, fldsrc, flddst, mapindex, routehandles, mapfile, rc)

    use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandlePrint, ESMF_Field, ESMF_MAXSTR
    use ESMF                  , only : ESMF_PoleMethod_Flag, ESMF_POLEMETHOD_ALLAVG, ESMF_POLEMETHOD_NONE
    use ESMF                  , only : ESMF_FieldSMMStore, ESMF_FieldRedistStore, ESMF_FieldRegridStore
    use ESMF                  , only : ESMF_RouteHandleIsCreated, ESMF_RouteHandleCreate
    use ESMF                  , only : ESMF_REGRIDMETHOD_BILINEAR, ESMF_REGRIDMETHOD_PATCH
    use ESMF                  , only : ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_DSTAREA, ESMF_NORMTYPE_FRACAREA
    use ESMF                  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_NEAREST_STOD
    use ESMF                  , only : ESMF_EXTRAPMETHOD_NEAREST_STOD
    use ESMF                  , only : ESMF_Mesh, ESMF_MeshLoc, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_I4
    use ESMF                  , only : ESMF_MeshGet, ESMF_DistGridGet, ESMF_DistGrid, ESMF_TYPEKIND_R8
    use ESMF                  , only : ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldWrite, ESMF_FieldDestroy
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mappatch_uv3d, mapbilnr_uv3d, mapfcopy
    use med_internalstate_mod , only : mapunset, mapnames, nmappers
    use med_internalstate_mod , only : mapnstod, mapnstod_consd, mapnstod_consf, mapnstod_consd
    use med_internalstate_mod , only : mapfillv_bilnr, mapbilnr_nstod, mapconsf_aofrac
    use med_internalstate_mod , only : compocn, compwav, complnd, compname
    use med_internalstate_mod , only : coupling_mode, dststatus_print
    use med_internalstate_mod , only : defaultMasks
    use med_constants_mod     , only : ispval_mask => med_constants_ispval_mask

    ! input/output variables
    integer                    , intent(in)    :: n1
    integer                    , intent(in)    :: n2
    type(ESMF_Field)           , intent(inout) :: fldsrc
    type(ESMF_Field)           , intent(inout) :: flddst
    integer                    , intent(in)    :: mapindex
    type(ESMF_RouteHandle)     , intent(inout) :: routehandles(:)
    character(len=*), optional , intent(in)    :: mapfile
    integer                    , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)            :: dstmesh
    type(ESMF_Field)           :: dststatusfield, doffield
    type(ESMF_DistGrid)        :: distgrid
    character(len=CS)          :: string
    character(len=CS)          :: mapname
    character(len=CL)          :: fname
    integer                    :: srcMaskValue
    integer                    :: dstMaskValue
    character(len=ESMF_MAXSTR) :: lmapfile
    logical                    :: rhprint = .false., ldstprint = .false.
    integer                    :: ns
    integer(I4), pointer       :: dof(:)
    integer                    :: srcTermProcessing_Value = 0
    type(ESMF_PoleMethod_Flag) :: polemethod
    character(len=*), parameter :: subname=' (module_med_map: med_map_routehandles_initfrom_field) '
    !---------------------------------------------

    lmapfile = 'unset'
    if (present(mapfile)) then
       lmapfile = trim(mapfile)
    end if

    mapname = trim(mapnames(mapindex))
    call ESMF_LogWrite(trim(subname)//": mapname "//trim(mapname), ESMF_LOGMSG_INFO)

    ! create a field to retrieve the dststatus field
    call ESMF_FieldGet(flddst, mesh=dstmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dststatusfield = ESMF_FieldCreate(dstmesh, ESMF_TYPEKIND_I4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! set local flag to false
    ldstprint = .false.

    ! set src and dst masking using defaults
    srcMaskValue = defaultMasks(n1,1)
    dstMaskValue = defaultMasks(n2,2)

    ! override defaults for specific cases
    if (trim(coupling_mode) == 'cesm') then
       if (n1 == compwav .and. n2 == compocn) then
         srcMaskValue = 0
         dstMaskValue = ispval_mask
      endif
    end if
    if (trim(coupling_mode(1:4)) == 'nems') then
       if (n1 == compatm .and. n2 == complnd) then
          srcMaskValue = ispval_mask
          dstMaskValue = ispval_mask
       end if
    end if
    if (trim(coupling_mode) == 'hafs') then
       if (n1 == compatm .and. n2 == compwav) then
          srcMaskValue = ispval_mask
       end if
    end if
    write(string,'(a,i10,a,i10)') trim(compname(n1))//' to '//trim(compname(n2))//' srcMask = ', &
               srcMaskValue,' dstMask = ',dstMaskValue
    call ESMF_LogWrite(trim(string), ESMF_LOGMSG_INFO)

    polemethod=ESMF_POLEMETHOD_ALLAVG
    if (trim(coupling_mode) == 'cesm' .or. trim(coupling_mode(1:4)) == 'nems') then
       if (n1 == compwav .or. n2 == compwav) then
         polemethod = ESMF_POLEMETHOD_NONE ! todo: remove this when ESMF tripolar mapping fix is in place.
       endif
    end if

    ! Create route handle
    if (mapindex == mapfcopy) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH redist for '//trim(string)
       end if
       call ESMF_FieldRedistStore(fldsrc, flddst, routehandle=routehandles(mapfcopy), &
            ignoreUnmatchedIndices = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (lmapfile /= 'unset') then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//&
               ' via input file '//trim(mapfile)//' for '//trim(string)
       end if
       call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, routehandle=routehandles(mapindex), &
            ignoreUnmatchedIndices=.true., &
            srcTermProcessing=srcTermProcessing_Value, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mapbilnr .or. mapindex == mapbilnr_uv3d) then
       if (.not. ESMF_RouteHandleIsCreated(routehandles(mapbilnr))) then
          if (mastertask) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          end if
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapbilnr), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
               polemethod=polemethod, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ldstprint = .true.
       end if
    else if (mapindex == mapfillv_bilnr) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       end if
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapfillv_bilnr), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapbilnr_nstod) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       end if
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapbilnr_nstod), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapconsf .or. mapindex == mapnstod_consf) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       end if
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsf), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_FRACAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapconsf_aofrac) then
       if (.not. ESMF_RouteHandleIsCreated(routehandles(mapconsf))) then
          if (mastertask) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          end if
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsf_aofrac), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
               normType=ESMF_NORMTYPE_FRACAREA, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ldstprint = .true.
       else
          ! Copy existing consf RH
          if (mastertask) then
             write(logunit,'(A)') trim(subname)//' copying RH(mapconsf) to '//trim(mapname)//' for '//trim(string)
          end if
          routehandles(mapconsf_aofrac) = ESMF_RouteHandleCreate(routehandles(mapconsf), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    else if (mapindex == mapconsd .or. mapindex == mapnstod_consd) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       end if
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsd), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_DSTAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mappatch .or. mapindex == mappatch_uv3d) then
       if (.not. ESMF_RouteHandleIsCreated(routehandles(mappatch))) then
          if (mastertask) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          end if
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mappatch), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_PATCH, &
               polemethod=polemethod, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ldstprint = .true.
       end if
    else
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' mapindex '//trim(mapname)//' not supported for '//trim(string)
       end if
       call ESMF_LogWrite(trim(subname)//' mapindex '//trim(mapname)//' not supported ', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

    ! Output destination status field to file if requested
    if (dststatus_print .and. ldstprint) then
      fname = 'dststatus.'//trim(compname(n1))//'.'//trim(compname(n2))//'.'//trim(mapname)//'.nc'
      call ESMF_LogWrite(trim(subname)//": writing dstStatusField to "//trim(fname), ESMF_LOGMSG_INFO)

      call ESMF_FieldWrite(dststatusfield, filename=trim(fname), variableName='dststatus', &
           overwrite=.true., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! the sequence index in order to sort the dststatus field
      call ESMF_MeshGet(dstmesh, elementDistgrid=distgrid, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(dof(ns))
      call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      doffield = ESMF_FieldCreate(dstmesh, dof, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldWrite(doffield, fileName='dof.'//trim(compname(n2))//'.nc', variableName='dof', &
           overwrite=.true., rc=rc)
      deallocate(dof)
      call ESMF_FieldDestroy(doffield, rc=rc, noGarbage=.true.)
    end if

    ! consd_nstod method requires a second routehandle
    if (mapindex == mapnstod .or. mapindex == mapnstod_consd .or. mapindex == mapnstod_consf) then
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapnstod), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.

       ! Output destination status field to file if requested
       if (dststatus_print .and. ldstprint) then
          fname = 'dststatus.'//trim(compname(n1))//'.'//trim(compname(n2))//'.'//trim(mapname)//'_2.nc'
          call ESMF_LogWrite(trim(subname)//": writing dstStatusField to "//trim(fname), ESMF_LOGMSG_INFO)

          call ESMF_FieldWrite(dststatusfield, filename=trim(fname), variableName='dststatus', overwrite=.true., rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Output route handle to file if requested
    if (rhprint) then
       if (mastertask) then
          write(logunit,'(a)') trim(subname)//trim(string)//": printing  RH for "//trim(mapname)
       end if
       call ESMF_RouteHandlePrint(routehandles(mapindex), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_FieldDestroy(dststatusfield, rc=rc, noGarbage=.true.)

  end subroutine med_map_routehandles_initfrom_field

  !================================================================================
  logical function med_map_RH_is_created_RH3d(RHs,n1,n2,mapindex,rc)

    use ESMF, only : ESMF_RouteHandle

    ! input/output variables
    type(ESMF_RouteHandle) , intent(in)    :: RHs(:,:,:)
    integer                , intent(in)    :: n1
    integer                , intent(in)    :: n2
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=' (module_MED_map:med_map_RH_is_created_RH3d) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS
    med_map_RH_is_created_RH3d = med_map_RH_is_created_RH1d(RHs(n1,n2,:),mapindex,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end function med_map_RH_is_created_RH3d

!================================================================================

  logical function med_map_RH_is_created_RH1d(RHs,mapindex,rc)

    use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated
    use med_internalstate_mod , only : mapconsd, mapconsf, mapnstod
    use med_internalstate_mod , only : mapnstod_consd, mapnstod_consf

    ! input/output varaibes
    type(ESMF_RouteHandle) , intent(in)    :: RHs(:)
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    integer :: rc1, rc2
    logical :: mapexists
    character(len=*), parameter :: subname=' (module_MED_map:med_map_RH_is_created_RH1d) '
    !-----------------------------------------------------------

    rc  = ESMF_SUCCESS
    rc1 = ESMF_SUCCESS
    rc2 = ESMF_SUCCESS

    mapexists = .false.
    if      (mapindex == mapnstod_consd .and. &
             ESMF_RouteHandleIsCreated(RHs(mapnstod), rc=rc1) .and. &
             ESMF_RouteHandleIsCreated(RHs(mapconsd), rc=rc2)) then
       rc = rc1
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rc = rc2
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       mapexists = .true.
    else if (mapindex == mapnstod_consf .and. &
             ESMF_RouteHandleIsCreated(RHs(mapnstod), rc=rc1) .and. &
             ESMF_RouteHandleIsCreated(RHs(mapconsf), rc=rc2)) then
       rc = rc1
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rc = rc2
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       mapexists = .true.
    else if (ESMF_RouteHandleIsCreated(RHs(mapindex), rc=rc1)) then
       rc = rc1
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       mapexists = .true.
    end if
    med_map_RH_is_created_RH1d = mapexists

  end function med_map_RH_is_created_RH1d

  !================================================================================
  subroutine med_map_packed_field_create(destcomp, flds_scalar_name, &
       fieldsSrc, FBSrc, FBDst, packed_data, rc)

    use ESMF
    use esmFlds               , only : med_fldList_entry_type, med_fldList_getNumFlds, med_fldList_type
    use esmFlds               , only : med_fld_getFldInfo
    use med_internalstate_mod , only : ncomps, compname, mapnames
    use med_internalstate_mod , only : packed_data_type

    ! input/output variables
    integer                      , intent(in)    :: destcomp
    character(len=*)             , intent(in)    :: flds_scalar_name
    type(med_fldList_type)       , intent(in), target    :: fieldsSrc  ! mapping types top of LL
    type(ESMF_FieldBundle)       , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)       , intent(inout) :: FBDst
    type(packed_data_type)       , intent(inout) :: packed_data(:) ! array over mapping types
    integer                      , intent(out)   :: rc

    ! local variables
    integer                    :: nf, nu
    integer, allocatable       :: npacked(:)
    integer                    :: fieldcount
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
    type(med_fldlist_entry_type), pointer :: fldptr
    character(CL)              :: shortname
    integer                    :: destindex
    character(CL), allocatable :: fieldNameList(:)
    character(CS)              :: mapnorm_mapindex
    character(len=CX)          :: tmpstr
    character(len=*), parameter :: subname=' (module_MED_map:med_packed_field_create) '
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
    call ESMF_FieldGet(fieldlist_src(1), mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, numOwnedElements=lsize_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine local size of destination fields
    ! Allocate a destination fortran pointer for the new packed field bundle
    call ESMF_FieldGet(fieldlist_dst(1), mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, numOwnedElements=lsize_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Gather all fields that will be mapped with a target map index into a packed field
    ! Calculated size of packed field based on the fact that some fields have
    ! ungridded dimensions and need to unwrap them into separate fields for the
    ! purposes of packing

    if (mastertask) write(logunit,*)

    ! Determine the normalization type for each packed_data mapping element
    ! Loop over mapping types
    do mapindex = 1,nmappers
       mapnorm_mapindex = 'not_set'
       ! Loop over source field bundle
       do nf = 1, fieldCount
          ! Loop over the fldsSrc types
          fldptr => fieldsSrc%fields
          do while(associated(fldptr))
             ! Note that fieldnamelist is an array of names for the source fields
             ! The assumption is that there is only one mapping normalization
             ! for any given mapping type
             call med_fld_GetFldInfo(fldptr, compsrc=destcomp, shortname=shortname, mapindex=destindex)
             if ( destindex == mapindex .and. &
                  trim(shortname) == trim(fieldnamelist(nf))) then
                ! Set the normalization to the input
                call med_Fld_GetFldInfo(fldptr, compsrc=destcomp, mapnorm=packed_data(mapindex)%mapnorm)
                if (mapnorm_mapindex == 'not_set') then
                   mapnorm_mapindex = packed_data(mapindex)%mapnorm
                   write(tmpstr,*)'Map type '//trim(mapnames(mapindex)) &
                      //', destcomp '//trim(compname(destcomp)) &
                      //',  mapnorm '//trim(mapnorm_mapindex) &
                      //'  '//trim(fieldnamelist(nf))
                   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
                else
                   if (mapnorm_mapindex /= packed_data(mapindex)%mapnorm) then
                     write(tmpstr,*)'Map type '//trim(mapnames(mapindex)) &
                        //', destcomp '//trim(compname(destcomp)) &
                        //',  mapnorm '//trim(mapnorm_mapindex) &
                        //' set; cannot set mapnorm to '//trim(packed_data(mapindex)%mapnorm) &
                        //'  '//trim(fieldnamelist(nf))
                     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_ERROR)
                     call ESMF_Finalize(endflag=ESMF_END_ABORT)
                   end if
                end if
             end if
             fldptr => fldptr%next
          end do
       end do
    end do

    ! Allocate memory to keep tracked of packing index for each mapping type
    allocate(npacked(nmappers))
    npacked(:) = 0

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
          fldptr => fieldsSrc%fields
          do while(associated(fldptr))
             call med_fld_GetFldInfo(fldptr, compsrc=destcomp, shortname=shortname, mapindex=destIndex)
             if ( destIndex == mapindex .and. &
                  trim(shortname) == trim(fieldnamelist(nf))) then

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
                        trim(packed_data(mapindex)%mapnorm), &
                        trim(fieldnamelist(nf)), &
                        packed_data(mapindex)%fldindex(nf)
                end if

             end if! end if source field is mapped to destination field with mapindex
             fldptr => fldptr%next
          end do ! end loop over FBSrc fields
       end do ! end loop over fldsSrc elements

       if (npacked(mapindex) > 0) then
          ! Create the packed source field bundle for mapindex
          allocate(ptrsrc_packed(npacked(mapindex), lsize_src))
          packed_data(mapindex)%field_src = ESMF_FieldCreate(lmesh_src, &
               ptrsrc_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Create the packed destination field bundle for mapindex
          allocate(ptrdst_packed(npacked(mapindex), lsize_dst))
          packed_data(mapindex)%field_dst = ESMF_FieldCreate(lmesh_dst, &
               ptrdst_packed, gridToFieldMap=(/2/),  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          packed_data(mapindex)%field_fracsrc = ESMF_FieldCreate(lmesh_src, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          packed_data(mapindex)%field_fracdst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end do ! end loop over mapindex

    deallocate(npacked)
    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

  end subroutine med_map_packed_field_create

  !================================================================================
  subroutine med_map_field_packed(FBSrc, FBDst, FBFracSrc, field_normOne, packed_data, routehandles, rc)

    ! -----------------------------------------------
    ! Do regridding via packed field bundles
    ! -----------------------------------------------

    use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldIsCreated
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_FieldRedist, ESMF_RouteHandle
    use med_internalstate_mod , only : nmappers, mapfcopy
    use med_internalstate_mod , only : mappatch_uv3d, mappatch, mapbilnr_uv3d, mapbilnr
    use med_internalstate_mod , only : packed_data_type

    ! input/output variables
    type(ESMF_FieldBundle)    , intent(in)    :: FBSrc
    type(ESMF_FieldBundle)    , intent(inout) :: FBDst
    type(ESMF_Field)          , intent(in)    :: field_normOne(:)  ! array over mapping types
    type(ESMF_FieldBundle)    , intent(in)    :: FBFracSrc         ! fraction field bundle for source
    type(packed_data_type)    , intent(inout) :: packed_data(:)    ! array over mapping types
    type(ESMF_RouteHandle)    , intent(inout) :: routehandles(:)
    integer                   , intent(out)   :: rc

    ! local variables
    integer                    :: nf, nu, np, n
    integer                    :: fieldcount
    integer                    :: mapindex
    integer                    :: ungriddedUBound(1)
    real(r8), pointer          :: dataptr1d(:)
    real(r8), pointer          :: dataptr2d(:,:)
    real(r8), pointer          :: dataptr2d_packed(:,:)
    type(ESMF_Field)           :: field_fracsrc
    type(ESMF_Field), pointer  :: fieldlist_src(:)
    type(ESMF_Field), pointer  :: fieldlist_dst(:)
    real(r8), pointer          :: data_norm(:)
    real(r8), pointer          :: data_dst(:,:)
    character(len=*), parameter  :: subname=' (module_MED_map:med_map_field_packed) '
    !-----------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    ! Get field count for both FBsrc and FBdst
    if (ESMF_FieldBundleIsCreated(FBsrc)) then
       call ESMF_FieldBundleGet(FBsrc, fieldCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(fieldlist_src(fieldcount))
       call ESMF_FieldBundleGet(FBsrc, fieldlist=fieldlist_src, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(fieldlist_dst(fieldcount))
       call ESMF_FieldBundleGet(FBdst, fieldlist=fieldlist_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       fieldcount=0
    endif


    ! Loop over mapping types
    do mapindex = 1,nmappers

       if (ESMF_FieldIsCreated(packed_data(mapindex)%field_src)) then

          if (mapindex == mappatch_uv3d) then

             ! For mappatch_uv3d do not use packed field bundles
             call med_map_uv_cart3d(FBsrc, FBdst, routehandles, mappatch, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          else if (mapindex == mapbilnr_uv3d) then

             ! For mapbilnr_uv3d do not use packed field bundles
             call med_map_uv_cart3d(FBsrc, FBdst, routehandles, mapbilnr, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

          else

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
                call ESMF_FieldRedist(&
                     packed_data(mapindex)%field_src, &
                     packed_data(mapindex)%field_dst, &
                     routehandles(mapindex), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

             else if ( trim(packed_data(mapindex)%mapnorm) /= 'unset' .and. &
                  trim(packed_data(mapindex)%mapnorm) /= 'one'   .and. &
                  trim(packed_data(mapindex)%mapnorm) /= 'none') then

                ! Normalized mapping - assume that  each packed field has only one normalization type
                call ESMF_FieldBundleGet(FBFracSrc, packed_data(mapindex)%mapnorm, field=field_fracsrc, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call med_map_field_normalized(&
                     field_src=packed_data(mapindex)%field_src, &
                     field_dst=packed_data(mapindex)%field_dst, &
                     routehandles=routehandles, &
                     maptype=mapindex, &
                     field_normsrc=field_fracsrc, &
                     field_normdst=packed_data(mapindex)%field_fracdst, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

             else if ( trim(packed_data(mapindex)%mapnorm) == 'one' .or. trim(packed_data(mapindex)%mapnorm) == 'none') then

                ! Mapping with no normalization that is not redistribution
                call med_map_field (&
                     field_src=packed_data(mapindex)%field_src, &
                     field_dst=packed_data(mapindex)%field_dst, &
                     routehandles=routehandles, &
                     maptype=mapindex, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Obtain unity normalization factor and multiply
                ! interpolated field by reciprocal of normalization factor
                if (trim(packed_data(mapindex)%mapnorm) == 'one') then
                   call ESMF_FieldGet(field_normOne(mapindex), farrayPtr=data_norm, rc=rc)
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
       end if
    end do ! end of loop over mapindex

    if (ESMF_FieldBundleIsCreated(FBsrc)) then
      deallocate(fieldlist_src)
      deallocate(fieldlist_dst)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_map_field_packed

  !================================================================================
  subroutine med_map_field_normalized(field_src, field_dst, routehandles, maptype, &
       field_normsrc, field_normdst, rc)

    ! -----------------------------------------------
    ! Map a normalized field
    ! -----------------------------------------------

    use ESMF        , only : ESMF_Field, ESMF_FieldGet, ESMF_RouteHandle
    use ESMF        , only : ESMF_SUCCESS

    ! input/output variables
    type(ESMF_Field)       , intent(in)    :: field_src
    type(ESMF_Field)       , intent(inout) :: field_dst
    type(ESMF_Field)       , intent(in)    :: field_normsrc
    type(ESMF_Field)       , intent(inout) :: field_normdst
    type(ESMF_RouteHandle) , intent(inout) :: routehandles(:)
    integer                , intent(in)    :: maptype
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n
    real(r8), pointer :: data_src2d(:,:)
    real(r8), pointer :: data_dst2d(:,:)
    real(r8), pointer :: data_srctmp2d(:,:)
    real(r8), pointer :: data_src1d(:)
    real(r8), pointer :: data_dst1d(:)
    real(r8), pointer :: data_srctmp1d(:)
    real(r8), pointer :: data_normsrc(:)
    real(r8), pointer :: data_normdst(:)
    integer           :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    integer           :: lsize_src
    integer           :: lsize_dst
    character(len=*), parameter  :: subname=' (module_MED_map:med_map_field_normalized) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! get a pointer (data_fracsrc) to the normalization array
    ! get a pointer (data_src) to source field data in FBSrc
    ! copy data_src to data_srctmp

    ! normalize data_src by data_fracsrc

    call ESMF_FieldGet(field_normsrc, farrayPtr=data_normsrc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize_src = size(data_normsrc)

    call ESMF_FieldGet(field_src, ungriddedUBound=ungriddedUBound, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (ungriddedUbound(1) > 0) then
       call ESMF_FieldGet(field_src, farrayPtr=data_src2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(data_srctmp2d(size(data_src2d,dim=1), lsize_src))
       data_srctmp2d(:,:) = data_src2d(:,:)
       do n = 1,lsize_src
          data_src2d(:,n) = data_src2d(:,n) * data_normsrc(n)
       end do
    else
       call ESMF_FieldGet(field_src, farrayPtr=data_src1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(data_srctmp1d(lsize_src))
       data_srctmp1d(:) = data_src1d(:)
       do n = 1,lsize_src
          data_src1d(n) = data_src1d(n) * data_normsrc(n)
       end do
    end if

    ! regrid normalized packed source field
    call med_map_field (field_src=field_src, field_dst=field_dst, routehandles=routehandles, maptype=maptype, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! restore original value to packed source field
    if (ungriddedUbound(1) > 0) then
       data_src2d(:,:) = data_srctmp2d(:,:)
       deallocate(data_srctmp2d)
    else
       data_src1d(:) = data_srctmp1d(:)
       deallocate(data_srctmp1d)
    end if

    ! regrid normalization field from source to destination
    call med_map_field(field_src=field_normsrc, field_dst=field_normdst, routehandles=routehandles, maptype=maptype, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get pointer to mapped fraction and normalize
    ! destination mapped values by the reciprocal of the mapped fraction
    call ESMF_FieldGet(field_normdst, farrayPtr=data_normdst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize_dst = size(data_normdst)

    if (ungriddedUbound(1) > 0) then
       call ESMF_FieldGet(field_dst, farrayPtr=data_dst2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize_dst
          if (data_normdst(n) == 0.0_r8) then
             data_dst2d(:,n) = 0.0_r8
          else
             data_dst2d(:,n) = data_dst2d(:,n)/data_normdst(n)
          end if
       end do
    else
       call ESMF_FieldGet(field_dst, farrayPtr=data_dst1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize_dst
          if (data_normdst(n) == 0.0_r8) then
             data_dst1d(n) = 0.0_r8
          else
             data_dst1d(n) = data_dst1d(n)/data_normdst(n)
          end if
       end do
    end if
  end subroutine med_map_field_normalized

  !================================================================================
  subroutine med_map_field(field_src, field_dst, routehandles, maptype, fldname, rc)

    !---------------------------------------------------
    ! map the source field to the destination field
    !---------------------------------------------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_MAXSTR
    use ESMF                  , only : ESMF_KIND_R8
    use ESMF                  , only : ESMF_Field, ESMF_FieldRegrid
    use ESMF                  , only : ESMF_FieldFill
    use ESMF                  , only : ESMF_TERMORDER_SRCSEQ, ESMF_Region_Flag, ESMF_REGION_TOTAL
    use ESMF                  , only : ESMF_REGION_SELECT
    use ESMF                  , only : ESMF_RouteHandle
    use med_internalstate_mod , only : mapnstod_consd, mapnstod_consf, mapnstod_consd, mapnstod
    use med_internalstate_mod , only : mapconsd, mapconsf
    use med_internalstate_mod , only : mapfillv_bilnr
    use med_methods_mod       , only : Field_diagnose => med_methods_Field_diagnose

    ! input/output variables
    type(ESMF_Field)       , intent(in)           :: field_src
    type(ESMF_Field)       , intent(inout)        :: field_dst
    type(ESMF_RouteHandle) , intent(inout)        :: routehandles(:)
    integer                , intent(in)           :: maptype
    character(len=*)       , intent(in), optional :: fldname
    integer                , intent(out)          :: rc

    ! local variables
    logical :: checkflag = .false.
    character(len=CS) :: lfldname
    real(ESMF_KIND_R8), parameter :: fillValue = 9.99e20_ESMF_KIND_R8
    character(len=*), parameter :: subname='(module_MED_map:med_map_field) '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

#ifdef DEBUG
    checkflag = .true.
#endif
    lfldname = 'unknown'
    if (present(fldname)) lfldname = trim(fldname)

    if (maptype == mapnstod_consd) then
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(mapnstod), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after nstod: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(mapconsd), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_SELECT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after consd: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    else if (maptype == mapnstod_consf) then
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(mapnstod), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after nstod: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(mapconsf), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_SELECT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after consf: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    else if (maptype == mapfillv_bilnr) then
       call ESMF_FieldFill(field_dst, dataFillScheme="const", const1=fillValue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after fillv: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(mapfillv_bilnr), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_SELECT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call Field_diagnose(field_dst, lfldname, " --> after bilnr: ", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_map_field

  !================================================================================
  subroutine med_map_uv_cart3d(FBsrc, FBdst, routehandles, mapindex, rc)

    use ESMF          , only : ESMF_Mesh, ESMF_MeshGet, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
    use ESMF          , only : ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
    use ESMF          , only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF          , only : ESMF_RouteHandle
    use med_constants_mod , only : shr_const_pi

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FBsrc
    type(ESMF_FieldBundle) , intent(inout) :: FBdst
    type(ESMF_RouteHandle) , intent(inout) :: routehandles(:)
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)    :: usrc
    type(ESMF_Field)    :: vsrc
    type(ESMF_Field)    :: udst
    type(ESMF_Field)    :: vdst
    integer             :: n
    real(r8)            :: lon,lat
    real(r8)            :: coslon,coslat
    real(r8)            :: sinlon,sinlat
    real(r8)            :: ux,uy,uz
    type(ESMF_Mesh)     :: lmesh_src
    type(ESMF_Mesh)     :: lmesh_dst
    real(r8), pointer   :: data_u_src(:)
    real(r8), pointer   :: data_u_dst(:)
    real(r8), pointer   :: data_v_src(:)
    real(r8), pointer   :: data_v_dst(:)
    real(r8), pointer   :: data2d_src(:,:)
    real(r8), pointer   :: data2d_dst(:,:)
    real(r8), pointer   :: ownedElemCoords_src(:)
    real(r8), pointer   :: ownedElemCoords_dst(:)
    integer             :: numOwnedElements
    integer             :: spatialDim
    real(r8), parameter :: deg2rad = shr_const_pi/180.0_R8  ! deg to rads
    logical             :: first_time = .true.
    character(len=*), parameter :: subname=' (module_MED_map:med_map_uv_cart3d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get fields for atm u,v velocities
    call ESMF_FieldBundleGet(FBSrc, fieldName='Sa_u', field=usrc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBDst, fieldName='Sa_u', field=udst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBSrc, fieldName='Sa_v', field=vsrc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBDst, fieldName='Sa_v', field=vdst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! GET pointer to input u and v data source field data
    call ESMF_FieldGet(usrc, farrayPtr=data_u_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(vsrc, farrayPtr=data_v_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get pointer to destination data that will be filled in after
    ! rotation back from cart3d
    call ESMF_FieldGet(udst, farrayPtr=data_u_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(vdst, farrayPtr=data_v_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get source mesh and coordinates
    call ESMF_FieldGet(usrc, mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords_src(spatialDim*numOwnedElements))
    call ESMF_MeshGet(lmesh_src, ownedElemCoords=ownedElemCoords_src)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get destination mesh and coordinates
    call ESMF_FieldGet(udst, mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords_dst(spatialDim*numOwnedElements))
    call ESMF_MeshGet(lmesh_dst, ownedElemCoords=ownedElemCoords_dst)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (first_time) then
       ! Create two module fields -  vec_src3d and vec_dst3d -  that contain
       ! all three fields with undistributed dimensions for each
       uv3d_src = ESMF_FieldCreate(lmesh_src, ESMF_TYPEKIND_R8, name='src3d', &
            ungriddedLbound=(/1/), ungriddedUbound=(/3/), gridToFieldMap=(/2/), &
            meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       uv3d_dst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, name='dst3d', &
            ungriddedLbound=(/1/), ungriddedUbound=(/3/), gridToFieldMap=(/2/), &
            meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       first_time = .false.
    end if

    ! get pointers to source and destination data that will be filled in with rotation to cart3d
    call ESMF_FieldGet(uv3d_src, farrayPtr=data2d_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(uv3d_dst, farrayPtr=data2d_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Rotate Source data to cart3d
    do n = 1,size(data_u_src)
       lon = ownedElemCoords_src(2*n-1)
       lat = ownedElemCoords_src(2*n)
       sinlon = sin(lon*deg2rad)
       coslon = cos(lon*deg2rad)
       sinlat = sin(lat*deg2rad)
       coslat = cos(lat*deg2rad)
       data2d_src(1,n) = -coslon*sinlat*data_v_src(n) - sinlon*data_u_src(n) ! x
       data2d_src(2,n) = -sinlon*sinlat*data_v_src(n) + coslon*data_u_src(n) ! y
       data2d_src(3,n) =  coslat*data_v_src(n)                               ! z
    enddo

    ! Map all thee vector fields at once from source to destination grid
    call med_map_field(field_src=uv3d_src, field_dst=uv3d_dst, &
         routehandles=routehandles, maptype=mapindex, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Rotate destination data back from cart3d to original
    do n = 1,size(data_u_dst)
       lon = ownedElemCoords_dst(2*n-1)
       lat = ownedElemCoords_dst(2*n)
       sinlon = sin(lon*deg2rad)
       coslon = cos(lon*deg2rad)
       sinlat = sin(lat*deg2rad)
       coslat = cos(lat*deg2rad)
       ux = data2d_dst(1,n)
       uy = data2d_dst(2,n)
       uz = data2d_dst(3,n)
       data_u_dst(n) = -sinlon*ux + coslon*uy
       data_v_dst(n) = -coslon*sinlat*ux - sinlon*sinlat*uy + coslat*uz
    enddo

    ! Deallocate data
    deallocate(ownedElemCoords_src)
    deallocate(ownedElemCoords_dst)

  end subroutine med_map_uv_cart3d

end module med_map_mod
