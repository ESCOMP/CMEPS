module med_map_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use shr_const_mod         , only : shr_const_pi
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use esmFlds               , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy
  use esmFlds               , only : mapunset, mapnames, nmappers
  use esmFlds               , only : mapnstod, mapnstod_consd, mapnstod_consf, mapnstod_consd
  use esmFlds               , only : ncomps, compatm, compice, compocn, compname
  use esmFlds               , only : mapfcopy, mapconsd, mapconsf, mapnstod
  use esmFlds               , only : mapuv_with_cart3d, fldListFr, coupling_mode
  use med_internalstate_mod , only : InternalState, logunit, mastertask
  use med_constants_mod     , only : ispval_mask       => med_constants_ispval_mask
  use med_constants_mod     , only : czero             => med_constants_czero
  use med_constants_mod     , only : dbug_flag         => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr            => med_utils_ChkErr
  use med_utils_mod         , only : memcheck          => med_memcheck
  use med_methods_mod       , only : FB_getFieldN      => med_methods_FB_getFieldN
  use med_methods_mod       , only : FB_init           => med_methods_FB_Init
  use med_methods_mod       , only : FB_reset          => med_methods_FB_Reset
  use med_methods_mod       , only : FB_Clean          => med_methods_FB_Clean
  use med_methods_mod       , only : FB_Field_diagnose => med_methods_FB_Field_diagnose
  use med_methods_mod       , only : FB_FldChk         => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFieldByName => med_methods_FB_GetFieldByName
  use med_methods_mod       , only : Field_diagnose    => med_methods_Field_diagnose
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  ! public routines
  public :: med_map_routehandles_init
  public :: med_map_rh_is_created
  public :: med_map_mapnorm_init
  public :: med_map_packed_field_create
  public :: med_map_field_packed
  public :: med_map_field_normalized
  public :: med_map_field
  public :: med_map_fb_field_regrid  ! TODO: do we still need this?

  interface med_map_routehandles_init
     module procedure med_map_routehandles_init_esmflds
     module procedure med_map_routehandles_init_field
  end interface

  interface med_map_RH_is_created
     module procedure med_map_RH_is_created_RH3d
     module procedure med_map_RH_is_created_RH1d
  end interface

  ! private module variables

  character(len=CL)       :: flds_scalar_name
  integer                 :: srcTermProcessing_Value = 0 ! should this be a module variable?
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_map_RouteHandles_init_esmflds(gcomp, llogunit, rc)

    !---------------------------------------------
    ! Initialize route handles in the mediator
    ! Assumptions:
    !   -  Route handles are created per target field bundles NOT
    !      per individual fields in the bundle
    !   -  ALL fields in the bundle are on identical grids
    !   -  MULTIPLE route handles are going to be generated for
    !      given field bundle source and destination grids
    !    - Route handles will ONLY be created if coupling is active
    !      between n1 and n2
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

    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush, ESMF_KIND_I4
    use ESMF  , only : ESMF_GridComp, ESMF_VM, ESMF_Field, ESMF_PoleMethod_Flag, ESMF_POLEMETHOD_ALLAVG
    use ESMF  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_FieldSMMStore
    use ESMF  , only : ESMF_FieldRedistStore, ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_BILINEAR
    use ESMF  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_FRACAREA
    use ESMF  , only : ESMF_REGRIDMETHOD_NEAREST_STOD
    use ESMF  , only : ESMF_NORMTYPE_DSTAREA, ESMF_REGRIDMETHOD_PATCH, ESMF_RouteHandlePrint
    use NUOPC , only : NUOPC_Write

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: llogunit
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Field)        :: fldsrc
    type(ESMF_Field)        :: flddst
    integer                 :: localPet
    integer                 :: n,n1,n2,m,nf,nflds,ncomp
    integer                 :: SrcMaskValue
    integer                 :: DstMaskValue
    character(len=128)      :: value
    character(len=128)      :: rhname
    character(len=128)      :: rhname_file
    character(len=CS)       :: mapname
    character(len=CX)       :: mapfile
    character(len=CS)       :: string
    integer                 :: mapindex
    logical                 :: rhprint_flag = .false.
    logical                 :: mapexists = .false.
    real(R8)      , pointer :: factorList(:)
    character(CL) , pointer :: fldnames(:)
    !integer(ESMF_KIND_I4), pointer :: unmappedDstList(:)
    character(len=128)      :: logMsg
    type(ESMF_PoleMethod_Flag), parameter :: polemethod=ESMF_POLEMETHOD_ALLAVG
    character(len=*), parameter :: subname=' (RouteHandles_init) '
    !-----------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite("Starting to initialize RHs", ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    rc = ESMF_SUCCESS

    ! Determine mastertask
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    mastertask = .false.
    if (localPet == 0) mastertask=.true.
    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create the necessary route handles
    if (mastertask) write(llogunit,*) ' '
    do n1 = 1, ncomps
       do n2 = 1, ncomps

          if (trim(coupling_mode) == 'cesm') then
             dstMaskValue = ispval_mask
             srcMaskValue = ispval_mask
             if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
             if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
          else if (coupling_mode(1:4) == 'nems') then
             if (n1 == compatm .and. (n2 == compocn .or. n2 == compice)) then
                srcMaskValue = 1
                dstMaskValue = 0
             else if (n2 == compatm .and. (n1 == compocn .or. n1 == compice)) then
                srcMaskValue = 0
                dstMaskValue = 1
             else if ((n1 == compocn .and. n2 == compice) .or. (n1 == compice .and. n2 == compocn)) then
                srcMaskValue = 0
                dstMaskValue = 0
             else
                ! TODO: what should the condition be here?
                dstMaskValue = ispval_mask
                srcMaskValue = ispval_mask
             end if
          else if (trim(coupling_mode) == 'hafs') then
             dstMaskValue = ispval_mask
             srcMaskValue = ispval_mask
             if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
             if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
          end if

          !--- get single fields from bundles
          !--- 1) ASSUMES all fields in the bundle are on identical grids
          !--- 2) MULTIPLE route handles are going to be generated for
          !---    given field bundle source and destination grids

          if (n1 /= n2) then
             ! Determine route handle names
             rhname = trim(compname(n1))//"2"//trim(compname(n2))

             if (is_local%wrap%med_coupling_active(n1,n2)) then ! If coupling is active between n1 and n2

                call FB_GetFieldN(is_local%wrap%FBImp(n1,n1), 1, fldsrc, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                call FB_GetFieldN(is_local%wrap%FBImp(n1,n2), 1, flddst, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Loop over fields
                do nf = 1,size(fldListFr(n1)%flds)

                   ! Determine the mapping type for mapping field nf from n1 to n2
                   mapindex = fldListFr(n1)%flds(nf)%mapindex(n2)

                   ! separate check first since Fortran does not have short-circuit evaluation
                   if (mapindex == mapunset) cycle

                   ! Create route handle for target mapindex if route handle is required
                   ! (i.e. mapindex /= mapunset) and route handle has not already been created
                   mapexists = med_map_RH_is_created(is_local%wrap%RH,n1,n2,mapindex,rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return

                   if (.not. mapexists) then
                      mapname  = trim(mapnames(mapindex))
                      mapfile  = trim(fldListFr(n1)%flds(nf)%mapfile(n2))
                      string   = trim(rhname)//'_weights'

                      if (mapindex == mapfcopy) then
                         if (mastertask) then
                            write(llogunit,'(3A)') subname,trim(string),' RH redist '
                         end if
                         call ESMF_LogWrite(subname // trim(string) // ' RH redist ', ESMF_LOGMSG_INFO)
                         call ESMF_FieldRedistStore(fldsrc, flddst, &
                              routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                              ignoreUnmatchedIndices = .true., rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                      else if (mapfile /= 'unset') then
                         if (mastertask) then
                            write(llogunit,'(4A)') subname,trim(string),' RH '//trim(mapname)//' via input file ',&
                                 trim(mapfile)
                         end if
                         call ESMF_LogWrite(subname // trim(string) //&
                              ' RH '//trim(mapname)//' via input file '//trim(mapfile), ESMF_LOGMSG_INFO)
                         call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, &
                              routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                              ignoreUnmatchedIndices=.true., &
                              srcTermProcessing=srcTermProcessing_Value, rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                      else
                         ! Create route handle on the fly
                         if (mastertask) write(llogunit,'(3A)') subname,trim(string),&
                              ' RH regrid for '//trim(mapname)//' computed on the fly'
                         call ESMF_LogWrite(subname // trim(string) //&
                              ' RH regrid for '//trim(mapname)//' computed on the fly', ESMF_LOGMSG_INFO)
                         if (mapindex == mapbilnr) then
                            srcTermProcessing_Value = 0
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
                                 polemethod=polemethod, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         else if ((mapindex == mapconsf .or. mapindex == mapnstod_consf) .and. &
                              .not. med_map_RH_is_created(is_local%wrap%RH(n1,n2,:),mapconsf,rc)) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapconsf), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
                                 normType=ESMF_NORMTYPE_FRACAREA, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                                 rc=rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         else if ((mapindex == mapconsd .or. mapindex == mapnstod_consd) .and. &
                              .not. med_map_RH_is_created(is_local%wrap%RH(n1,n2,:),mapconsd,rc)) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapconsd), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
                                 normType=ESMF_NORMTYPE_DSTAREA, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                                 rc=rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         else if (mapindex == mappatch) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_PATCH, &
                                 polemethod=polemethod, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         end if
                         ! consd_nstod method requires a second routehandle
                         if ((mapindex == mapnstod .or. mapindex == mapnstod_consd .or. mapindex == mapnstod_consf) .and. &
                              .not. med_map_RH_is_created(is_local%wrap%RH(n1,n2,:),mapnstod,rc)) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapnstod), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                                 rc=rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         end if
                         if (rhprint_flag .and. mapindex /= mapnstod_consd .and. mapindex /= mapnstod_consf) then
                            call NUOPC_Write(factorList, "array_med_"//trim(string)//"_consf.nc", rc)
                            if (chkerr(rc,__LINE__,u_FILE_u)) return
                         end if
                         !if (associated(unmappedDstList)) then
                         !   write(logMsg,*) trim(subname),trim(string),&
                         !      number of unmapped dest points = ', size(unmappedDstList)
                         !   call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
                         !end if
                      end if
                      if (rhprint_flag .and. mapindex /= mapnstod_consd .and. mapindex /= mapnstod_consf) then
                         call ESMF_LogWrite(trim(subname)//trim(string)//": printing  RH for "//trim(mapname), &
                              ESMF_LOGMSG_INFO)
                         call ESMF_RouteHandlePrint(is_local%wrap%RH(n1,n2,mapindex), rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                      endif
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      ! Check that a valid route handle has been created
                      if (.not.med_map_RH_is_created(is_local%wrap%RH,n1,n2,mapindex,rc=rc)) then
                         call ESMF_LogWrite(trim(subname)//trim(string)//": failed   RH "//trim(mapname), &
                              ESMF_LOGMSG_INFO)
                      endif
                   end if
                end do ! loop over fields
             end if ! if coupling is active between n1 and n2
          end if ! if n1 not equal to n2
       end do ! loop over n2
    end do ! loop over n1

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_RouteHandles_init_esmflds

  !================================================================================
  subroutine med_map_routehandles_init_field(n1, n2, FBsrc, FBdst, mapindex, RouteHandle, rc)

    !---------------------------------------------
    ! Initialize initialize additional route handles
    ! for mapping fractions
    !
    !---------------------------------------------

    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF , only : ESMF_GridComp, ESMF_FieldBundle, ESMF_RouteHandle, ESMF_Field
    use ESMF , only : ESMF_GridComp, ESMF_VM, ESMF_Field, ESMF_PoleMethod_Flag, ESMF_POLEMETHOD_ALLAVG
    use ESMF , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_FieldSMMStore
    use ESMF , only : ESMF_FieldRedistStore, ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_BILINEAR
    use ESMF , only : ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_DSTAREA, ESMF_NORMTYPE_FRACAREA
    use ESMF , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_NEAREST_STOD
    use ESMF , only : ESMF_REGRIDMETHOD_PATCH

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
    character(len=CS)  :: string
    character(len=CS)  :: mapname
    integer            :: SrcMaskValue
    integer            :: DstMaskValue
    type(ESMF_PoleMethod_Flag), parameter :: polemethod=ESMF_POLEMETHOD_ALLAVG
    character(len=*), parameter :: subname=' (med_map_routehandles_init) '
    !---------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS
    if (dbug_flag > 1) then
       call ESMF_LogWrite("Initializing RHs not yet created and needed", &
            ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    mapname = trim(mapnames(mapindex))
    if (mastertask) then
       string = trim(compname(n1))//"2"//trim(compname(n2))//'_weights'
       write(logunit,'(3A)') subname, trim(string),&
            ' RH regrid for '//trim(mapname)//' computed on the fly'
    end if

    if (trim(coupling_mode) == 'cesm') then
       dstMaskValue = ispval_mask
       srcMaskValue = ispval_mask
       if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
       if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
    else if (coupling_mode(1:4) == 'nems') then
       if (n1 == compatm .and. (n2 == compocn .or. n2 == compice)) then
          srcMaskValue = 1
          dstMaskValue = 0
       else if (n2 == compatm .and. (n1 == compocn .or. n1 == compice)) then
          srcMaskValue = 0
          dstMaskValue = 1
       else if ((n1 == compocn .and. n2 == compice) .or. (n1 == compice .and. n2 == compocn)) then
          srcMaskValue = 0
          dstMaskValue = 0
       else
          ! TODO: what should the condition be here?
          dstMaskValue = ispval_mask
          srcMaskValue = ispval_mask
       end if
    else if (trim(coupling_mode) == 'hafs') then
       dstMaskValue = ispval_mask
       srcMaskValue = ispval_mask
       if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
       if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
    end if

    call FB_getFieldN(FBsrc, 1, fldsrc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFieldN(FBDst, 1, flddst, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create route handle on the fly
    if (mastertask) write(logunit,'(3A)') subname,trim(string),&
         ' RH regrid for '//trim(mapname)//' computed on the fly'
    call ESMF_LogWrite(subname // trim(string) //&
         ' RH regrid for '//trim(mapname)//' computed on the fly', ESMF_LOGMSG_INFO)

    if (mapindex == mapfcopy) then
       if (mastertask) then
          write(logunit,'(3A)') subname,trim(string),' RH redist '
       end if
       call ESMF_LogWrite(subname // trim(string) // ' RH redist ', ESMF_LOGMSG_INFO)
       call ESMF_FieldRedistStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapindex), &
            ignoreUnmatchedIndices = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mapbilnr) then
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapindex), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mapconsf .or. mapindex == mapnstod_consf) then
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapconsf), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_FRACAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mapconsd .or. mapindex == mapnstod_consd) then
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapconsd), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_DSTAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mappatch) then
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapindex), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_PATCH, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//' mapindex '//trim(mapnames(mapindex))//&
            ' not supported for fraction mapping', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

    ! consd_nstod method requires a second routehandle
    if (mapindex == mapnstod .or. mapindex == mapnstod_consd .or. mapindex == mapnstod_consf) then
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=routehandle(n1,n2,mapnstod), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_routehandles_init_field

  !================================================================================
  logical function med_map_RH_is_created_RH3d(RHs,n1,n2,mapindex,rc)

    use ESMF  , only : ESMF_RouteHandle

    ! input/output variables
    type(ESMF_RouteHandle) , intent(in)    :: RHs(:,:,:)
    integer                , intent(in)    :: n1
    integer                , intent(in)    :: n2
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    integer :: rc1, rc2
    logical :: mapexists
    character(len=*), parameter :: subname=' (med_map_RH_is_created) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS
    med_map_RH_is_created_RH3d = med_map_RH_is_created_RH1d(RHs(n1,n2,:),mapindex,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end function med_map_RH_is_created_RH3d

!================================================================================

  logical function med_map_RH_is_created_RH1d(RHs,mapindex,rc)

    use ESMF  , only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated

    ! input/output varaibes
    type(ESMF_RouteHandle) , intent(in)    :: RHs(:)
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    integer :: rc1, rc2
    logical :: mapexists
    character(len=*), parameter :: subname=' (med_map_RH_is_created_RH1d) '
    !-----------------------------------------------------------

    rc  = ESMF_SUCCESS
    rc1 = ESMF_SUCCESS
    rc2 = ESMF_SUCCESS

    mapexists = .false.
    if      (mapindex == mapnstod_consd .and. &
             ESMF_RouteHandleIsCreated(RHs(mapnstod), rc=rc1) .and. &
             ESMF_RouteHandleIsCreated(RHs(mapconsd), rc=rc2)) then
       mapexists = .true.
    else if (mapindex == mapnstod_consf .and. &
             ESMF_RouteHandleIsCreated(RHs(mapnstod), rc=rc1) .and. &
             ESMF_RouteHandleIsCreated(RHs(mapconsf), rc=rc2)) then
       mapexists = .true.
    else if (ESMF_RouteHandleIsCreated(RHs(mapindex), rc=rc1)) then
       mapexists = .true.
    end if

    med_map_RH_is_created_RH1d = mapexists

    rc = rc1
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    rc = rc2
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end function med_map_RH_is_created_RH1d

  !================================================================================
  subroutine med_map_mapnorm_init(gcomp, rc)

    !---------------------------------------
    ! Initialize unity normalization fields and do the mapping for unity normalization up front
    !---------------------------------------

    use ESMF , only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF , only: ESMF_GridComp
    use ESMF , only: ESMF_Mesh, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT
    use ESMF , only: ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleCreate
    use ESMF , only: ESMF_FieldBundleIsCreated
    use ESMF , only: ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldDestroy

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)       :: is_local
    integer                   :: n1, n2, m
    character(len=1)          :: cn1,cn2,cm
    real(R8), pointer         :: dataptr(:) => null()
    integer                   :: fieldCount
    type(ESMF_Field), pointer :: fieldlist(:) => null()
    type(ESMF_Field)          :: field_src
    type(ESMF_Mesh)           :: mesh_src
    type(ESMF_Mesh)           :: mesh_dst
    character(len=*),parameter :: subname=' (med_map_mapnorm_init) '
    !-----------------------------------------------------------
    call t_startf('MED:'//subname)

    rc = ESMF_SUCCESS

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": start", ESMF_LOGMSG_INFO)
    endif
    if (mastertask) then
       write(logunit,*)
       write(logunit,'(a)') trim(subname)//"Initializing unity map normalizations"
    endif

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create the destination normalization field
    do n1 = 1,ncomps

       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(n1,n1))) then
          ! Get source mesh
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldlist(fieldcount))
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldlist=fieldlist, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(fieldlist(1), mesh=mesh_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_src = ESMF_FieldCreate(mesh_src, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_src, farrayptr=dataPtr, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr(:) = 1.0_R8

          do n2 = 1,ncomps
             if ( n1 /= n2 .and. &
                  ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(n1,n2)) .and. &
                  is_local%wrap%med_coupling_active(n1,n2) ) then

                ! Get destination mesh
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n2), fieldlist=fieldlist, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(fieldlist(1), mesh=mesh_dst, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Createis_local%wrap%field_NormOne(n1,n2,m)
                do m = 1,nmappers
                   if (med_map_RH_is_created(is_local%wrap%RH,n1,n2,m,rc=rc)) then
                      is_local%wrap%field_NormOne(n1,n2,m) = ESMF_FieldCreate(mesh_dst, &
                           ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      call ESMF_FieldGet(is_local%wrap%field_NormOne(n1,n2,m), farrayptr=dataptr, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      dataptr(:) = czero
                      call med_map_field( &
                           field_src=field_src, &
                           field_dst=is_local%wrap%field_NormOne(n1,n2,m), &
                           routehandles=is_local%wrap%RH(n1,n2,:), &
                           maptype=m, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      if (mastertask) then
                         write(cn1,'(i1)') n1; write(cn2,'(i1)') n2; write(cm ,'(i1)') m
                         write(logunit,'(a)') trim(subname)//' created field_NormOne for '&
                              //compname(n1)//'->'//compname(n2)//' with mapping '//mapnames(m)
                      endif
                   end if
                end do ! end of loop over m mappers
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

  end subroutine med_map_mapnorm_init

  !================================================================================
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
       ! Loop over source field bundle
       do nf = 1, fieldCount
          ! Loop over the fldsSrc types
          do ns = 1,size(fldsSrc)
             ! Note that fieldnamelist is an array of names for the source fields
             ! The assumption is that there is only one mapping normalization
             ! for any given mapping type
             if ( fldsSrc(ns)%mapindex(destcomp) == mapindex .and. &
                  trim(fldsSrc(ns)%shortname) == trim(fieldnamelist(nf))) then
                ! Set the normalization to the input 
                packed_data(mapindex)%mapnorm = fldsSrc(ns)%mapnorm(destcomp)
             end if
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
          do ns = 1,size(fldsSrc)

             if ( fldsSrc(ns)%mapindex(destcomp) == mapindex .and. &
                  trim(fldsSrc(ns)%shortname) == trim(fieldnamelist(nf))) then

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
          packed_data(mapindex)%field_fracdst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
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
    use ESMF                  , only : ESMF_FieldRedist, ESMF_RouteHandle
    use esmFlds               , only : nmappers, mapfcopy
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
    integer                    :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
    real(r8), pointer          :: dataptr1d(:) => null()
    real(r8), pointer          :: dataptr2d(:,:) => null()
    real(r8), pointer          :: dataptr2d_packed(:,:) => null()
    type(ESMF_Field)           :: lfield
    type(ESMF_Field)           :: field_fracsrc
    type(ESMF_Field), pointer  :: fieldlist_src(:) => null()
    type(ESMF_Field), pointer  :: fieldlist_dst(:) => null()
    character(CL), allocatable :: fieldNameList(:)
    real(r8), pointer          :: data_norm(:) => null()
    real(r8), pointer          :: data_dst(:,:) => null()
    character(len=*), parameter  :: subname=' (med_map_field_packed) '
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
    end do ! end of loop over mapindex

    deallocate(fieldlist_src)
    deallocate(fieldlist_dst)

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
    real(r8), pointer :: data_src2d(:,:)    => null()
    real(r8), pointer :: data_dst2d(:,:)    => null()
    real(r8), pointer :: data_srctmp2d(:,:) => null()
    real(r8), pointer :: data_src1d(:)      => null()
    real(r8), pointer :: data_dst1d(:)      => null()
    real(r8), pointer :: data_srctmp1d(:)   => null()
    real(r8), pointer :: data_normsrc(:)    => null()
    real(r8), pointer :: data_normdst(:)    => null()
    integer           :: ungriddedUBound(1)     ! currently the size must equal 1 for rank 2 fields
    integer           :: lsize_src
    integer           :: lsize_dst
    character(len=*), parameter  :: subname=' (med_map_field_normalized) '
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

    use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF            , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_MAXSTR
    use ESMF            , only : ESMF_Field, ESMF_FieldRegrid
    use ESMF            , only : ESMF_TERMORDER_SRCSEQ, ESMF_Region_Flag, ESMF_REGION_TOTAL
    use ESMF            , only : ESMF_REGION_SELECT
    use ESMF            , only : ESMF_RouteHandle
    use esmFlds         , only : mapnstod_consd, mapnstod_consf, mapnstod_consd, mapnstod
    use esmFlds         , only : mapconsd, mapconsf
    use med_methods_mod , only : Field_diagnose => med_methods_Field_diagnose

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
    character(len=*), parameter :: subname='(med_map_field) '
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
    else
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RouteHandles(maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_map_field

  !================================================================================
  subroutine med_map_fb_field_regrid(FBin,fldin,FBout,fldout,RouteHandles,mapindex,rc)

    ! ----------------------------------------------
    ! Regrid a field in a field bundle to another field in a field bundle
    ! ----------------------------------------------

    use ESMF  , only : ESMF_FieldBundle, ESMF_RouteHandle, ESMF_Field
    use perf_mod , only : t_startf, t_stopf

    type(ESMF_FieldBundle), intent(in)           :: FBin
    character(len=*)      , intent(in)           :: fldin
    type(ESMF_FieldBundle), intent(inout)        :: FBout
    character(len=*)      , intent(in)           :: fldout
    type(ESMF_RouteHandle), intent(inout)        :: RouteHandles(:)
    integer               , intent(in)           :: mapindex
    integer               , intent(out)          :: rc
    ! ----------------------------------------------

    ! local
    type(ESMF_Field)       :: field1, field2
    character(CS)          :: lfldname
    character(len=*),parameter :: subname='(med_map_fb_field_regrid)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": start", ESMF_LOGMSG_INFO)
    endif

    call t_startf(subname)
    rc = ESMF_SUCCESS

    lfldname=trim(fldin)//'->'//trim(fldout)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    if (FB_FldChk(FBin , trim(fldin) , rc=rc) .and. FB_FldChk(FBout, trim(fldout), rc=rc)) then

       call FB_GetFieldByName(FBin, trim(fldin), field1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call FB_GetFieldByName(FBout, trim(fldout), field2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call med_map_field(field_src=field1, field_dst=field2, routehandles=routehandles, maptype=mapindex, &
            fldname=trim(lfldname), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//" field not found: "//&
            trim(fldin)//","//trim(fldout), ESMF_LOGMSG_INFO)
    endif

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf(subname)

  end subroutine med_map_fb_field_regrid

  !================================================================================
  subroutine med_map_uv_cart3d(usrc, vsrc, udst, vdst, RouteHandles, mapindex, rc)

    use ESMF, only : ESMF_Mesh, ESMF_MeshGet, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
    use ESMF, only : ESMF_Field, ESMF_FieldGet
    use ESMF, only : ESMF_FieldCreate, ESMF_FieldDestroy, ESMF_FieldRegrid
    use ESMF, only : ESMF_RouteHandle, ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL

    ! input/output variables
    type(ESMF_Field)       , intent(in)    :: usrc
    type(ESMF_Field)       , intent(in)    :: vsrc
    type(ESMF_Field)       , intent(inout) :: udst
    type(ESMF_Field)       , intent(inout) :: vdst
    type(ESMF_RouteHandle) , intent(inout) :: RouteHandles(:)
    integer                , intent(in)    :: mapindex
    integer                , intent(out)   :: rc

    ! local variables
    integer             :: n
    real(r8)            :: lon,lat
    real(r8)            :: coslon,coslat
    real(r8)            :: sinlon,sinlat
    real(r8)            :: ux,uy,uz
    type(ESMF_Mesh)     :: lmesh_src
    type(ESMF_Mesh)     :: lmesh_dst
    type(ESMF_Field)    :: field3d_src
    type(ESMF_Field)    :: field3d_dst
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
    logical             :: checkflag = .false.
    real(r8), parameter :: deg2rad = shr_const_pi/180.0_R8  ! deg to rads
    character(len=*), parameter :: subname=' med_map_uv_cart3d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Create two field bundles vec_src3d and vec_dst3d that contain
    ! all three fields with undistributed dimensions for each

    ! Get pointer to input u and v data source field data
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

    ! get source mesh and coordinates
    call ESMF_FieldGet(usrc, mesh=lmesh_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_src, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords_src(spatialDim*numOwnedElements))
    call ESMF_MeshGet(lmesh_src, ownedElemCoords=ownedElemCoords_src)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get destination mesh and coordinates
    call ESMF_FieldGet(udst, mesh=lmesh_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh_dst, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords_dst(spatialDim*numOwnedElements))
    call ESMF_MeshGet(lmesh_dst, ownedElemCoords=ownedElemCoords_dst)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! create source field and destination fields
    field3d_src = ESMF_FieldCreate(lmesh_src, ESMF_TYPEKIND_R8, name='src3d', &
         ungriddedLbound=(/1/), ungriddedUbound=(/3/), gridToFieldMap=(/2/), &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field3d_dst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, name='dst3d', &
         ungriddedLbound=(/1/), ungriddedUbound=(/3/), gridToFieldMap=(/2/), &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get pointers to source and destination data that will be filled in with rotation to cart3d
    call ESMF_FieldGet(field3d_src, farrayPtr=data2d_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field3d_dst, farrayPtr=data2d_dst, rc=rc)
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
    call med_map_field(&
         field_src=field3d_src, &
         field_dst=field3d_dst, &
         routehandles=routehandles, &
         maptype=mapindex, rc=rc)
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
    call ESMF_FieldDestroy(field3d_src, noGarbage=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field3d_dst, noGarbage=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_map_uv_cart3d

end module med_map_mod
