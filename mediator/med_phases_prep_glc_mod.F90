module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing glc export from mediator
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_FieldBundle
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleIsCreated
  use ESMF                  , only : ESMF_MeshLoc, ESMF_Mesh, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use esmFlds               , only : compglc, complnd, mapconsf, mapbilnr
  use esmFlds               , only : shr_nuopc_fldlist_type
  use med_constants_mod     , only : R8, CS
  use med_constants_mod     , only : czero => med_constants_czero
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, logunit
  use med_map_lnd2glc_mod   , only : med_map_lnd2glc
  use shr_nuopc_utils_mod   , only : chkerr => shr_nuopc_utils_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_prep_glc_init
  public :: med_phases_prep_glc_accum
  public :: med_phases_prep_glc_avg

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the FBLndAccum field bundle.

  type(ESMF_FieldBundle) :: FBLndAccum_lnd
  type(ESMF_FieldBundle) :: FBLndAccum_glc
  integer                :: FBLndAccumCnt
  type(shr_nuopc_fldlist_type) :: fldlist_lnd2glc

  character(len=14) :: fldnames_fr_lnd(3) = (/'Flgl_qice_elev','Sl_tsrf_elev  ','Sl_topo_elev  '/)
  character(len=14) :: fldnames_to_glc(2) = (/'Flgl_qice     ','Sl_tsrf       '/)

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_glc_init(gcomp, rc)

    !---------------------------------------
    ! Create land accumulation field bundles on and and glc grid and initialize accumulation count
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,n,ncnt
    type(ESMF_Mesh)     :: lmesh
    type(ESMF_Field)    :: lfield
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    integer             :: nec ! Number of elevation classes
    character(len=*),parameter  :: subname='(med_phases_prep_glc_init)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Count the number of fields outside of scalar data
    !---------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,complnd))) then
       ncnt = 0
       call ESMF_LogWrite(trim(subname)//": FBImp(complnd,complnd) is not created", ESMF_LOGMSG_INFO)
    else
       ! The scalar field has been removed from all mediator field bundles - so determine ncnt for below
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=ncnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Create field bundles for the fldnames_fr_lnd that have an
    ! undistributed dimension corresponding to elevation classes
    ! nec is the size of the undistributed dimension (number of elevation classes)
    !---------------------------------------

    if (ncnt > 0) then

       ! Create accumulation field bundle from land on the land grid
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(1), fldptr2=data2d_in, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       nec = size(data2d_in, dim=1)
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBLndAccum_lnd = ESMF_FieldBundleCreate(name='FBLndAccum_lnd', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/nec/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum_lnd, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
       call shr_nuopc_methods_FB_reset(FBLndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create accumulation field bundle from land on the glc grid
       ! Determine glc mesh from the mesh from the first export field to glc
       ! However FBlndAccum_glc has the fields fldnames_fr_lnd BUT ON the glc grid

       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fldnames_to_glc(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBlndAccum_glc = ESMF_FieldBundleCreate(name='FBLndAccum_glc', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/nec/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum_glc, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
       call shr_nuopc_methods_FB_reset(FBLndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBLndAccumCnt = 0

       allocate(fldlist_lnd2glc%flds(3))
       fldlist_lnd2glc%flds(1)%shortname = trim(fldnames_fr_lnd(1))
       fldlist_lnd2glc%flds(2)%shortname = trim(fldnames_fr_lnd(2))
       fldlist_lnd2glc%flds(3)%shortname = trim(fldnames_fr_lnd(3))
       fldlist_lnd2glc%flds(1)%mapindex(compglc) = mapconsf
       fldlist_lnd2glc%flds(2)%mapindex(compglc) = mapbilnr
       fldlist_lnd2glc%flds(3)%mapindex(compglc) = mapbilnr
       fldlist_lnd2glc%flds(1)%mapnorm(compglc) = 'lfrac'
       fldlist_lnd2glc%flds(2)%mapnorm(compglc) = 'none'
       fldlist_lnd2glc%flds(3)%mapnorm(compglc) = 'none'

    end if


    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_init

  !================================================================================================

  subroutine med_phases_prep_glc_accum(gcomp, rc)

    ! Carry out accumulation for the land-ice (glc) component
    ! Accumulation and averaging is done on the land input to the river component on the land grid
    ! Mapping from the land to the glc grid is then done with the time averaged fields

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,n,ncnt
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    character(len=*),parameter  :: subname='(med_phases_prep_glc_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Count the number of fields outside of scalar data
    !---------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,complnd))) then
       ncnt = 0
       call ESMF_LogWrite(trim(subname)//": FBImp(complnd,complnd) is not created", ESMF_LOGMSG_INFO)
    else
       ! The scalar field has been removed from all mediator field bundles - so determine ncnt for below
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=ncnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! accumulator land input to glc on land grid
    !---------------------------------------

    if (ncnt > 0) then

       do n = 1, size(fldnames_fr_lnd)
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), &
               fldnames_fr_lnd(n), fldptr2=data2d_in, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(FBLndAccum_lnd, fldnames_fr_lnd(n), fldptr2=data2d_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do i = 1,size(data2d_out, dim=2)
             data2d_out(:,i) = data2d_out(:,i) + data2d_in(:,i)
          end do
       end do

       FBLndAccumCnt = FBLndAccumCnt + 1

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(FBLndAccum_lnd, string=trim(subname)// ' FBLndAccum_lnd ',  rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum

  !================================================================================================

  subroutine med_phases_prep_glc_avg(gcomp, rc)

    ! Prepares the GLC export Fields from the mediator

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    integer                :: i, n, ncnt            ! counters
    real(r8), pointer      :: data2d(:,:)
    character(len=*) , parameter   :: subname='(med_phases_prep_glc_avg)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=ncnt, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt ==  0) then
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBExp(compglc), returning", &
            ESMF_LOGMSG_INFO)

    else

       !---------------------------------------
       ! Average import from accumulated land import FB
       !---------------------------------------

       do n = 1, size(fldnames_fr_lnd)
          call shr_nuopc_methods_FB_GetFldPtr(FBLndAccum_lnd, fldnames_fr_lnd(n), fldptr2=data2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          data2d(:,:) = data2d(:,:) / real(FBLndAccumCnt)
       end do

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(FBlndAccum_lnd, string=trim(subname)//&
               ' FBlndAccum for after avg for field bundle ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       ! Map accumulated field bundle from land grid (with elevation classes) to glc grid (without elevation classes)
       ! and set FBExp(compglc) data
       !---------------------------------------

       call shr_nuopc_methods_FB_reset(FBLndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call med_map_lnd2glc(fldlist_lnd2glc, fldnames_fr_lnd, fldnames_to_glc, FBLndAccum_lnd, FBlndAccum_glc, gcomp, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       ! zero accumulator and accumulated field bundles
       !---------------------------------------

       FBLndAccumCnt = 0

       call shr_nuopc_methods_FB_reset(FBLndAccum_lnd, value=czero, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_reset(FBLndAccum_glc, value=czero, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! update local scalar data - set valid input flag to .true.  TODO:
       !---------------------------------------

       !call med_infodata_set_valid_glc_input(.true., med_infodata, rc=rc) ! TODO: fix this
       !if (chkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

end module med_phases_prep_glc_mod
