module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing glc export from mediator
  !-----------------------------------------------------------------------------

  use ESMF                , only : ESMF_FieldBundle
  use shr_nuopc_utils_mod , only : chkerr => shr_nuopc_utils_ChkErr
  use med_constants_mod   , only : R8, CS
  use med_constants_mod   , only : czero => med_constants_czero
  use med_constants_mod   , only : dbug_flag=>med_constants_dbug_flag

  implicit none
  private

  public  :: med_phases_prep_glc_avg
  public  :: med_phases_prep_glc_accum

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the FBLndAccum field bundle.

  type(ESMF_FieldBundle) :: FBLndAccum_lnd
  type(ESMF_FieldBundle) :: FBLndAccum_glc
  integer                :: FBLndAccumCnt

  ! Number of elevation classes
  integer :: nEC                     

  character(len=14) :: fldnames_fr_lnd(3) = (/'Flgl_qice_elev','Sl_tsrf_elev  ','Sl_topo_elev  '/)   
  character(len=14) :: fldnames_to_glc(3) = (/'Flgl_qice     ','Sl_tsrf       ','Sl_topo       '/)   

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_glc_accum(gcomp, rc)

    ! Carry out accumulation for the land-ice (glc) component
    ! Accumulation and averaging is done on the land input to the river component on the land grid
    ! Mapping from the land to the glc grid is then done with the time averaged fields

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_MeshLoc, ESMF_Mesh, ESMF_Field
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use esmFlds               , only : compglc, complnd
    use perf_mod              , only : t_startf, t_stopf
    use glc_elevclass_mod     , only : glc_get_num_elevation_classes

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,j,n,ncnt
    type(ESMF_MeshLoc)  :: meshloc
    type(ESMF_Mesh)     :: lmesh
    type(ESMF_Field)    :: lfield
    logical             :: first_call = .true.
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
    ! Create land accumulation field bundles on and and glc grid and initialize accumulation count
    !---------------------------------------

    if (first_call) then
       
       ! Determine number of elevation classes
       nEC = glc_get_num_elevation_classes()

       ! Create field bundles for the fldnames_fr_lnd that have an
       ! undistributed dimension corresponding to elevation classes

       call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(complnd,complnd), 1, lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh, meshloc=meshloc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       FBLndAccum = ESMF_FieldBundleCreate(name='FBLndAccum_lnd', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=flds_fr_lnd, &
            ungriddedLbound=1, ungriddedUbound=nec+1, gridToFieldMap=(/2/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBlndAccum_lnd, (/lfield/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_reset(FBLndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(compglc,compglc), 1, lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh, meshloc=meshloc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       FBGlcAccum = ESMF_FieldBundleCreate(name='FBLndAccum_glc', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=flds_fr_lnd, &
            ungriddedLbound=1, ungriddedUbound=nec+1, gridToFieldMap=(/2/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBlndAccum_glc, (/lfield/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_reset(FBLndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBLndAccumCnt = 0
      
       first_call = .false.
    end if

    !---------------------------------------
    ! Count the number of fields outside of scalar data
    !---------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,complnd))) then
       ncnt = 0
       call ESMF_LogWrite(trim(subname)//": FBImp(complnd,complnd) is not created", ESMF_LOGMSG_INFO)
    else 
       ! The scalar field has been removed from all mediator field bundles - so check if the fieldCount is
       ! 0 and not 1 here
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=ncnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBimp(complnd), returning", &
            ESMF_LOGMSG_INFO)
    end if

    !---------------------------------------
    ! accumulator land input to glc on land grid
    !---------------------------------------

    if (ncnt > 0) then

       call shr_nuopc_methods_FB_accum(FBLndAccum_lnd, is_local%wrap%FBImp(complnd,complnd), rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       FBLndAccumCnt = FBLndAccumCnt + 1

       if (dbug_flag > 10) then
          call shr_nuopc_methods_FB_diagnose(FBLndAccum_lnd, string=trim(subname)// ' FBLndAccum_lnd ',  rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       ! TODO: update local scalar data - set valid input flag to .false.
       !---------------------------------------

       call med_infodata_set_valid_glc_input(.false., med_infodata, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum

  !================================================================================================

  subroutine med_phases_prep_glc_avg(gcomp, rc)

    ! Prepares the GLC export Fields from the mediator

    use NUOPC                 , only : NUOPC_IsConnected
    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_Array, ESMF_ArrayGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundle, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_Mesh, ESMF_MeshGet
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_map_lnd2glc_mod   , only : med_map_lnd2glc
    use med_internalstate_mod , only : InternalState, logunit

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    real    , pointer      :: topolnd_g_ec(:,:)     ! topo in elevation classes
    real(r8), pointer      :: dataptr_g(:)          ! temporary data pointer for one elevation class
    real(r8), pointer      :: topoglc_g(:)          ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer      :: data_ice_covered_g(:) ! data for ice-covered regions on the GLC grid
    integer                :: nfld
    integer                :: ec
    integer                :: i,j,n,g,ncnt
    integer                :: lsize_g
    real(r8), pointer      :: glc_ice_covered(:)    ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer      :: glc_elevclass(:)      ! elevation classes glc grid
    real(r8), pointer      :: dataexp_g(:)          ! pointer into
    real(r8)               :: elev_l, elev_u        ! lower and upper elevations in interpolation range
    real(r8)               :: d_elev                ! elev_u - elev_l
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
       ! Average import from land accumuled FB
       !---------------------------------------

       call shr_nuopc_methods_FB_average(FBlndAccum_lnd, FBlndAccumCnt, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_diagnose(FBlndAccum_lnd, string=trim(subname)//&
            ' FBlndAccum for after avg for field bundle ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! Map accumulated field bundle from land grid (with elevation classes) to glc grid (without elevation classes)
       ! and set FBExp(compglc) data
       !---------------------------------------
       
       call med_map_lnd2glc(fldnames_fr_lnd, fldnames_to_glc, FBLndAcum_lnd, FBlndAccum_glc, gcomp, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! zero accumulator
       !---------------------------------------
       
       FBLndAccumCnt = 0

       call shr_nuopc_methods_FB_reset(FBLndAccum, value=czero, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! diagnostic output
       !---------------------------------------

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       ! update local scalar data - set valid input flag to .true.  TODO:
       !---------------------------------------

       call med_infodata_set_valid_glc_input(.true., med_infodata, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

end module med_phases_prep_glc_mod
