module med_map_lnd2glc_mod

  !---------------------------------------------------------------------
  ! Purpose:
  !
  ! This module contains routines for mapping fields from the LND grid (separated by GLC
  ! elevation class) onto the GLC grid
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
  !---------------------------------------------------------------------

#include "shr_assert.h"

  use ESMF                  , only : ESMF_FieldBundle
  use shr_sys_mod           , only : shr_sys_abort
  use shr_nuopc_utils_mod   , only : chkerr => shr_nuopc_utils_ChkErr
  use med_constants_mod     , only : R8, CS
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
  use esmFlds               , only : compglc, complnd, mapbilnr, mapconsf  

  implicit none
  private

  public  :: med_map_lnd2glc  ! map one field from LND -> GLC grid

  private :: get_glc_elevation_classes ! get the elevation class of each point on the glc grid
  private :: do_renormalize_smb  ! renormalizes surface mass balance

  ! needed for smb renormalization
  type(ESMF_FieldBundle) :: FBlnd_elev  ! contains fields with undistributed dimensions for elev classes
  type(ESMF_FieldBundle) :: FBglc_norm
  type(ESMF_FieldBundle) :: FBglc_icemask
  type(ESMF_FieldBundle) :: FBlnd_icemask
  type(ESMF_FieldBundle) :: FBglc_frac
  type(ESMF_FieldBundle) :: FBlnd_frac
  type(ESMF_FieldBundle) :: FBglc_other
  type(ESMF_FieldBundle) :: FBlnd_other

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize

  ! Number of elevation classes
  integer :: nEC                     

  ! Field names without elevation class (_elev) suffix
  character(len=*), parameter :: qice_fieldname   = 'Flgl_qice' ! Name of flux field giving surface mass balance
  character(len=*), parameter :: Sg_frac_field    = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_field    = 'Sg_topo'
  character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_map_lnd2glc(fldnames_fr_lnd, fldnames_to_glc, FBlndAccum, FBglcAccum, gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use glc_elevclass_mod     , only : glc_get_num_elevation_classes

    ! input/output variables
    character(len=*)       , intent(in)    :: fldnames_fr_lnd(:)
    character(len=*)       , intent(in)    :: fldnames_to_glc(:)
    type(ESMF_FieldBundle) , intent(inout) :: FBlndAccum
    type(ESMF_FieldBundle) , intent(inout) :: FBglcAccum
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    integer                , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    real(r8), pointer   :: topolnd_g_ec(:,:)     ! topo in elevation classes
    real(r8), pointer   :: dataptr_g(:)          ! temporary data pointer for one elevation class
    real(r8), pointer   :: topoglc_g(:)          ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer   :: data_ice_covered_g(:) ! data for ice-covered regions on the GLC grid
    real(r8), pointer   :: glc_ice_covered(:)    ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer   :: glc_elevclass(:)      ! elevation classes glc grid
    real(r8), pointer   :: dataexp_g(:)          ! pointer into
    real(r8), pointer   :: dataptr2d(:,:)
    real(r8), pointer   :: dataptr1d(:)
    real(r8), pointer   :: qice_g(:)             ! SMB (Flgl_qice) on glc grid without elev classes
    real(r8), pointer   :: qice_l(:,:)           ! SMB (Flgl_qice) on land grid with elev classes
    real(r8)            :: elev_l, elev_u        ! lower and upper elevations in interpolation range
    real(r8)            :: d_elev                ! elev_u - elev_l
    integer             :: nfld, ec
    integer             :: i,j,n,g,ncnt, lsize_g
    character(len=CS)   :: glc_renormalize_smb
    logical             :: glc_coupled_fluxes
    logical             :: lnd_prognostic
    logical             :: smb_renormalize
    logical             :: first_call = .true.
    character(len=*) , parameter   :: subname='(med_map_lnd2glc)'
    !---------------------------------------

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Initialize module variables on first call
    ! ------------------------------------------------------------------------

    if (first_call) then

       ! Determine number of elevation classes
       nEC = glc_get_num_elevation_classes()

       ! For now hardwire these
       ! TODO (mvertens, 2018-11-25) : put in the correct logic for the following rather than hardwiring
       glc_renormalize_smb = 'off'
       glc_coupled_fluxes  = .false.
       lnd_prognostic      = .false.
       smb_renormalize = do_renormalize_smb(glc_renormalize_smb, glc_coupled_fluxes, lnd_prognostic)

       if (smb_renormalize) then
          call shr_nuopc_methods_FB_init(FBlnd_elev, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), FBflds=is_local%wrap%FBImp(complnd,complnd), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! The following are fields with no undistributed dimension for the elevation classes

          call shr_nuopc_methods_FB_init(FBglc_icemask, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(compglc,compglc), fieldnameList=(/Sg_icemask_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return
          call shr_nuopc_methods_FB_init(FBlnd_icemask, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=(/Sg_icemask_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return

          call shr_nuopc_methods_FB_init(FBglc_frac, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(compglc,compglc), fieldNameList=(/Sg_frac_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return
          call shr_nuopc_methods_FB_init(FBlnd_frac, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=(/Sg_frac_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return

          call shr_nuopc_methods_FB_init(FBglc_other, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(compglc,compglc), fieldNameList=(/Sg_topo_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return
          call shr_nuopc_methods_FB_init(FBlnd_other, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=(/Sg_topo_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return

          ! Create a field bundle that will contain
          ! (fraction in this elevation class) x (icemask) for a given elevation class
          call shr_nuopc_methods_FB_init(FBglc_norm, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=(/'Sg_frac_times_icemask'/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return
       end if
       first_call = .false.
    end if

    ! ------------------------------------------------------------------------
    ! Map the accumulate land field from the land grid (in multiple elevation classes)
    ! to the glc grid (in multiple elevation classes) using bilinear interpolation
    ! ------------------------------------------------------------------------
    
    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle
    ! this case.
    
    call med_map_FB_Regrid_Norm(fldnames=fldnames_fr_lnd, &
         FBSrc=FBlndAccum, FBDst=FBglcAccum, &
         FBFrac=is_local%wrap%FBFrac(complnd), mapnorm='lfrac', &
         RouteHandle=is_local%wrap%RH(complnd,compglc,mapbilnr), &
         string='mapping normalized elevation class data from lnd to to glc', rc=rc)

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point on glc grid (output is topoglc_g)
    ! ------------------------------------------------------------------------

    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_frac_field, fldptr1=glc_ice_covered, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_topo_field, fldptr1=topoglc_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    lsize_g = size(glc_ice_covered)
    allocate(glc_elevclass(lsize_g))
    call get_glc_elevation_classes(glc_ice_covered, topoglc_g, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Determine topo field in multiple elevation classes on the glc grid
    ! ------------------------------------------------------------------------

    call shr_nuopc_methods_FB_getFldPtr(FBglcAccum, 'Sl_topo_elev', fldptr2=topolnd_g_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Loop over fields in export field bundle to GLC
    ! ------------------------------------------------------------------------

    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle this case.

    allocate(data_ice_covered_g (lsize_g))
    do nfld = 1, size(fldnames_to_glc)

       ! ------------------------------------------------------------------------
       ! Perform vertical interpolation of data onto ice sheet topography
       ! This maps all of the input elevation classes into an export to glc without elevation classes
       ! ------------------------------------------------------------------------

       ! Get a pointer to the land data in multiple elevation classes on the glc grid
       call shr_nuopc_methods_FB_getFldPtr(FBglcAccum, fldnames_fr_lnd(nfld), fldptr2=dataptr2d, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! Get a pointer to the data for the field that will be sent to glc (without elevation classes)
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBExp(compglc), fldnames_to_glc(nfld), fldptr1=dataexp_g, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       data_ice_covered_g(:) = 0._r8
       do n = 1, lsize_g

          ! For each ice sheet point, find bounding EC values...

          if (topoglc_g(n) < topolnd_g_ec(1,n)) then

             ! lower than lowest mean EC elevation value
             data_ice_covered_g(n) = dataptr2d(1,n)

          else if (topoglc_g(n) >= topolnd_g_ec(nEC, n)) then

             ! higher than highest mean EC elevation value
             data_ice_covered_g(n) = dataptr2d(nEC,n)

          else

             ! do linear interpolation of data in the vertical
             do ec = 2, nEC ! TODO: are these the right bounds?
                if (topoglc_g(n) < topolnd_g_EC(ec,n)) then
                   elev_l = topolnd_g_EC(ec-1,n)
                   elev_u = topolnd_g_EC(ec  ,n)
                   d_elev = elev_u - elev_l
                   if (d_elev <= 0) then
                      ! This shouldn't happen, but handle it in case it does. In this case,
                      ! let's arbitrarily use the mean of the two elevation classes, rather
                      ! than the weighted mean.
                      write(logunit,*) subname//' WARNING: topo diff between elevation classes <= 0'
                      write(logunit,*) 'n, ec, elev_l, elev_u = ', n, ec, elev_l, elev_u
                      write(logunit,*) 'Simply using mean of the two elevation classes,'
                      write(logunit,*) 'rather than the weighted mean.'
                      data_ice_covered_g(n) = dataptr2d(ec-1,n) * 0.5_r8 &
                           + dataptr2d(ec  ,n) * 0.5_r8
                   else

                      data_ice_covered_g(n) =  dataptr2d(ec-1,n) * (elev_u - topoglc_g(n)) / d_elev  &
                                             + dataptr2d(ec  ,n) * (topoglc_g(n) - elev_l) / d_elev
                   end if

                   exit
                end if
             end do
          end if  ! topoglc_g(n)

          if (glc_elevclass(n) /= 0) then
             ! ice-covered cells have interpolated values
             dataexp_g(n) = data_ice_covered_g(n)
          else
             ! non ice-covered cells have bare land value
             dataexp_g(n) = real(dataptr2d(1,n))
          end if
       end do  ! lsize_g

       ! ------------------------------------------------------------------------
       ! Renormalize surface mass balance (smb, here named dataexp_g) so that the global
       ! integral on the glc grid is equal to the global integral on the land grid.
       ! ------------------------------------------------------------------------

       ! No longer need to make a preemptive adjustment to qice_g to account for area differences 
       ! between CISM and the coupler. In NUOPC, the area correction is done in! the cap not in the 
       ! mediator, so to preserve the bilinear mapping values, do not need to do any area correction 
       ! scaling in the CISM NUOPC cap

       if (trim(fldnames_to_glc(nfld)) == trim(qice_fieldname) .and. smb_renormalize) then
          call shr_nuopc_methods_FB_getFldPtr(FBlndAccum, trim(qice_fieldname)//'_elev', fldptr2=qice_l, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBExp(compglc), qice_fieldname, fldptr1=qice_g, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return

          call prep_glc_renormalize_smb(gcomp, qice_l, qice_g, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end do  ! end of loop over fields
    deallocate(data_ice_covered_g)

  end subroutine med_map_lnd2glc

  !================================================================================================

  subroutine prep_glc_renormalize_smb(gcomp, qice_l, qice_g, rc)

    !------------------
    ! Renormalizes surface mass balance (smb, here named qice_g) so that the global
    ! integral on the glc grid is equal to the global integral on the land grid.
    !
    ! (1) Map Sg_icemask from the glc grid to the land grid.
    !     Because of coupler lags, the current Sg_icemask_l might not be up to date with Sg_icemask_g.
    ! (2) Map Sg_ice_covered from the glc grid to the land grid
    !     This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    !
    ! This is required for conservation - although conservation is only necessary if we
    ! are running with a fully-interactive, two-way-coupled glc.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
    !------------------

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_VMAllReduce, ESMF_REDUCE_SUM
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet
    use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_Array, ESMF_ArrayGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_map_glc2lnd_mod   , only : med_map_glc2lnd_elevclass
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)    :: gcomp
    real(r8) , pointer     :: qice_g(:)   ! SMB (Flgl_qice) on glc grid without elev classes
    real(r8) , pointer     :: qice_l(:,:) ! SMB (Flgl_qice) on land grid with elev classes
    integer  , intent(out) :: rc          ! return error code

    ! local variables
    ! Note: Sg_icemask defines where the ice sheet model can receive a nonzero SMB from the land model.
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    type(ESMF_Mesh)     :: lmesh
    type(ESMF_Field)    :: lfield
    type(ESMF_Array)    :: larray
    real(r8), pointer   :: aream_l(:)      ! cell areas on land grid, for mapping
    real(r8), pointer   :: aream_g(:)      ! cell areas on glc grid, for mapping
    real(r8), pointer   :: frac_l(:,:)     ! EC fractions (Sg_ice_covered) on land grid
    real(r8), pointer   :: Sg_icemask_g(:) ! icemask on glc grid
    real(r8), pointer   :: Sg_icemask_l(:) ! icemask on land grid
    real(r8), pointer   :: lfrac(:)        ! land fraction on land grid
    integer             :: ec              ! loop index over elevation classes
    integer             :: n

    ! local and global sums of accumulation and ablation; used to compute renormalization factors
    real(r8), target :: local_accum_on_lnd_grid(1), global_accum_on_lnd_grid(1)
    real(r8), target :: local_accum_on_glc_grid(1), global_accum_on_glc_grid(1)
    real(r8), target :: local_ablat_on_lnd_grid(1), global_ablat_on_lnd_grid(1)
    real(r8), target :: local_ablat_on_glc_grid(1), global_ablat_on_glc_grid(1)

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor ! ratio between global ablation on the two grids

    real(r8) :: effective_area      ! grid cell area multiplied by min(lfrac,Sg_icemask_l).
    ! This is the area that can contribute SMB to the ice sheet model.
    character(len=*), parameter  :: subname='(prep_glc_renormalize_smb)'
    !---------------------------------------------------------------

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
    ! Get vm
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Determine areas on glc and land grid and land fraction on land grid
    !---------------------------------------

    ! determine fraction on land grid
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBFrac(complnd), 'lfrac', fldptr1=lfrac, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! determine areas on land grid 
    call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(complnd,complnd), fieldnum=1, field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh, elemAreaArray=larray, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ArrayGet(larray, farrayptr=aream_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine areas on glc grid 
    call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(compglc,compglc), fieldnum=1, field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh, elemAreaArray=larray, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ArrayGet(larray, farrayptr=aream_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map Sg_icemask from the glc grid to the land grid.
    !---------------------------------------

    ! This may not be necessary, if Sg_icemask_l has already been mapped from Sg_icemask_g.
    ! It is done here for two reasons:
    ! (1) The mapping will *not* have been done if we are running with dlnd (e.g., a TG case).
    ! (2) Because of coupler lags, the current Sg_icemask_l might not be up to date with
    !     Sg_icemask_g. This probably isn't a problem in practice, but doing the mapping
    !     here ensures the mask is up to date.
    !
    ! This mapping uses the same options as the standard glc -> lnd mapping done in
    ! prep_lnd_calc_g2x_lx. If that mapping ever changed (e.g., changing norm to
    ! .false.), then we should change this mapping, too.
    !
    ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but this
    ! requires some more thought

    call med_map_FB_Regrid_Norm((/Sg_icemask_field/), &
         is_local%wrap%FBImp(compglc,compglc), FBlnd_icemask, &
         is_local%wrap%FBNormOne(compglc,compglc,mapconsf ), 'one', &
         is_local%wrap%RH(compglc, complnd, mapconsf), &
         string='mapping Sg_imask_g to Sg_imask_l (from glc to land)', rc=rc)

    ! Note that Sg_icemask_l will not have elevation classes
    call shr_nuopc_methods_FB_getFldPtr(FBlnd_icemask, Sg_icemask_field, fldptr1=Sg_icemask_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map Sg_ice_covered from the glc grid (no elevation classes) to the land grid (with elevation classes)
    !---------------------------------------

    ! This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    ! These fields are needed to integrate the total SMB on the land grid, for conservation purposes.
    ! As above, the mapping may not be necessary, because Sg_ice_covered might already have been mapped.
    ! However, the mapping will not have been done in a TG case with dlnd, and it might not
    ! be up to date because of coupler lags (though the latter probably isn't a problem
    ! in practice).
    !
    ! Note that, for a case with full two-way coupling, we will only conserve if the
    ! actual land cover used over the course of the year matches these currently-remapped
    ! values. This should generally be the case with the current coupling setup.
    !
    ! One could argue that it would be safer (for conservation purposes) if LND sent its
    ! grid cell average SMB values, or if it sent its own notion of the area in each
    ! elevation class for the purpose of creating grid cell average SMB values here. But
    ! these options cause problems if we're not doing full two-way coupling (e.g., in a TG
    ! case with dlnd, or in the common case where GLC is a diagnostic component that
    ! doesn't cause updates in the glacier areas in LND). In these cases without full
    ! two-way coupling, if we use the LND's notion of the area in each elevation class,
    ! then the conservation corrections would end up correcting for discrepancies in
    ! elevation class areas between LND and GLC, rather than just correcting for
    ! discrepancies arising from the remapping of SMB. (And before you get worried: It
    ! doesn't matter that we are not conserving in these cases without full two-way
    ! coupling, because GLC isn't connected with the rest of the system in terms of energy
    ! and mass in these cases. So in these cases, it's okay that the LND integral computed
    ! here differs from the integral that LND itself would compute.)

    ! Map frac_field on glc grid without elevation classes to frac_field on land grid with elevation classes
    call med_map_glc2lnd_elevclass(FBglc_icemask, FBglc_norm, &
         FBglc_frac, FBlnd_frac, FBglc_other, FBlnd_other, &
         icemask_field=Sg_icemask_field, &
         frac_field=Sg_frac_field, &
         topo_field=Sg_topo_field, &
         fields_to_map=(/Sg_topo_field/), &
         ntopo_field=1, &
         RouteHandle=is_local%wrap%RH(compglc,complnd,mapconsf), &
         FBlnd_elev=FBlnd_elev, rc=rc)

    ! Note: frac_l comes from med_map_glc2lnd_elevclass
    call shr_nuopc_methods_FB_getFldPtr(FBlnd_elev, trim(Sg_frac_field)//'_elev', fldptr2=frac_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Note: qice_l comes from FBlndAccum

    ! Sum qice_l over all elevation classes for each local land grid cell then do a global sum

    local_accum_on_lnd_grid(1) = 0.0_r8
    local_ablat_on_lnd_grid(1) = 0.0_r8
    do n = 1, size(lfrac)
       ! Calculate effective area for sum -  need the mapped Sg_icemask_l
       effective_area = min(lfrac(n),Sg_icemask_l(n)) * aream_l(n)

       do ec = 1, nEC+1  
          if (qice_l(ec,n) >= 0.0_r8) then
             local_accum_on_lnd_grid(1) = local_accum_on_lnd_grid(1) + effective_area * frac_l(ec,n) * qice_l(ec,n)
          else
             local_ablat_on_lnd_grid(1) = local_ablat_on_lnd_grid(1) + effective_area * frac_l(ec,n) * qice_l(ec,n)
          endif
       enddo  ! ec
    enddo  ! n
    call ESMF_VMAllreduce(vm, senddata=local_accum_on_lnd_grid, recvdata=global_accum_on_lnd_grid,&
         count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    call ESMF_VMAllreduce(vm, senddata=local_ablat_on_lnd_grid, recvdata=global_ablat_on_lnd_grid,&
         count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)

    ! Sum qice_g over local glc grid cells.  

    ! TODO: is the following a problem
    ! Note: This sum uses the coupler areas (aream_g), which differ from the native CISM areas.
    ! But since the original qice_g (from bilinear remapping) has been multiplied by
    ! area_g/aream_g above, this calculation is equivalent to multiplying the original qice_g
    ! by the native CISM areas (area_g).
    ! If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    ! then it would be appropriate to use the native CISM areas in this sum.

    ! Determine Sg_icemask_g from the export FB to glc
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBexp(compglc), Sg_icemask_field, fldptr1=Sg_icemask_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    local_accum_on_glc_grid(1) = 0.0_r8
    local_ablat_on_glc_grid(1) = 0.0_r8
    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          local_accum_on_glc_grid(1) = local_accum_on_glc_grid(1) + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       else
          local_ablat_on_glc_grid(1) = local_ablat_on_glc_grid(1) + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       endif
    enddo  ! n
    call ESMF_VMAllreduce(vm, senddata=local_accum_on_glc_grid, recvdata=global_accum_on_glc_grid,&
         count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    call ESMF_VMAllreduce(vm, senddata=local_ablat_on_glc_grid, recvdata=global_ablat_on_glc_grid,&
         count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)

    ! Renormalize
    if (global_accum_on_glc_grid(1) > 0.0_r8) then
       accum_renorm_factor = global_accum_on_lnd_grid(1) / global_accum_on_glc_grid(1)
    else
       accum_renorm_factor = 0.0_r8
    endif

    if (global_ablat_on_glc_grid(1) < 0.0_r8) then  ! negative by definition
       ablat_renorm_factor = global_ablat_on_lnd_grid(1) / global_ablat_on_glc_grid(1)
    else
       ablat_renorm_factor = 0.0_r8
    endif

    if (mastertask) then
       write(logunit,*) 'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          qice_g(n) = qice_g(n) * accum_renorm_factor
       else
          qice_g(n) = qice_g(n) * ablat_renorm_factor
       endif
    enddo

  end subroutine prep_glc_renormalize_smb

  !-----------------------------------------------------------------------
  subroutine get_glc_elevation_classes(glc_ice_covered, glc_topo, glc_elevclass)

    !------------------
    ! Get the elevation class of each point on the glc grid.
    ! For grid cells that are ice-free, the elevation class is set to 0.
    ! All arguments (glc_ice_covered, glc_topo and glc_elevclass) must be the same size.
    !------------------

    use med_internalstate_mod , only : logunit
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
    use glc_elevclass_mod     , only : glc_get_elevation_class

    ! input/output variables
    real(r8), intent(in)  :: glc_ice_covered(:) ! ice-covered (1) vs. ice-free (0)
    real(r8), intent(in)  :: glc_topo(:)        ! ice topographic height
    integer , intent(out) :: glc_elevclass(:)   ! elevation class

    ! local variables
    integer :: npts
    integer :: glc_pt
    integer :: err_code

    ! Tolerance for checking whether ice_covered is 0 or 1
    real(r8), parameter :: ice_covered_tol = 1.e-13

    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT_FL((size(glc_ice_covered) == npts), __FILE__, __LINE__)
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

    do glc_pt = 1, npts
       if (abs(glc_ice_covered(glc_pt) - 1._r8) < ice_covered_tol) then
          ! This is an ice-covered point

          call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)
          if ( err_code == GLC_ELEVCLASS_ERR_NONE .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_LOW .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_HIGH) then
             ! These are all acceptable "errors" - it is even okay for these purposes if
             ! the elevation is lower than the lower bound of elevation class 1, or
             ! higher than the upper bound of the top elevation class.

             ! Do nothing
          else
             write(logunit,*) subname, ': ERROR getting elevation class for ', glc_pt
             write(logunit,*) glc_errcode_to_string(err_code)
             call shr_sys_abort(subname//': ERROR getting elevation class')
          end if
       else if (abs(glc_ice_covered(glc_pt) - 0._r8) < ice_covered_tol) then
          ! This is a bare land point (no ice)
          glc_elevclass(glc_pt) = 0
       else
          ! glc_ice_covered is some value other than 0 or 1
          ! The lnd -> glc downscaling code would need to be reworked if we wanted to
          ! handle a continuous fraction between 0 and 1.
          write(logunit,*) subname, ': ERROR: glc_ice_covered must be 0 or 1'
          write(logunit,*) 'glc_pt, glc_ice_covered = ', glc_pt, glc_ice_covered(glc_pt)
          call shr_sys_abort(subname//': ERROR: glc_ice_covered must be 0 or 1')
       end if
    end do

  end subroutine get_glc_elevation_classes

  !================================================================================================

  function do_renormalize_smb(glc_renormalize_smb, glc_coupled_fluxes, lnd_prognostic) &
       result(smb_renormalize)

    ! Returns a logical saying whether we should do the smb renormalization

    ! function return
    logical :: smb_renormalize   ! function return value

    ! input/output variables
    character(len=*), intent(in) :: glc_renormalize_smb  ! namelist option saying whether to do smb renormalization
    logical         , intent(in) :: glc_coupled_fluxes   ! does glc send fluxes to other components?
    logical         , intent(in) :: lnd_prognostic       ! is lnd a prognostic component?

    ! local variables
    character(len=*), parameter :: subname = '(do_renormalize_smb)'
    !---------------------------------------------------------------

    select case (glc_renormalize_smb)
    case ('on')
       smb_renormalize = .true.

    case ('off')
       smb_renormalize = .false.

    case ('on_if_glc_coupled_fluxes')
       if (.not. lnd_prognostic) then
          ! Do not renormalize if running glc with dlnd (T compsets): In this case
          ! there is no feedback from glc to lnd, and conservation is not important
          smb_renormalize = .false.
       else if (.not. glc_coupled_fluxes) then
          ! Do not renormalize if glc does not send fluxes to other components: In this
          ! case conservation is not important
          smb_renormalize = .false.
       else
          ! lnd_prognostic is true and glc_coupled_fluxes is true
          smb_renormalize = .true.
       end if

    case default
       write(logunit,*) subname,' ERROR: unknown value for glc_renormalize_smb: ', &
            trim(glc_renormalize_smb)
       call shr_sys_abort(subname//' ERROR: unknown value for glc_renormalize_smb')

    end select

  end function do_renormalize_smb

end module med_map_lnd2glc_mod
