module med_map_glc2lnd_mod

  !---------------------------------------------------------------------
  ! This module contains routines for mapping fields from the GLC grid onto the LND grid
  ! (separated by GLC elevation class)
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1sjsaiPYsPJ9A7dVGJIHGg4rVIY2qF5aRXbNzSXVAafU/edit?usp=sharing
  !---------------------------------------------------------------------

#include "shr_assert.h"

  use med_constants_mod   , only : CS
  use med_constants_mod   , only : R8
  use med_constants_mod   , only : dbug_flag=>med_constants_dbug_flag
  use shr_nuopc_utils_mod , only : chkerr=>shr_nuopc_utils_ChkErr

  implicit none
  private

  public  :: med_map_glc2lnd_elevclass ! map all fields from GLC -> LND grid will be separated by elevation class
  private :: get_glc_elevation_classes ! get elevation class of each glc cell

  ! private module variables
  integer :: num_elevation_classes
  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_map_glc2lnd_elevclass(FBglc_icemask, FBglc_norm, &
       FBglc_frac, FBlnd_frac, FBglc_other, FBlnd_other, &
       icemask_field, frac_field, topo_field, fields_to_map, ntopo_field, RouteHandle, FBlnd_elev, rc)

    !------------------
    ! Maps fields from the GLC grid to the LND grid.
    ! On the GLC grid the fields will not have elevation classes.
    ! On the LND grid they will have elevation classes.
    !
    ! Maps frac_field, topo_field, plus all fields defined in extra_fields (a character array)
    ! without the '_elev' suffix.
    !
    ! Assumes that 
    ! - FBglc_icemask contains icemask_field (NOT mapped here, but needed as an input to the mapping)
    ! - FBglc_frac  and FBlnd_frac  contain frac_field
    ! - FBglc_other and FBlnd_other contain fields_to_map (topo_field plus each field in extra_fields)
    !
    ! Assumes that FBlnd (on land grid) contains:
    ! - frac_field (with multiple elevation classes in undistributed dimension)
    ! - topo_field (with multiple elevation classes in undistributed dimension)
    ! - And similarly for each field in extra_fields
    !------------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_RouteHandle
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use glc_elevclass_mod     , only : glc_mean_elevation_virtual
    use glc_elevclass_mod     , only : glc_get_num_elevation_classes
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)           :: FBglc_icemask
    type(ESMF_FieldBundle) , intent(inout)        :: FBglc_norm 
    type(ESMF_FieldBundle) , intent(inout)        :: FBglc_frac 
    type(ESMF_FieldBundle) , intent(inout)        :: FBlnd_frac 
    type(ESMF_FieldBundle) , intent(inout)        :: FBglc_other
    type(ESMF_FieldBundle) , intent(inout)        :: FBlnd_other
    character(len=*)       , intent(in)           :: icemask_field    ! name  of field FBglc_icemask containing ice mask
    character(len=*)       , intent(in)           :: frac_field       ! name  of field FBglc_frac    containing ice fraction
    character(len=*)       , intent(in)           :: topo_field       ! name  of field FBglc_topo    containing topo
    character(len=*)       , intent(in)           :: fields_to_map(:) ! names of fields_to_map other than frac_field
    integer                , intent(in)           :: ntopo_field      ! index of topo_field in fields_to_map(:)
    type(ESMF_RouteHandle) , intent(inout)        :: RouteHandle      ! glc->lnd mapping route handle
    type(ESMF_FieldBundle) , intent(inout)        :: FBlnd_elev       ! output field bundle with 2d fields for elev classes
    integer                , intent(out)          :: rc

    ! local variables
    integer                :: ec, nfld, n
    real(r8)               :: topo_virtual
    integer  , allocatable :: glc_elevclass(:)                    ! elevation class of each glc cell (assuming cell is ice-covered)
    real(r8) , pointer     :: glc_icemask_g(:)                    ! glc ice mask field on glc grid
    real(r8) , pointer     :: glc_frac_g(:)                       ! total ice fraction in each glc cell
    real(r8) , pointer     :: glc_frac_g_ec(:)                    ! glc fractions in one elev class, on the glc grid
    real(r8) , pointer     :: glc_frac_l_ec(:)                    ! glc fractions in one elev class mapped to the land grid
    real(r8) , pointer     :: glc_frac_times_icemask_g_ec(:)      ! (glc fraction in one elev class) x (icemask), on the glc grid
    real(r8) , pointer     :: glc_topo_g(:)                       ! topographic height of each glc cell
    real(r8) , pointer     :: dataptr2d(:,:)
    real(r8) , pointer     :: dataptr1d(:)
    integer                :: lsize_g 
    character(len=3)       :: cvalue
    logical                :: first_call = .true.
    character(len=*), parameter :: subname = 'map_glc2lnd_elevclass'
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    ! ------------------------------------------------------------------------
    ! Set module variables on first call
    ! ------------------------------------------------------------------------

    if (first_call) then
       num_elevation_classes = glc_get_num_elevation_classes()
       first_call = .false.
    end if

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point on the glc grid (glc_elevclass)
    ! ------------------------------------------------------------------------

    ! glc_frac is the total ice fraction in each glc gridcell
    call shr_nuopc_methods_FB_getFldPtr(FBglc_frac, trim(frac_field), fldptr1=glc_frac_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! glc_topo is the topographic height of each glc gridcell
    call shr_nuopc_methods_FB_getFldPtr(FBglc_other, trim(topo_field), fldptr1=glc_topo_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! given glc_topo for each glc grid cell, set the glc elevation classes, glc_elevclass
    ! for each grid cell
    lsize_g = size(glc_frac_g)

    allocate(glc_elevclass(lsize_g))
    call get_glc_elevation_classes(glc_topo_g, glc_elevclass)

    !---------------------------------
    ! Loop over all elevation classes
    !---------------------------------
    do ec = 0, num_elevation_classes

       !---------------------------------
       ! Determine fractional ice coverage for each elevation class on the land grid
       !---------------------------------

       ! set pointer to array holding glc fraction in one elev class, on the glc grid
       ! the values in frac_field_glc will be reset for each loop iteration over elevation classes
       call shr_nuopc_methods_FB_getFldPtr(FBglc_frac, frac_field, fldptr1=glc_frac_g_ec, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! set fractional ice coverage for a given elevation class (THIS IS WHAT WILL GET MAPPED)
       ! this will set the array glc_frac_this_ec_g, ice fractions in this elevation class, for each glc gridcell.
       ! note that glc_elevclass gives the elevation class of each glc grid cell, assuming that
       ! the grid cell is ice-covered.
       if (ec == 0) then
          glc_frac_g_ec(:) = 1._r8 - glc_frac_g(:)
       else
          where (glc_elevclass == ec)
             glc_frac_g_ec = glc_frac_g
          elsewhere
             glc_frac_g_ec = 0._r8
          end where
       end if

       ! map fraction in each elevation class from the glc grid to the land grid
       ! Note that FBlnd_frac does not contain an undistributed dimension
       write(cvalue,*) ec
       call med_map_FB_Regrid_Norm( &
            fldnames=(/frac_field/),&
            FBSrc=FBglc_frac, &
            FBDst=FBlnd_frac, &
            FBfrac=FBglc_icemask, mapnorm=trim(icemask_field), &
            RouteHandle=RouteHandle, &
            string='mapping elevation class fraction from glc to land elev class '//trim(cvalue), rc=rc)

       ! set pointer to array holding glc fraction in one elev class, on the land grid
       ! the values in glc_frac_l_ec will be reset for each loop iteration over elevation classes
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_frac, trim(frac_field), fldptr1=glc_frac_l_ec, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! frac_field_lnd2d and topo_field_lnd2d are 2d arrays where the elevation class is the 
       ! inner most index
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_elev, trim(frac_field), fldptr2=dataptr2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       dataptr2d(ec+1,:) = glc_frac_l_ec(:) 
       
       !---------------------------------
       ! Map other fields to the land grid
       !---------------------------------

       ! Fill in data for FBglc_norm
       ! Only grid cells that are both (a) within the icemask and (b) in this elevation class
       ! will be included in the following mapping.

       call shr_nuopc_methods_FB_getFldPtr(FBglc_norm, 'Sg_frac_times_icemask', fldptr1=glc_frac_times_icemask_g_ec, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Determine ice mask field on glc grid (glc_icemask)

       call shr_nuopc_methods_FB_getFldPtr(FBglc_icemask, trim(icemask_field), fldptr1=glc_icemask_g)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       glc_frac_times_icemask_g_ec(:) = glc_frac_g(:) * glc_icemask_g(:)

       ! Map all other fields normalizing by glc_frac_this_ec_times_icemask_g - this is what introduces
       ! elevation class information from the glc grid (without elevation classes) to the 
       ! land grid (with elevation classes)
       ! Note that bare land values are mapped in the same way as ice-covered values
       ! Note that FBlnd_other does not contain fields with undistributed dimensions

       call med_map_FB_Regrid_Norm(&
            fldnames=fields_to_map(2:size(fields_to_map)), &
            FBSrc=FBglc_other, &
            FBDst=FBlnd_other, &
            FBfrac=FBglc_norm, mapnorm='Sg_frac_times_icemask', &
            RouteHandle=RouteHandle, &
            string='mapping other fields from glc to land elevation class '//trim(cvalue), rc=rc)

       ! Fill in FBlnd for other fields
       do nfld = 1,size(fields_to_map)
          call shr_nuopc_methods_FB_getFldPtr(FBlnd_other, trim(fields_to_map(nfld)), fldptr1=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_getFldPtr(FBlnd_elev, trim(fields_to_map(nfld))//'_elev', fldptr2=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          dataptr2d(ec+1,:) = dataptr1d(:)
       end do

       ! ------------------------------------------------------------------------
       ! set the topo field for virtual columns, in a given elevation class.
       ! ------------------------------------------------------------------------
       
       ! This is needed because virtual columns (i.e., elevation classes that have no
       ! contributing glc grid cells) won't have any topographic information mapped onto
       ! them, so would otherwise end up with an elevation of 0.
       
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_elev, trim(fields_to_map(ntopo_field))//'_elev', fldptr2=dataptr2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       topo_virtual = glc_mean_elevation_virtual(ec)
       do n = 1,size(dataptr2d, dim=2)
          if (glc_frac_l_ec(n) <= 0._r8) then
             dataptr2d(ec+1,n) = topo_virtual
          end if
       end do

    end do  ! loop over elevation classes

    ! clean up
    deallocate(glc_elevclass)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_map_glc2lnd_elevclass

  !================================================================================================

  subroutine get_glc_elevation_classes(glc_topo, glc_elevclass)

    !------------------
    ! Get elevation class of each grid cell on the glc grid.
    !
    ! This does not consider glc_frac: it simply gives the elevation class that the grid
    ! cell would be in if it were ice-covered. So it never returns an elevation class of
    ! 0 (bare land). (This design would allow us, in the future, to have glc grid cells
    ! that are part ice-covered, part ice-free.)
    !------------------

    use med_internalstate_mod , only : logunit
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
    use glc_elevclass_mod     , only : glc_get_elevation_class
    use shr_sys_mod           , only : shr_sys_abort

    ! input/output variables
    real(r8), intent(in)  :: glc_topo(:)      ! topographic height
    integer , intent(out) :: glc_elevclass(:) ! elevation class
    !
    ! local variables
    integer :: npts
    integer :: glc_pt
    integer :: err_code
    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

    do glc_pt = 1, npts
       call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)
       select case (err_code)
       case (GLC_ELEVCLASS_ERR_NONE)
          ! Do nothing
       case (GLC_ELEVCLASS_ERR_TOO_LOW, GLC_ELEVCLASS_ERR_TOO_HIGH)
          write(logunit,*) subname, ': WARNING, for glc_pt, topo = ', glc_pt, glc_topo(glc_pt)
          write(logunit,*) glc_errcode_to_string(err_code)
       case default
          write(logunit,*) subname, ': ERROR getting elevation class for glc_pt = ', glc_pt
          write(logunit,*) glc_errcode_to_string(err_code)
          call shr_sys_abort(subname//': ERROR getting elevation class')
       end select
    end do

  end subroutine get_glc_elevation_classes

end module med_map_glc2lnd_mod
