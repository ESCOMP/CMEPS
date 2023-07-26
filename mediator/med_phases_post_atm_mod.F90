module med_phases_post_atm_mod

  !-----------------------------------------------------------------------------
  ! Mediator phase for post atm calculations, maps atm->ice, atm->lnd, atm->ocn
  ! and atm->wav
  !-----------------------------------------------------------------------------

  implicit none
  private

  public :: med_phases_post_atm

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

   real(ESMF_KIND_R8) function field_integral(field, frac)
      use ESMF 

      ! inputs
      type(ESMF_Field)            :: field
      real(ESMF_KIND_R8), pointer :: frac(:)

      ! local vars
      type(ESMF_Mesh)             :: mesh
      type(ESMF_Grid)             :: grid
      real(ESMF_KIND_R8), pointer :: area1(:), fld_ptr1(:), area2(:, :), fld_ptr2(:, :)
      integer(ESMF_KIND_I4), pointer :: mask1(:), mask2(:, :)
      integer                     :: n, elementCount, i, j, lbnd1, ubnd1, lbnd2, ubnd2
      real(ESMF_KIND_R8)          :: local_integral(1), mask_integral(1), norm
      type(ESMF_VM)               :: vm  
      type(ESMF_GeomType_Flag)    :: geomtype
      logical                     :: is_masked
      character(len=200)               :: tmpstr
   

      local_integral(1) = 0.0
      mask_integral(1) = 0.0

      call ESMF_FieldGet(field, geomtype=geomtype)

      if (geomtype == ESMF_GEOMTYPE_MESH) then
         call ESMF_FieldGet(field, farrayptr=fld_ptr1)
         call ESMF_FieldGet(field, mesh=mesh)
         call ESMF_MeshGet(mesh, elementCount=elementCount)
         
         allocate(area1(elementCount))
         allocate(mask1(elementCount))
         
         call ESMF_MeshGet(mesh, elementArea=area1)

         call ESMF_MeshGet(mesh, elementMaskIsPresent=is_masked)

         if (is_masked) then
            call ESMF_MeshGet(mesh, elementMask=mask1)
         else
            mask1(:) = 1
         end if

         do n=1, elementCount
            local_integral(1) = local_integral(1) + area1(n) * fld_ptr1(n) * mask1(n) * frac(n)
            mask_integral(1) = mask_integral(1) + area1(n) * mask1(n) * frac(n)
         end do
      else

         call ESMF_FieldGet(field, farrayptr=fld_ptr2)
         call ESMF_FieldGet(field, grid=grid)

         call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=area2)

         call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=mask2)

         lbnd1 = lbound(area2, 1)
         ubnd1 = ubound(area2, 1)
         lbnd2 = lbound(area2, 2)
         ubnd2 = ubound(area2, 2)

         do j = lbnd2, ubnd2
            do i = lbnd1, ubnd1
               local_integral(1) = local_integral(1) + area2(i, j) * fld_ptr2(i, j) * mask2(i, j)
               mask_integral(1) = mask_integral(1) + area2(i, j) * mask2(i, j)
            end do
         end do

      end if

      field_integral = -1.0
      norm = -1.0
      call ESMF_VMGetCurrent(vm)
      call ESMF_VMAllFullReduce(vm, local_integral, field_integral, 1, ESMF_REDUCE_SUM)
      call ESMF_VMAllFullReduce(vm, mask_integral, norm, 1, ESMF_REDUCE_SUM)

      field_integral = field_integral / norm

      return

   end function field_integral

  subroutine med_phases_post_atm(gcomp, rc)

    !---------------------------------------
    ! map atm to ocn and atm to ice and atm to land
    !---------------------------------------
    use ESMF
    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated, ESMF_Field
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_phases_history_mod, only : med_phases_history_write_comp
    use med_map_mod           , only : med_map_field_packed
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr    => med_utils_ChkErr
    use med_internalstate_mod , only : compocn, compatm, compice, complnd, compwav
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*), parameter :: subname='(med_phases_post_atm)'
    
    integer                     :: fieldCount, i
    character(len=CL), pointer  :: fieldNameList(:)
    character(len=200)          :: tmpString
    character(len=1000)         :: msgString
    type(ESMF_Field), pointer   :: atmFieldList(:), ocnFieldList(:)
    real                        :: integral_1, integral_2, error
    real(ESMF_KIND_R8), pointer :: frac_1(:), frac_2(:), fld_ptr(:)
    type(ESMF_Field)            :: lfield
    
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map atm to ocn
    if (is_local%wrap%med_coupling_active(compatm,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compocn,:), &
            packed_data=is_local%wrap%packed_data(compatm,compocn,:), &
            routehandles=is_local%wrap%RH(compatm,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ocn')
       
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fieldCount=fieldCount)
       allocate(fieldNameList(fieldcount))
       allocate(atmFieldList(fieldcount))
       allocate(ocnFieldList(fieldcount))
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fieldNameList=fieldNameList, fieldList=atmFieldList)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compocn), fieldList=ocnFieldList)
       
       ! get fractions
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ofrac', field=lfield, rc=rc)
       call ESMF_FieldGet(lfield, farrayPtr=fld_ptr, rc=rc)
       allocate(frac_1(size(fld_ptr)))
       
       do i = 1,size(fld_ptr)
         frac_1(i) = fld_ptr(i)
       end do
       
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), fieldName='ifrac', field=lfield, rc=rc)
       call ESMF_FieldGet(lfield, farrayPtr=fld_ptr, rc=rc)
       do i = 1,size(fld_ptr)
         frac_1(i) = frac_1(i) + fld_ptr(i)
       end do

       call ESMF_FieldGet(ocnFieldList(1), farrayPtr=fld_ptr, rc=rc)
       allocate(frac_2(size(fld_ptr)))
       frac_2(:) = 1.0

       do i=1,fieldCount
         integral_1 = field_integral(atmFieldList(i), frac_1)
         integral_2 = field_integral(ocnFieldList(i), frac_2)
         error = (integral_2 - integral_1) / abs(integral_1)

         write (tmpString, *) error
         msgString = 'Field: ' // trim(fieldNameList(i)) // '. Error: ' //  trim(tmpString)
         write (tmpString, *) integral_1
         msgString = trim(msgString) // '. Atm grid integral: ' // trim(tmpString)
         write (tmpString, *) integral_2
         msgString = trim(msgString) // '. Ocn grid integral: ' // trim(tmpString) // '.'

         if (abs(integral_2) > 0) then
            call ESMF_Logwrite(trim(msgString))
         end if
       end do
       
       deallocate(fieldNameList)
       deallocate(atmFieldList)
       deallocate(ocnFieldList)
       deallocate(frac_1)
       deallocate(frac_2)

    end if
    ! map atm->ice
    if (is_local%wrap%med_coupling_active(compatm,compice)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ice')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compice), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compice,:), &
            packed_data=is_local%wrap%packed_data(compatm,compice,:), &
            routehandles=is_local%wrap%RH(compatm,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ice')
    end if
    ! map atm->lnd
    if (is_local%wrap%med_coupling_active(compatm,complnd)) then
       call t_startf('MED:'//trim(subname)//' map_atm2lnd')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,complnd), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,complnd,:), &
            packed_data=is_local%wrap%packed_data(compatm,complnd,:), &
            routehandles=is_local%wrap%RH(compatm,complnd,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2lnd')
    end if
    ! map atm->wav
    if (is_local%wrap%med_coupling_active(compatm,compwav)) then
       call t_startf('MED:'//trim(subname)//' map_atm2wav')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compwav), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compwav,:), &
            packed_data=is_local%wrap%packed_data(compatm,compwav,:), &
            routehandles=is_local%wrap%RH(compatm,compwav,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2wav')
    end if

    ! Write atm inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_comp(gcomp, compatm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_atm

end module med_phases_post_atm_mod
