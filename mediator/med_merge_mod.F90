module med_merge_mod

  !-----------------------------------------------------------------------------
  ! Performs merges from source field bundles to destination field bundle
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : logunit, compmed, compname
  use med_constants_mod     , only : dbug_flag         => med_constants_dbug_flag
  use med_constants_mod     , only : czero             => med_constants_czero
  use med_utils_mod         , only : ChkErr            => med_utils_ChkErr
  use med_methods_mod       , only : FB_FldChk         => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr      => med_methods_FB_GetFldPtr
  use esmFlds               , only : med_fldList_type
  use esmFlds               , only : med_fld_GetFldInfo
  use esmFlds               , only : med_fldList_entry_type
  use esmFlds               , only : med_fldList_findName
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_merge_auto
  public  :: med_merge_field

  interface med_merge_auto ; module procedure &
       med_merge_auto_single_fldbun, &
       med_merge_auto_multi_fldbuns
  end interface

  interface med_merge_field ; module procedure &
       med_merge_field_1D
  end interface

  private :: med_merge_auto_field

  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_merge_auto_multi_fldbuns(coupling_active, FBOut, FBfrac, FBImp, fldListTo, FBMed1, FBMed2, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_LogSetError

    ! ----------------------------------------------
    ! Auto merge based on fldListTo info
    ! ----------------------------------------------

    ! input/output variables
    logical                , intent(in)            :: coupling_active(:) ! true => coupling is active
    type(ESMF_FieldBundle) , intent(inout)         :: FBOut              ! Merged output field bundle
    type(ESMF_FieldBundle) , intent(inout)         :: FBfrac             ! Fraction data for FBOut
    type(ESMF_FieldBundle) , intent(in)            :: FBImp(:)          ! Array of field bundles each mapping to the FBOut mesh
    type(med_fldList_type) , intent(in) , target   :: fldListTo          ! Information for merging
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed1             ! mediator field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed2             ! mediator field bundle
    integer                , intent(out)           :: rc

    ! local variables
    type(med_fldList_entry_type), pointer :: fldptr
    integer                    :: nfld_out,nfld_in,nm
    integer                    :: compsrc
    integer                    :: num_merge_fields
    integer                    :: num_merge_colon_fields
    character(CL)              :: merge_fields
    character(CL)              :: merge_field
    character(CS)              :: merge_type
    character(CS)              :: merge_fracname
    character(CS), pointer     :: merge_field_names(:)
    logical                    :: error_check = .false.  ! TODO: make this an input argument
    integer                    :: ungriddedUBound_out(1) ! size of ungridded dimension
    integer                    :: fieldcount
    character(CL)   , pointer  :: fieldnamelist(:)
    type(ESMF_Field), pointer  :: fieldlist(:)
    real(r8), pointer          :: dataptr1d(:)
    real(r8), pointer          :: dataptr2d(:,:)
    logical                    :: zero_output
    logical                    :: found
    character(len=*),parameter :: subname=' (module_med_merge_mod: med_merge_auto)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    rc = ESMF_SUCCESS

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
       rc = ESMF_SUCCESS
    end if

    call ESMF_FieldBundleGet(FBOut, fieldCount=fieldcount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnamelist(fieldcount))
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(FBOut, fieldnamelist=fieldnamelist, fieldlist=fieldlist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Want to loop over all of the fields in FBout here - and find the corresponding index in fldListTo(compxxx)
    ! for that field name - then call the corresponding merge routine below appropriately

    ! Loop over all fields in field bundle FBOut
    do nfld_out = 1,fieldcount
       zero_output = .true.
       call ESMF_FieldGet(fieldlist(nfld_out), ungriddedUBound=ungriddedUbound_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Find the next fieldname
       call med_fldList_findName(fldListTo%fields, fieldnamelist(nfld_out), found, fldptr)
       if (found) then
          ! Loop over all possible source components in the merging arrays returned from the above call
          ! If the merge field name from the source components is not set, then simply go to the next component
          do compsrc = 1,size(FBImp)
             ! Cycle if coupling is not active or mediator input is not present and compsrc is mediator
             if (compsrc == compmed) then
                if (.not. present(FBMed1) .and. .not. present(FBMed2)) then
                   CYCLE
                end if
             else if (.not. coupling_active(compsrc)) then
                CYCLE
             end if
                
             ! Determine the merge information for the import field
             call med_fld_GetFldInfo(fldptr, compsrc=compsrc, merge_fields=merge_fields, merge_type=merge_type, merge_fracname=merge_fracname, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             if (merge_type /= 'unset' .and. merge_field /= 'unset') then
                   ! If merge_field is a colon delimited string then cycle through every field - otherwise by default nm
                   ! will only equal 1
                num_merge_colon_fields = merge_listGetNum(merge_fields)
                do nm = 1,num_merge_colon_fields
                   ! Determine merge field name from source field
                   call merge_listGetName(merge_fields, nm, merge_field, rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   ! Perform error checks
                   if (error_check) then
                      call med_merge_auto_errcheck(compsrc, fieldnamelist(nfld_out), fieldlist(nfld_out), &
                           ungriddedUBound_out, trim(merge_field), FBImp(compsrc), &
                           FBMed1=FBMed1, FBMed2=FBMed2, rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   end if ! end of error check

                   ! Initialize initial output field data to zero before doing merge
                   if (zero_output) then
                      if (ungriddedUBound_out(1) > 0) then
                         call ESMF_FieldGet(fieldlist(nfld_out), farrayPtr=dataptr2d, rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                         dataptr2d(:,:) = czero
                      else
                         call ESMF_FieldGet(fieldlist(nfld_out), farrayPtr=dataptr1d, rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                         dataptr1d(:) = czero
                      end if
                      zero_output = .false.
                   end if

                   ! Perform merge
                   if ((present(FBMed1) .or. present(FBMed2)) .and. compsrc == compmed) then
                      if (FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                         call med_merge_auto_field(trim(merge_type), fieldlist(nfld_out), ungriddedUBound_out, &
                              FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      else if (FB_FldChk(FBMed2, trim(merge_field), rc=rc)) then
                         call med_merge_auto_field(trim(merge_type), fieldlist(nfld_out), ungriddedUBound_out, &
                              FB=FBMed2, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if
                   else
                      call med_merge_auto_field(trim(merge_type), fieldlist(nfld_out), ungriddedUBound_out, &
                           FB=FBImp(compsrc), FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   end if

                end do ! end of nm loop
             end if ! end of check of merge_type and merge_field not unset
          end do  ! end of compsrc loop
       end if ! end if found
    end do ! end of loop over fields in FBOut

    deallocate(fieldnamelist)
    deallocate(fieldlist)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_merge_auto_multi_fldbuns

  !===============================================================================
  subroutine med_merge_auto_single_fldbun(compsrc, FBOut, FBfrac, FBIn, fldListTo, rc)

    ! ----------------------------------------------
    ! Auto merge from one import field bundle based on fldListTo info.
    ! Want to loop over all of the fields in FBout here - and find the
    ! corresponding index in fldListTo for that field name - then call
    ! the corresponding merge routine below appropriately.
    ! ----------------------------------------------

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF , only : ESMF_LogSetError
    ! input/output variables
    integer                , intent(in)    :: compsrc
    type(ESMF_FieldBundle) , intent(inout) :: FBOut     ! Merged output field bundle
    type(ESMF_FieldBundle) , intent(inout) :: FBfrac    ! Fraction data for FBOut
    type(ESMF_FieldBundle) , intent(in)    :: FBIn      ! Single field bundle to merge to the FBOut mesh
    type(med_fldList_type) , intent(in), target :: fldListTo ! Information for merging
    integer                , intent(out)   :: rc

    ! local variables
    type(med_fldList_entry_type), pointer :: fldptr
    integer                    :: nfld_out,nfld_in,nm
    integer                    :: num_merge_fields
    integer                    :: num_merge_colon_fields
    character(CL)              :: merge_fields
    character(CL)              :: merge_field
    character(CS)              :: merge_type
    character(CS)              :: merge_fracname
    character(CS)              :: merge_field_name
    integer                    :: ungriddedUBound_out(1) ! size of ungridded dimension
    integer                    :: fieldcount
    character(CL)   , pointer  :: fieldnamelist(:)
    type(ESMF_Field), pointer  :: fieldlist(:)
    real(r8), pointer          :: dataptr1d(:)
    real(r8), pointer          :: dataptr2d(:,:)
    logical                    :: zero_output
    logical                    :: found
    character(len=*),parameter :: subname=' (module_med_merge_mod: med_merge_auto)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    rc = ESMF_SUCCESS

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
       rc = ESMF_SUCCESS
    end if

    call ESMF_FieldBundleGet(FBOut, fieldCount=fieldcount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnamelist(fieldcount))
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(FBOut, fieldnamelist=fieldnamelist, fieldlist=fieldlist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Loop over all fields in output field bundle FBOut
    do nfld_out = 1,fieldcount
       zero_output = .true.
       call ESMF_FieldGet(fieldlist(nfld_out), ungriddedUBound=ungriddedUbound_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Find the next fieldname
       call med_fldList_findName(fldListTo%fields, fieldnamelist(nfld_out), found, fldptr)
       if(found) then
          ! Determine the merge information for the import field
          call med_fld_GetFldInfo(fldptr, compsrc=compsrc, merge_fields=merge_fields, merge_type=merge_type, merge_fracname=merge_fracname)
          if (merge_type /= 'unset' .and. merge_fields /= 'unset') then

             ! If merge_field is a colon delimited string then cycle through every field - otherwise by default nm
             ! will only equal 1
             num_merge_colon_fields = merge_listGetNum(merge_fields)
             do nm = 1,num_merge_colon_fields
                ! Determine merge field name from source field
                call merge_listGetName(merge_fields, nm, merge_field, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                
                ! Initialize initial output field data to zero before doing merge
                if (zero_output) then
                   if (ungriddedUBound_out(1) > 0) then
                      call ESMF_FieldGet(fieldlist(nfld_out), farrayPtr=dataptr2d, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      dataptr2d(:,:) = czero
                   else
                      call ESMF_FieldGet(fieldlist(nfld_out), farrayPtr=dataptr1d, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                      dataptr1d(:) = czero
                   end if
                   zero_output = .false.
                end if
                
                ! Perform merge
                call med_merge_auto_field(trim(merge_type), fieldlist(nfld_out), ungriddedUBound_out, &
                     FB=FBIn, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                
             end do ! end of nm loop
          end if ! end of check of merge_type and merge_field not unset
       end if ! end of check if stdname and fldname are the same
    end do ! end of loop over fields in FBOut

    deallocate(fieldnamelist)
    deallocate(fieldlist)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_merge_auto_single_fldbun

  !===============================================================================
  subroutine med_merge_auto_field(merge_type, field_out, ungriddedUBound_out,  &
       FB, FBfld, FBw, fldw, rc)

    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogMsg_Error
    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF , only : ESMF_LogWrite, ESMF_LogMsg_Info
    use ESMF , only : ESMF_FieldGet, ESMF_Field

    ! input/output variables
    character(len=*)      ,intent(in)    :: merge_type
    type(ESMF_Field)      ,intent(inout) :: field_out
    integer               ,intent(in)    :: ungriddedUBound_out(1)
    type(ESMF_FieldBundle),intent(in)    :: FB
    character(len=*)      ,intent(in)    :: FBfld
    type(ESMF_FieldBundle),intent(inout) :: FBw     ! field bundle with weights
    character(len=*)      ,intent(in)    :: fldw    ! name of weight field to use in FBw
    integer               ,intent(out)   :: rc

    ! local variables
    integer           :: n
    type(ESMF_Field)  :: field_wgt
    type(ESMF_Field)  :: field_in
    real(R8), pointer :: dp1 (:)
    real(R8), pointer :: dp2(:,:)   ! output pointers to 1d and 2d fields
    real(R8), pointer :: dpf1(:)
    real(R8), pointer :: dpf2(:,:)  ! intput pointers to 1d and 2d fields
    real(R8), pointer :: dpw1(:)    ! weight pointer
    character(CL) :: name
    character(len=*),parameter :: subname=' (med_merge_mod: med_merge_auto_field)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    !-------------------------
    ! Error checks
    !-------------------------

    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge') then
       if (trim(fldw) == 'unset') then
          call ESMF_LogWrite(trim(subname)//": error required merge_fracname is not set", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       end if
       if (.not. FB_FldChk(FBw, trim(fldw), rc=rc)) then
          call ESMF_LogWrite(trim(subname)//": error "//trim(fldw)//"is not in FBw", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       end if
    end if

    !-------------------------
    ! Get appropriate field pointers
    !-------------------------

    ! Get input field
    call ESMF_FieldBundleGet(FB, FBfld, field=field_in, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get field pointer to output and input fields
    ! Assume that input and output ungridded upper bounds are the same - this is checked in error check

    if (ungriddedUBound_out(1) > 0) then
       call ESMF_FieldGet(field_in, farrayPtr=dpf2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_out, farrayPtr=dp2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_FieldGet(field_in, farrayPtr=dpf1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_out, farrayPtr=dp1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Get pointer to weights that weights have no ungridded dimensions
    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
       call ESMF_FieldBundleGet(FBw, fieldName=trim(fldw), field=field_wgt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_wgt, farrayPtr=dpw1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    ! Do supported merges
    if (trim(merge_type)  == 'copy') then
       if (ungriddedUBound_out(1) > 0) then
          dp2(:,:) = dpf2(:,:)
       else
          dp1(:) = dpf1(:)
       endif
    else if (trim(merge_type)  == 'copy_with_weights') then
       if (ungriddedUBound_out(1) > 0) then
          do n = 1,ungriddedUBound_out(1)
             dp2(n,:) = dpf2(n,:)*dpw1(:)
          end do
       else
          dp1(:) = dpf1(:)*dpw1(:)
       endif
    else if (trim(merge_type)  == 'merge' .or. trim(merge_type) == 'sum_with_weights') then
       if (ungriddedUBound_out(1) > 0) then
          do n = 1,ungriddedUBound_out(1)
             dp2(n,:) = dp2(n,:) + dpf2(n,:)*dpw1(:)
          end do
       else
          dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
       endif
    else if (trim(merge_type) == 'sum') then
       if (ungriddedUBound_out(1) > 0) then
          dp2(:,:) = dp2(:,:) + dpf2(:,:)
       else
          dp1(:) = dp1(:) + dpf1(:)
       endif
    else
       call ESMF_LogWrite(trim(subname)//": merge type "//trim(merge_type)//" not supported", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine med_merge_auto_field

  !===============================================================================
  subroutine med_merge_auto_errcheck(compsrc, fldname_out, field_out, &
       ungriddedUBound_out, merge_fldname, FBImp, FBMed1, FBMed2, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF , only : ESMF_LogSetError, ESMF_RC_OBJ_NOT_CREATED

    ! input/output variables
    integer                , intent(in)            :: compsrc                ! source component index
    character(len=*)       , intent(in)            :: fldname_out            ! output field name
    type(ESMF_Field)       , intent(in)            :: field_out              ! output field
    integer                , intent(in)            :: ungriddedUBound_out(1) ! ungridded upper bound
    character(len=*)       , intent(in)            :: merge_fldname          ! source  merge fieldname
    type(ESMF_FieldBundle) , intent(in)            :: FBImp                  ! source field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed1                 ! mediator field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed2                 ! mediator field bundle
    integer                , intent(out)           :: rc

    ! local variables
    type(ESMF_Field)  :: field_in
    integer           :: ungriddedUBound_in(1)  ! size of ungridded dimension, if any
    character(len=CL) :: errmsg
    character(len=*),parameter :: subname=' (module_med_merge_mod: med_merge_errcheck)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    if (compsrc == compmed) then
       if (present(FBMed1) .and. present(FBMed2)) then
          if (.not. ESMF_FieldBundleIsCreated(FBMed1)) then
             call ESMF_LogSetError(ESMF_RC_OBJ_NOT_CREATED, msg="Field bundle FBMed1 not created.", &
                  line=__LINE__, file=u_FILE_u, rcToReturn=rc)
             return
          endif
          if (.not. ESMF_FieldBundleIsCreated(FBMed2)) then
             call ESMF_LogSetError(ESMF_RC_OBJ_NOT_CREATED, msg="Field bundle FBMed2 not created.", &
                  line=__LINE__, file=u_FILE_u, rcToReturn=rc)
             return
          endif
          if (FB_FldChk(FBMed1, trim(merge_fldname), rc=rc)) then
             call ESMF_FieldBundleGet(FBMed1, trim(merge_fldname), field=field_in, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else if (FB_FldChk(FBMed2, trim(merge_fldname), rc=rc)) then
             call ESMF_FieldBundleGet(FBMed2, trim(merge_fldname), field=field_in, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//": ERROR merge_fldname = "//trim(merge_fldname)//" not found", &
                  ESMF_LOGMSG_ERROR, rc=rc)
             rc = ESMF_FAILURE
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    endif

    call ESMF_FieldBundleGet(FBImp, trim(merge_fldname), field=field_in, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_in, ungriddedUBound=ungriddedUBound_in, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (ungriddedUbound_in(1) /= UngriddedUbound_out(1)) then
       write(errmsg,*) trim(subname),' input field ungriddedUbound ',ungriddedUbound_in(1),&
            ' for '//trim(merge_fldname), &
            ' not equal to output field ungriddedUbound ',ungriddedUbound_out,' for '//trim(fldname_out)
       call ESMF_LogWrite(errmsg, ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
   endif

  end subroutine med_merge_auto_errcheck

  !===============================================================================
  subroutine med_merge_field_1D(FBout, fnameout, &
                                FBinA, fnameA, wgtA, &
                                FBinB, fnameB, wgtB, &
                                FBinC, fnameC, wgtC, &
                                FBinD, fnameD, wgtD, &
                                FBinE, fnameE, wgtE, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_LogWrite
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_ERROR
    use ESMF , only : ESMF_LOGMSG_WARNING, ESMF_LOGMSG_INFO

    ! ----------------------------------------------
    ! Supports up to a five way merge
    ! ----------------------------------------------

    ! input/output variabes
    type(ESMF_FieldBundle) , intent(inout)                 :: FBout
    character(len=*)       , intent(in)                    :: fnameout
    type(ESMF_FieldBundle) , intent(in)                    :: FBinA
    character(len=*)       , intent(in)                    :: fnameA
    real(R8)               , intent(in), pointer           :: wgtA(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinB
    character(len=*)       , intent(in), optional          :: fnameB
    real(R8)               , intent(in), optional, pointer :: wgtB(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinC
    character(len=*)       , intent(in), optional          :: fnameC
    real(R8)               , intent(in), optional, pointer :: wgtC(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinD
    character(len=*)       , intent(in), optional          :: fnameD
    real(R8)               , intent(in), optional, pointer :: wgtD(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinE
    character(len=*)       , intent(in), optional          :: fnameE
    real(R8)               , intent(in), optional, pointer :: wgtE(:)
    integer                , intent(out)                   :: rc

    ! local variables
    real(R8), pointer          :: dataOut(:)
    real(R8), pointer          :: dataPtr(:)
    real(R8), pointer          :: wgt(:)
    integer                    :: lb1,ub1,i,j,n
    logical                    :: wgtfound, FBinfound
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_merge_field_1D)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc=ESMF_SUCCESS

    ! check each field has a fieldname passed in
    if ((present(FBinB) .and. .not.present(fnameB)) .or. &
        (present(FBinC) .and. .not.present(fnameC)) .or. &
        (present(FBinD) .and. .not.present(fnameD)) .or. &
        (present(FBinE) .and. .not.present(fnameE))) then

       call ESMF_LogWrite(trim(subname)//": ERROR fname not present with FBin", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    endif

    if (.not. FB_FldChk(FBout, trim(fnameout), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    call FB_GetFldPtr(FBout, trim(fnameout), fldptr1=dataOut, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)

    dataOut = czero

    ! check that each field passed in actually exists, if not DO NOT do any merge
    FBinfound = .true.
    if (present(FBinB)) then
       if (.not. FB_FldChk(FBinB, trim(fnameB), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinC)) then
       if (.not. FB_FldChk(FBinC, trim(fnameC), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinD)) then
       if (.not. FB_FldChk(FBinD, trim(fnameD), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinE)) then
       if (.not. FB_FldChk(FBinE, trim(fnameE), rc=rc)) FBinfound = .false.
    endif
    if (.not. FBinfound) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not found in FBin, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    ! n=1,5 represents adding A to E inputs if they exist
    do n = 1,5
       FBinfound = .false.
       wgtfound = .false.

       if (n == 1) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinA, trim(fnameA), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          wgtfound = .true.
          wgt => wgtA

       elseif (n == 2 .and. present(FBinB)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinB, trim(fnameB), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtB)) then
             wgtfound = .true.
             wgt => wgtB
          endif

       elseif (n == 3 .and. present(FBinC)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinC, trim(fnameC), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtC)) then
             wgtfound = .true.
             wgt => wgtC
          endif

       elseif (n == 4 .and. present(FBinD)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinD, trim(fnameD), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtD)) then
             wgtfound = .true.
             wgt => wgtD
          endif

       elseif (n == 5 .and. present(FBinE)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinE, trim(fnameE), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtE)) then
             wgtfound = .true.
             wgt => wgtE
          endif

       endif

       if (FBinfound) then
          if (lbound(dataPtr,1) /= lbound(dataOut,1) .or. ubound(dataPtr,1) /= ubound(dataOut,1)) then
             call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif
          if (wgtfound) then
             if (lbound(dataPtr,1) /= lbound(wgt,1) .or. ubound(dataPtr,1) /= ubound(wgt,1)) then
                call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
             do i = lb1,ub1
                dataOut(i) = dataOut(i) + dataPtr(i) * wgt(i)
             enddo
          else
             do i = lb1,ub1
                dataOut(i) = dataOut(i) + dataPtr(i)
             enddo
          endif  ! wgtfound

       endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_merge_field_1D

  !===============================================================================
  integer function merge_listGetNum(str)

    !  return number of fields in a colon delimited string list

    ! input/output variables
    character(*),intent(in) :: str   ! string to search

    ! local variables
    integer          :: n
    integer          :: count          ! counts occurances of char
    character(len=1) :: listDel  = ":" ! note single exec implications
    !---------------------------------------

    merge_listGetNum = 0
    if (len_trim(str) > 0) then
       count = 0
       do n = 1, len_trim(str)
          if (str(n:n) == listDel) count = count + 1
       end do
       merge_listGetNum = count + 1
    endif

  end function merge_listGetNum

  !===============================================================================
  subroutine merge_listGetName(list, k, name, rc)

    ! Get name of k-th field in colon deliminted list

    use ESMF, only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite, ESMF_LOGMSG_INFO

    ! input/output variables
    character(len=*)  ,intent(in)  :: list    ! list/string
    integer           ,intent(in)  :: k       ! index of field
    character(len=*)  ,intent(out) :: name    ! k-th name in list
    integer, optional ,intent(out) :: rc      ! return code

    ! local variables
    integer          :: i,n   ! generic indecies
    integer          :: kFlds ! number of fields in list
    integer          :: i0,i1 ! name = list(i0:i1)
    integer          :: nChar
    logical          :: valid_list
    character(len=1) :: listDel  = ':'
    character(len=2) :: listDel2 = '::'
    !---------------------------------------

    rc = ESMF_SUCCESS
    if(k==1) then
       name = trim(list)
       return
    endif
    ! check that this is a valid list
    valid_list = .true.
    nChar = len_trim(list)
    if (nChar < 1) then                           ! list is an empty string
       valid_list = .false.
    else if (    list(1:1)     == listDel  ) then ! first char is delimiter
       valid_list = .false.
    else if (list(nChar:nChar) == listDel  ) then ! last  char is delimiter
       valid_list = .false.
    else if (index(trim(list)," " )     > 0) then ! white-space in a field name
       valid_list = .false.
    else if (index(trim(list),listDel2) > 0) then ! found zero length field
       valid_list = .false.
    end if
    if (.not. valid_list) then
       write(logunit,*) "ERROR: invalid list = ",trim(list)
       call ESMF_LogWrite("ERROR: invalid list = "//trim(list), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    !--- check that this is a valid index ---
    kFlds = merge_listGetNum(list)
    if (k<1 .or. kFlds<k) then
       write(logunit,*) "ERROR: invalid index = ",k
       write(logunit,*) "ERROR:          list = ",trim(list)
       call ESMF_LogWrite("ERROR: invalid index = "//trim(list), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    ! start with whole list, then remove fields before and after desired field ---
    i0 = 1
    i1 = len_trim(list)

    ! remove field names before desired field
    do n=2,k
       i = index(list(i0:i1),listDel)
       i0 = i0 + i
    end do

    ! remove field names after desired field
    if ( k < kFlds ) then
       i = index(list(i0:i1),listDel)
       i1 = i0 + i - 2
    end if

    ! copy result into output variable
    name = list(i0:i1)//"   "

  end subroutine merge_listGetName

end module med_merge_mod
