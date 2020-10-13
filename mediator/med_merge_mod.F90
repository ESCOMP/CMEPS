module med_merge_mod

  !-----------------------------------------------------------------------------
  ! Performs merges from source field bundles to destination field bundle
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : logunit
  use med_constants_mod     , only : dbug_flag         => med_constants_dbug_flag
  use med_constants_mod     , only : spval_init        => med_constants_spval_init
  use med_constants_mod     , only : spval             => med_constants_spval
  use med_constants_mod     , only : czero             => med_constants_czero
  use med_utils_mod         , only : ChkErr            => med_utils_ChkErr
  use med_methods_mod       , only : FB_FldChk         => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr      => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FieldPtr_Compare  => med_methods_FieldPtr_Compare
  use esmFlds               , only : compmed, compname
  use esmFlds               , only : med_fldList_type
  use esmFlds               , only : med_fldList_GetNumFlds
  use esmFlds               , only : med_fldList_GetFldInfo
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_merge_auto
  public  :: med_merge_field

  interface med_merge_field ; module procedure &
       med_merge_field_1D, &
       med_merge_field_2D
  end interface

  private :: med_merge_auto_field

  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================
  
  subroutine med_merge_auto(compout, coupling_active, FBOut, FBfrac, FBImp, fldListTo, &
       FBMed1, FBMed2, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF , only : ESMF_LogSetError, ESMF_RC_OBJ_NOT_CREATED

    ! ----------------------------------------------
    ! Auto merge based on fldListTo info
    ! ----------------------------------------------

    ! input/output variables
    integer                , intent(in)            :: compout            ! component index for FBOut
    logical                , intent(in)            :: coupling_active(:) ! true => coupling is active
    type(ESMF_FieldBundle) , intent(inout)         :: FBOut              ! Merged output field bundle
    type(ESMF_FieldBundle) , intent(inout)         :: FBfrac             ! Fraction data for FBOut
    type(ESMF_FieldBundle) , intent(in)            :: FBImp(:)           ! Array of field bundles each mapping to the FBOut mesh
    type(med_fldList_type) , intent(in)            :: fldListTo          ! Information for merging
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed1             ! mediator field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed2             ! mediator field bundle
    integer                , intent(out)           :: rc

    ! local variables
    integer                    :: nfld_out,nfld_in,nm
    integer                    :: compsrc
    integer                    :: num_merge_fields
    integer                    :: num_merge_colon_fields 
    character(CL)              :: fldname_out
    character(CL)              :: merge_fields
    character(CL)              :: merge_field
    character(CS)              :: merge_type
    character(CS)              :: merge_fracname
    character(CS), allocatable :: merge_field_names(:)
    integer                    :: rank_out
    type(ESMF_Field)           :: field_out
    real(r8), pointer          :: dataptr1d(:)
    real(r8), pointer          :: dataptr2d(:,:)
    logical                    :: error_check = .false. ! TODO: make this an input argument
    integer                    :: ungriddedUBound_out(1) ! currently the size must equal 1 for rank 2 fieldds
    integer                    :: gridToFieldMap_out(1)  ! currently the size must equal 1 for rank 2 fieldds
    integer                    :: fieldnamecount
    character(CL), pointer     :: fieldnamelist(:)
    character(CL)              :: msg
    character(len=*),parameter :: subname=' (module_med_merge_mod: med_merge_auto)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
       rc = ESMF_SUCCESS
    end if

    call ESMF_FieldBundleGet(FBOut, fieldCount=fieldnamecount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnamelist(fieldnamecount))
    call ESMF_FieldBundleGet(FBOut, fieldNameList=fieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    num_merge_fields = med_fldList_GetNumFlds(fldListTo)
    allocate(merge_field_names(num_merge_fields))
    do nfld_in = 1,num_merge_fields
       call med_fldList_GetFldInfo(fldListTo, nfld_in, merge_field_names(nfld_in))
    end do

    ! Want to loop over all of the fields in FBout here - and find the corresponding index in fldListTo(compxxx)
    ! for that field name - then call the corresponding merge routine below appropriately

    ! Loop over all fields in field bundle FBOut
    do nfld_out = 1,fieldnamecount

       ! Get the nth field in FBOut
       fldname_out = trim(fieldnamelist(nfld_out))
       call ESMF_FieldBundleGet(FBOut, trim(fldname_out), field=field_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Set nth field data to zero
       call ESMF_FieldGet(field_out, rank=rank_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          write(msg,*)trim(subname),'output field ',trim(fldname_out),' has rank ',rank_out
          call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
       end if

       if (rank_out == 1) then
          call ESMF_FieldGet(field_out, farrayPtr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d(:) = czero
       else if (rank_out == 2) then
          call ESMF_FieldGet(field_out, farrayptr=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr2d(:,:) = czero
          call ESMF_FieldGet(field_out, ungriddedUBound=ungriddedUBound_out, &
               gridToFieldMap=gridToFieldMap_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Loop over the field in fldListTo
       do nfld_in = 1,med_fldList_GetNumFlds(fldListTo)

          if (trim(merge_field_names(nfld_in)) == trim(fldname_out)) then

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
                call med_fldList_GetFldInfo(fldListTo, nfld_in, compsrc, merge_fields, merge_type, merge_fracname)

                if (merge_type /= 'unset' .and. merge_field /= 'unset') then

                   ! If merge_field is a colon delimited string then cycle through every field - otherwise by default nm
                   ! will only equal 1
                   num_merge_colon_fields = merge_listGetNum(merge_fields)
                   do nm = 1,num_merge_colon_fields

                      ! Determine merge field name from source field
                      if (num_merge_fields == 1) then
                         merge_field = trim(merge_fields)
                      else
                         call merge_listGetName(merge_fields, nm, merge_field, rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if

                      ! Perform error checks  
                      if (error_check) then
                         call med_merge_errcheck(compsrc, fldname_out, field_out, rank_out, &
                              ungriddedUBound_out, gridToFieldMap_out, trim(merge_field), &
                              FBImp(compsrc), FBMed1=FBMed1, FBMed2=FBMed2, rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if ! end of error check

                      ! Perform merge
                      if ((present(FBMed1) .or. present(FBMed2)) .and. compsrc == compmed) then
                         if (FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                            call med_merge_auto_field(trim(merge_type), field_out, &
                                 rank_out, ungriddedUBound_out, gridToFieldMap_out, &
                                 FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                            if (ChkErr(rc,__LINE__,u_FILE_u)) return
                         else if (FB_FldChk(FBMed2, trim(merge_field), rc=rc)) then
                            call med_merge_auto_field(trim(merge_type), field_out, &
                                 rank_out, ungriddedUBound_out, gridToFieldMap_out, &
                                 FB=FBMed2, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                            if (ChkErr(rc,__LINE__,u_FILE_u)) return
                         end if
                      else
                         call med_merge_auto_field(trim(merge_type), field_out, &
                              rank_out, ungriddedUBound_out, gridToFieldMap_out, &
                              FB=FBImp(compsrc), FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if

                   end do ! end of nm loop
                end if ! end of check of merge_type and merge_field not unset
             end do  ! end of compsrc loop
          end if ! end of check if stdname and fldname are the same
       end do ! end of loop over fldsListTo
    end do ! end of loop over fields in FBOut

    deallocate(fieldnamelist)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_merge_auto

  !===============================================================================
  subroutine med_merge_errcheck(compsrc, fldname_out, field_out, rank_out, &
       ungriddedUBound_out, gridToFieldMap_out, merge_fldname, FBImp, FBMed1, FBMed2, rc)

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF , only : ESMF_LogSetError, ESMF_RC_OBJ_NOT_CREATED

    ! input/output variables
    integer                , intent(in)            :: compsrc                ! source component index 
    character(len=*)       , intent(in)            :: fldname_out            ! output field name 
    type(ESMF_Field)       , intent(in)            :: field_out              ! output field
    integer                , intent(in)            :: rank_out               ! rank of output field
    integer                , intent(in)            :: ungriddedUBound_out(1) ! currently the size must equal 1 for rank 2 fields
    integer                , intent(in)            :: gridToFieldMap_out(1)  ! currently the size must equal 1 for rank 2 fields
    character(len=*)       , intent(in)            :: merge_fldname          ! source  merge fieldname 
    type(ESMF_FieldBundle) , intent(in)            :: FBImp                  ! source field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed1                 ! mediator field bundle
    type(ESMF_FieldBundle) , intent(in) , optional :: FBMed2                 ! mediator field bundle
    integer                , intent(out)           :: rc

    ! local variables
    type(ESMF_Field)  :: field_in
    integer           :: rank_in
    integer           :: ungriddedUBound_in(1)  ! currently the size must equal 1 for rank 2 fieldds
    integer           :: gridToFieldMap_in(1)   ! currently the size must equal 1 for rank 2 fieldds
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
    call ESMF_FieldGet(field_in, rank=rank_in, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (rank_in /= rank_out) then
       write(errmsg,*) trim(subname),' input field rank ',rank_in,' for '//trim(merge_fldname), &
            ' not equal to output field rank ',rank_out,' for '//trim(fldname_out)
       call ESMF_LogWrite(errmsg, ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    else if (rank_out == 2) then
       call ESMF_FieldGet(field_in, ungriddedUBound=ungriddedUBound_in, &
            gridToFieldMap=gridToFieldMap_in, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (ungriddedUBound_out(1) /= ungriddedUBound_in(1)) then
          write(errmsg,*) trim(subname),"ungriddedUBound_in (",ungriddedUBound_in(1),&
               ") not equal to ungriddedUBound_out (",ungriddedUBound_out(1),") for "//trim(fldname_out)
          call ESMF_LogWrite(errmsg, ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       else if (gridToFieldMap_in(1) /= gridToFieldMap_out(1)) then
          write(errmsg,*) trim(subname),"gridtofieldmap_in (",gridtofieldmap_in(1),&
               ") not equal to gridtofieldmap_out (",gridtofieldmap_out(1),") for "//trim(fldname_out)
          call ESMF_LogWrite(errmsg, ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
    endif

  end subroutine med_merge_errcheck

  !===============================================================================
  subroutine med_merge_auto_field(merge_type, field_out, &
       rank_out, ungriddedUBound_out, gridToFieldMap_out, &
       FB, FBfld, FBw, fldw, rc)

    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogMsg_Error
    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF , only : ESMF_LogWrite, ESMF_LogMsg_Info
    use ESMF , only : ESMF_FieldGet, ESMF_Field

    ! input/output variables
    integer               ,intent(in)    :: rank_out  ! rank of output array
    character(len=*)      ,intent(in)    :: merge_type
    type(ESMF_Field)      ,intent(inout) :: field_out
    integer               ,intent(in)    :: ungriddedUBound_out(1)
    integer               ,intent(in)    :: gridToFieldMap_out(1)
    type(ESMF_FieldBundle),intent(in)    :: FB
    character(len=*)      ,intent(in)    :: FBfld
    type(ESMF_FieldBundle),intent(inout) :: FBw     ! field bundle with weights
    character(len=*)      ,intent(in)    :: fldw    ! name of weight field to use in FBw
    integer               ,intent(out)   :: rc

    ! local variables
    integer           :: n
    type(ESMF_Field)  :: field_wgt
    type(ESMF_Field)  :: field_in
    real(R8), pointer :: dp1 (:), dp2(:,:)         ! output pointers to 1d and 2d fields
    real(R8), pointer :: dpf1(:), dpf2(:,:)        ! intput pointers to 1d and 2d fields
    real(R8), pointer :: dpw1(:)                   ! weight pointer
    character(len=*),parameter :: subname=' (med_merge_mod: med_merge)'
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
    ! Assume that input and output ranks are the same - this is checked in error check
    if (rank_out == 1) then
       call ESMF_FieldGet(field_in, farrayPtr=dpf1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_out, farrayPtr=dp1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (rank_out == 2) then
       call ESMF_FieldGet(field_in, farrayPtr=dpf2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_out, farrayPtr=dp2, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Get pointer to weights that weights are only rank 1
    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
       call ESMF_FieldBundleGet(FBw, fieldName=trim(fldw), field=field_wgt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_wgt, farrayPtr=dpw1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    ! Do supported merges
    if (trim(merge_type)  == 'copy') then
       if (rank_out == 1) then
          dp1(:) = dpf1(:)
       else
          dp2(:,:) = dpf2(:,:)
       endif
    else if (trim(merge_type)  == 'copy_with_weights') then
       if (rank_out == 1) then
          dp1(:) = dpf1(:)*dpw1(:)
       else
          do n = 1,ungriddedUBound_out(1)
             if (gridToFieldMap_out(1) == 1) then
                dp2(:,n) = dpf2(:,n)*dpw1(:)
             else if (gridToFieldMap_out(1) == 2) then
                dp2(n,:) = dpf2(n,:)*dpw1(:)
             end if
          end do
       endif
    else if (trim(merge_type)  == 'merge' .or. trim(merge_type) == 'sum_with_weights') then
       if (rank_out == 1) then
          dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
       else
          do n = 1,ungriddedUBound_out(1)
             if (gridToFieldMap_out(1) == 1) then
                dp2(:,n) = dp2(:,n) + dpf2(:,n)*dpw1(:)
             else if (gridToFieldMap_out(1) == 2) then
                dp2(n,:) = dp2(n,:) + dpf2(n,:)*dpw1(:)
             end if
          end do
       endif
    else if (trim(merge_type) == 'sum') then
       if (rank_out == 1) then
          dp1(:) = dp1(:) + dpf1(:)
       else
          dp2(:,:) = dp2(:,:) + dpf2(:,:)
       endif
    else
       call ESMF_LogWrite(trim(subname)//": merge type "//trim(merge_type)//" not supported", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine med_merge_auto_field

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
    character(len=*),parameter :: subname='(med_merge_fieldo_1d)'
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
          if (.not.FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
             call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif

          if (wgtfound) then
             if (.not.FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
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

  subroutine med_merge_field_2D(FBout, fnameout,     &
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

    ! input/output arguments
    type(ESMF_FieldBundle) , intent(inout)                 :: FBout
    character(len=*)       , intent(in)                    :: fnameout
    type(ESMF_FieldBundle) , intent(in)                    :: FBinA
    character(len=*)       , intent(in)                    :: fnameA
    real(R8)               , intent(in), pointer           :: wgtA(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinB
    character(len=*)       , intent(in), optional          :: fnameB
    real(R8)               , intent(in), optional, pointer :: wgtB(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinC
    character(len=*)       , intent(in), optional          :: fnameC
    real(R8)               , intent(in), optional, pointer :: wgtC(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinD
    character(len=*)       , intent(in), optional          :: fnameD
    real(R8)               , intent(in), optional, pointer :: wgtD(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinE
    character(len=*)       , intent(in), optional          :: fnameE
    real(R8)               , intent(in), optional, pointer :: wgtE(:,:)
    integer                , intent(out)                   :: rc

    ! local variables
    real(R8), pointer          :: dataOut(:,:)
    real(R8), pointer          :: dataPtr(:,:)
    real(R8), pointer          :: wgt(:,:)
    integer                    :: lb1,ub1,lb2,ub2,i,j,n
    logical                    :: wgtfound, FBinfound
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_merge_field_2d)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc=ESMF_SUCCESS

    if (.not. FB_FldChk(FBout, trim(fnameout), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//&
            trim(fnameout), ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    call FB_GetFldPtr(FBout, trim(fnameout), fldptr2=dataOut, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)
    lb2 = lbound(dataOut,2)
    ub2 = ubound(dataOut,2)

    dataOut = czero

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
          call FB_GetFldPtr(FBinA, trim(fnameA), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          wgtfound = .true.
          wgt => wgtA

       elseif (n == 2 .and. present(FBinB)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinB, trim(fnameB), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtB)) then
             wgtfound = .true.
             wgt => wgtB
          endif

       elseif (n == 3 .and. present(FBinC)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinC, trim(fnameC), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtC)) then
             wgtfound = .true.
             wgt => wgtC
          endif

       elseif (n == 4 .and. present(FBinD)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinD, trim(fnameD), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtD)) then
             wgtfound = .true.
             wgt => wgtD
          endif

       elseif (n == 5 .and. present(FBinE)) then
          FBinfound = .true.
          call FB_GetFldPtr(FBinE, trim(fnameE), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtE)) then
             wgtfound = .true.
             wgt => wgtE
          endif

       endif

       if (FBinfound) then
          if (.not.FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
             call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif

          if (wgtfound) then
             if (.not. FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
                call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
             do j = lb2,ub2
                do i = lb1,ub1
                   dataOut(i,j) = dataOut(i,j) + dataPtr(i,j) * wgt(i,j)
                enddo
             enddo
          else
             do j = lb2,ub2
                do i = lb1,ub1
                   dataOut(i,j) = dataOut(i,j) + dataPtr(i,j)
                enddo
             enddo
          endif  ! wgtfound

       endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_merge_field_2D

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
