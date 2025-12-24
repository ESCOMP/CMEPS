module med_field_info_mod

  !-----------------------------------------------------------------------------
  ! Defines a type and related operations for storing metadata about fields that can be
  ! used to create an ESMF FieldBundle.
  !-----------------------------------------------------------------------------

  use ESMF, only : ESMF_MAXSTR, ESMF_SUCCESS
  use ESMF, only : ESMF_Field, ESMF_State, ESMF_AttributeGet, ESMF_StateGet
  use med_utils_mod, only : ChkErr => med_utils_ChkErr
  use shr_log_mod, only : shr_log_error
  use shr_string_mod, only : shr_string_withoutSuffix
  use shr_wtracers_mod, only : WTRACERS_SUFFIX

  implicit none
  private

  !-----------------------------------------------
  ! Public methods
  !-----------------------------------------------

  public :: med_field_info_create  ! Create a single field
  public :: med_field_info_array_from_names_wtracers_ungridded  ! Create an array of field_info objects based on an array of names, where water tracers are given an ungridded dimension
  public :: med_field_info_array_from_state  ! Create an array of field_info objects based on the fields in an ESMF State

  !-----------------------------------------------
  ! Types
  !-----------------------------------------------

  type, public :: med_field_info_type
     character(ESMF_MAXSTR) :: name
     integer :: n_ungridded  ! number of ungridded dimensions

     ! These arrays will be allocated to be of size ungridded_count
     integer, allocatable :: ungridded_lbound(:)
     integer, allocatable :: ungridded_ubound(:)
  end type med_field_info_type

  character(len=*),parameter :: u_FILE_u = &
     __FILE__

!================================================================================
contains
!================================================================================

  function med_field_info_create(name, ungridded_lbound, ungridded_ubound, rc) result(field_info)
    ! Create a single field

    ! input/output variables
    character(len=*), intent(in) :: name

    ! ungridded_lbound and ungridded_ubound must either both be present or both be absent;
    ! if present, they must be the same size
    integer, intent(in), optional :: ungridded_lbound(:)
    integer, intent(in), optional :: ungridded_ubound(:)

    integer, intent(out) :: rc
    type(med_field_info_type) :: field_info  ! function result

    ! local variables
    integer :: n_ungridded
    character(len=*), parameter :: subname = '(med_field_info_create)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(ungridded_lbound) .neqv. present(ungridded_ubound)) then
       call shr_log_error( &
            subname//": ERROR: ungridded_lbound and ungridded_ubound must both be present or both absent.", &
            line=__LINE__, file=u_FILE_u, rc=rc)
       return
    end if

    field_info%name = name

    if (present(ungridded_lbound)) then
       n_ungridded = size(ungridded_lbound)
       if (size(ungridded_ubound) /= n_ungridded) then
          call shr_log_error( &
               subname//": ERROR: ungridded_lbound and ungridded_ubound must have the same size.", &
               line=__LINE__, file=u_FILE_u, rc=rc)
          return
       end if
       field_info%n_ungridded = n_ungridded
       allocate(field_info%ungridded_lbound(n_ungridded))
       allocate(field_info%ungridded_ubound(n_ungridded))
       field_info%ungridded_lbound = ungridded_lbound
       field_info%ungridded_ubound = ungridded_ubound
    else
       field_info%n_ungridded = 0
    end if

  end function med_field_info_create

  !-----------------------------------------------------------------------------

  subroutine med_field_info_array_from_names_wtracers_ungridded(field_names, field_info_array, rc)
    ! Create an array of field_info objects based on an array of names, where water
    ! tracers are given an ungridded dimension.
    !
    ! It is assumed that fields generally have no ungridded dimensions. However, for
    ! fields ending with the water tracer suffix, it is instead assumed that they have a
    ! single ungridded dimension of size given by shr_wtracers_get_num_tracers.

    ! input/output variables
    character(len=*), intent(in) :: field_names(:)
    type(med_field_info_type), allocatable, intent(out) :: field_info_array(:)
    integer, intent(out) :: rc

    ! local variables
    integer :: i, n_fields
    logical :: is_tracer
    integer :: n_tracers
    integer :: localrc
    character(len=*), parameter :: subname = '(med_field_info_array_from_names_wtracers_ungridded)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    n_fields = size(field_names)
    allocate(field_info_array(n_fields))
    ! For now, hard-code n_tracers, since we haven't set up the tracer information; we'll
    ! fix this in an upcoming set of changes
    n_tracers = 0

    do i = 1, n_fields
       call shr_string_withoutSuffix( &
            in_str = field_names(i), &
            suffix = WTRACERS_SUFFIX, &
            has_suffix = is_tracer, &
            rc = localrc)
       if (localrc /= 0) then
          call shr_log_error(subname//": ERROR in shr_string_withoutSuffix", rc=rc)
          return
       end if

       if (is_tracer) then
          ! Field is a water tracer; assume a single ungridded dimension
          field_info_array(i) = med_field_info_create( &
               name=field_names(i), &
               ungridded_lbound=[1], &
               ungridded_ubound=[n_tracers], &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          ! Not a water tracer; assume no ungridded dimensions
          field_info_array(i) = med_field_info_create( &
               name=field_names(i), &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

  end subroutine med_field_info_array_from_names_wtracers_ungridded

  subroutine med_field_info_array_from_state(state, field_info_array, rc)
    ! Create an array of field_info objects based on the Fields in an ESMF State

    ! input/output variables
    type(ESMF_State), intent(in) :: state
    type(med_field_info_type), allocatable, intent(out) :: field_info_array(:)
    integer, intent(out) :: rc

    ! local variables
    integer :: i, n_fields
    character(ESMF_MAXSTR), allocatable :: field_names(:)
    type(ESMF_Field) :: field
    logical :: is_present
    integer :: n_ungridded
    integer, allocatable :: ungridded_lbound(:)
    integer, allocatable :: ungridded_ubound(:)
    character(len=*), parameter :: subname = '(med_field_info_array_from_state)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount=n_fields, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(field_names(n_fields))
    allocate(field_info_array(n_fields))
    call ESMF_StateGet(state, itemNameList=field_names, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do i = 1, n_fields
       call ESMF_StateGet(state, itemName=trim(field_names(i)), field=field, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AttributeGet(field, name="UngriddedUBound", convention="NUOPC", &
            purpose="Instance", itemCount=n_ungridded,  isPresent=is_present, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (.not. is_present) then
          n_ungridded = 0
       end if

       if (n_ungridded == 0) then
          field_info_array(i) = med_field_info_create( &
               name=field_names(i), &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          allocate(ungridded_lbound(n_ungridded))
          allocate(ungridded_ubound(n_ungridded))
          call ESMF_AttributeGet(field, name="UngriddedLBound", convention="NUOPC", &
               purpose="Instance", valueList=ungridded_lbound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeGet(field, name="UngriddedUBound", convention="NUOPC", &
               purpose="Instance", valueList=ungridded_ubound, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_info_array(i) = med_field_info_create( &
               name=field_names(i), &
               ungridded_lbound=ungridded_lbound, &
               ungridded_ubound=ungridded_ubound, &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          deallocate(ungridded_lbound)
          deallocate(ungridded_ubound)
       end if
    end do

  end subroutine med_field_info_array_from_state

end module med_field_info_mod
