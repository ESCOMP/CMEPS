module med_field_info_mod

  !-----------------------------------------------------------------------------
  ! Defines a type and related operations for storing metadata about fields that can be
  ! used to create an ESMF FieldBundle.
  !-----------------------------------------------------------------------------

  use ESMF            , only : ESMF_MAXSTR, ESMF_SUCCESS, ESMF_TYPEKIND_R8
  use ESMF            , only : ESMF_Field, ESMF_State, ESMF_StateGet
  use ESMF            , only : ESMF_Mesh, ESMF_MeshLoc
  use ESMF            , only : ESMF_FieldCreate, ESMF_FieldGet
  use med_utils_mod   , only : ChkErr => med_utils_ChkErr
  use shr_log_mod     , only : shr_log_error
  use wtracers_mod    , only : wtracers_is_wtracer_field

  implicit none
  private

  !-----------------------------------------------
  ! Public methods
  !-----------------------------------------------

  ! Create a single field_info object from direct specification of values
  public :: med_field_info_create_directly

  ! Create a single field_info object from information in an ESMF_Field
  public :: med_field_info_create_from_field

  ! Create an array of field_info objects based on an array of names, where water tracers
  ! are treated specially (being given an ungridded dimension)
  public :: med_field_info_array_from_names_wtracers

  ! Create an array of field_info objects based on the fields in an ESMF State
  public :: med_field_info_array_from_state

  ! Create an ESMF Field (using ESMF_FieldCreate) based on a field_info object
  public :: med_field_info_esmf_fieldcreate

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

  function med_field_info_create_directly(name, ungridded_lbound, ungridded_ubound, rc) result(field_info)
    ! Create a single field_info object from direct specification of values

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
    character(len=*), parameter :: subname = '(med_field_info_create_directly)'
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

  end function med_field_info_create_directly

  !-----------------------------------------------------------------------------

  function med_field_info_create_from_field(field, name, rc) result(field_info)
    ! Create a single field_info object from information in an ESMF_Field

    ! input/output variables
    ! We get information other than the name from this ESMF_Field object
    type(ESMF_Field), intent(in) :: field

    ! We should be able to get the name from the field, but in all current uses of this
    ! function, we already have the name available, so it's easy enough to just pass it in
    ! rather than making this function query it again. If future users did not already
    ! have the name readily available, we could either change this to optional or remove
    ! it entirely and just always get the name from querying the field.
    character(len=*), intent(in) :: name

    integer, intent(out) :: rc
    type(med_field_info_type) :: field_info  ! function result

    ! local variables
    integer :: n_ungridded
    integer, allocatable :: ungridded_lbound(:)
    integer, allocatable :: ungridded_ubound(:)

    character(len=*), parameter :: subname = '(med_field_info_create_from_field)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, ungriddedDimCount=n_ungridded, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (n_ungridded == 0) then
       field_info = med_field_info_create_directly( &
            name=name, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(ungridded_lbound(n_ungridded))
       allocate(ungridded_ubound(n_ungridded))
       call ESMF_FieldGet(field, &
            ungriddedLBound=ungridded_lbound, ungriddedUBound=ungridded_ubound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       field_info = med_field_info_create_directly( &
            name=name, &
            ungridded_lbound=ungridded_lbound, &
            ungridded_ubound=ungridded_ubound, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       deallocate(ungridded_lbound)
       deallocate(ungridded_ubound)
    end if
  end function med_field_info_create_from_field

  !-----------------------------------------------------------------------------

  subroutine med_field_info_array_from_names_wtracers(field_names, field_info_array, rc)
    ! Create an array of field_info objects based on an array of names, where water
    ! tracers are treated specially (being given an ungridded dimension).
    !
    ! It is assumed that fields generally have no ungridded dimensions. However, for
    ! fields ending with the water tracer suffix, it is instead assumed that they have a
    ! single ungridded dimension of size given by shr_wtracers_get_num_tracers.
    !
    ! field_info_array is allocated here (and, since it has intent(out), it is
    ! automatically deallocated if it is already allocated on entry to this subroutine)

    ! input/output variables
    character(len=*), intent(in) :: field_names(:)
    type(med_field_info_type), allocatable, intent(out) :: field_info_array(:)
    integer, intent(out) :: rc

    ! local variables
    integer :: i, n_fields
    logical :: is_tracer
    integer :: n_tracers
    character(len=*), parameter :: subname = '(med_field_info_array_from_names_wtracers)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    n_fields = size(field_names)
    allocate(field_info_array(n_fields))
    ! For now, hard-code n_tracers, since we haven't set up the tracer information; we'll
    ! fix this in an upcoming set of changes
    n_tracers = 0

    do i = 1, n_fields
       is_tracer = wtracers_is_wtracer_field(field_names(i))
       if (is_tracer) then
          ! Field is a water tracer; assume a single ungridded dimension
          field_info_array(i) = med_field_info_create_directly( &
               name=field_names(i), &
               ungridded_lbound=[1], &
               ungridded_ubound=[n_tracers], &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          ! Not a water tracer; assume no ungridded dimensions
          field_info_array(i) = med_field_info_create_directly( &
               name=field_names(i), &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

  end subroutine med_field_info_array_from_names_wtracers

  !-----------------------------------------------------------------------------

  subroutine med_field_info_array_from_state(state, field_info_array, rc)
    ! Create an array of field_info objects based on the Fields in an ESMF State
    !
    ! field_info_array is allocated here (and, since it has intent(out), it is
    ! automatically deallocated if it is already allocated on entry to this subroutine)

    ! input/output variables
    type(ESMF_State), intent(in) :: state
    type(med_field_info_type), allocatable, intent(out) :: field_info_array(:)
    integer, intent(out) :: rc

    ! local variables
    integer :: i, n_fields
    character(ESMF_MAXSTR), allocatable :: field_names(:)
    type(ESMF_Field) :: field
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

       field_info_array(i) = med_field_info_create_from_field( &
            field=field, &
            name=field_names(i), &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine med_field_info_array_from_state

  !-----------------------------------------------------------------------------

  subroutine med_field_info_esmf_fieldcreate(field_info, mesh, meshloc, field, rc)
    ! Create an ESMF Field (using ESMF_FieldCreate) based on a field_info object

    ! input/output variables
    type(med_field_info_type), intent(in) :: field_info
    type(ESMF_Mesh), intent(in) :: mesh
    type(ESMF_MeshLoc), intent(in) :: meshloc
    type(ESMF_Field), intent(out) :: field
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = '(med_field_info_esmf_fieldcreate)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (field_info%n_ungridded > 0) then
       field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=meshloc, &
            name=field_info%name, &
            ungriddedLbound=field_info%ungridded_lbound, &
            ungriddedUbound=field_info%ungridded_ubound, &
            gridToFieldMap=[field_info%n_ungridded+1], &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=meshloc, &
            name=field_info%name, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_field_info_esmf_fieldcreate

end module med_field_info_mod
