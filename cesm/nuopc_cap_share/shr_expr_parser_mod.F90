!=============================================================================
! expression parser utility --
!   for parsing simple linear mathematical expressions of the form
!   X = a*R + b*S + c*(X + Y + Z) ...
!
!=============================================================================
module shr_expr_parser_mod
  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : CXX => shr_kind_cxx

  implicit none
  private

  public :: shr_exp_parse     ! parses simple strings which contain expressions
  public :: shr_exp_item_t    ! user defined type which contains an expression component
  public :: shr_exp_list_destroy ! destroy the linked list returned by shr_exp_parse

  ! contains componets of expression
  type shr_exp_item_t
     character(len=64) :: name
     character(len=64),pointer :: vars(:) => null()
     real(r8)         ,pointer :: coeffs(:) => null()
     integer           :: n_terms = 0
     type(shr_exp_item_t), pointer :: next_item => null()
  end type shr_exp_item_t

contains

  ! -----------------------------------------------------------------
  ! parses expressions provided in array of strings
  ! -----------------------------------------------------------------
  function shr_exp_parse( exp_array, nitems ) result(exp_items_list)

    character(len=*), intent(in)   :: exp_array(:) ! contains a expressions
    integer, optional, intent(out) :: nitems       ! number of expressions parsed
    type(shr_exp_item_t), pointer  :: exp_items_list ! linked list of items returned

    integer :: i,j, n_exp_items
    type(shr_exp_item_t), pointer :: exp_item, list_item
    integer :: ndxs(512)
    integer :: nelem, j1,j2,k
    character(len=CXX) :: tmp_str, tmp_name
    character(len=8) :: xchr ! multipler
    real(r8) :: xdbl
    real(r8) :: coeff0
    logical :: more_to_come
    character(len=CXX), allocatable :: sums_grps(:)
    character(len=CXX) :: sum_string

    allocate(sums_grps(size(exp_array)))

    nullify( exp_items_list )
    nullify( exp_item )
    nullify( list_item )

    sums_grps(:) = ' '

    ! combine lines that have a trailing "+" with the next line
    i=1
    j=1
    loop1: do while( len_trim(exp_array(i)) > 0 )

       k = scan(exp_array(i), '+', back=.true. )
       more_to_come = k == len_trim(exp_array(i)) ! line ends with "+"

       if ( more_to_come ) then
          sums_grps(j) = trim(sums_grps(j)) // trim(adjustl(exp_array(i)))
       else
          sums_grps(j) = trim(sums_grps(j)) // trim(adjustl(exp_array(i)))
          j = j+1
       endif

       i = i+1
       if ( i > size(exp_array) ) exit loop1

    end do loop1

    n_exp_items = j-1

    ! a group is  a summation of terms

    ! parse the individual sum strings...  and form the groupings
    has_grps: if (n_exp_items>0) then

       ! from shr_megan_mod ... should be generalized and shared...
       grploop: do i = 1,n_exp_items

          ! parse out the term names
          ! from first parsing out the terms in the summation equation ("+" separates the terms)

          sum_string = sums_grps(i)
          j = scan( sum_string, '=' )
          nelem = 1
          ndxs(nelem) = j ! ndxs stores the index of each term of the equation

          ! find indices of all the terms in the equation
          tmp_str = trim( sum_string(j+1:) )
          j = scan( tmp_str, '+' )
          do while(j>0)
             nelem = nelem+1
             ndxs(nelem) = ndxs(nelem-1) + j
             tmp_str = tmp_str(j+1:)
             j = scan( tmp_str, '+' )
          enddo
          ndxs(nelem+1) = len(sum_string)+1

          allocate( exp_item )

          exp_item%n_terms = nelem ! number of terms

          exp_item%name = trim(adjustl( sum_string(:ndxs(1)-1))) ! thing to the left of the "=" is used as the name of the group

          ! now that we have the number of terms in the summation allocate memory for the terms
          allocate( exp_item%vars(nelem) )
          allocate( exp_item%coeffs(nelem) )

          coeff0 = 1._r8 ! default multiplier

          ! now parse out the multiplier from the terms
          elmloop: do k = 1,nelem

             exp_item%coeffs(k) = coeff0

             ! get the term name which follows the '*' operator if the is one
             tmp_name = adjustl(sum_string(ndxs(k)+1:ndxs(k+1)-1))

             j = scan( tmp_name, '*' )
             if (j>0) then

                xchr = tmp_name(1:j-1) ! get the multipler (left of the '*')
                read( xchr, * ) xdbl   ! convert the string to a real
                exp_item%coeffs(k) = xdbl ! store the multiplier

                j1 = scan( tmp_name, '(' )
                if (j1>0) then
                   coeff0 = xdbl
                   tmp_name = trim(adjustl(tmp_name(j1+1:))) ! get the term name (right of the '*')
                else
                   coeff0 = 1._r8
                   tmp_name = trim(adjustl(tmp_name(j+1:))) ! get the term name (right of the '*')
                endif

             endif

             j2 = scan( tmp_name, ')' )
             if (j2>0) then
                coeff0 = 1._r8
                tmp_name = tmp_name(1:j2-1)
             endif

             exp_item%vars(k) = trim(tmp_name)

          enddo elmloop

          if (associated(exp_item)) then
             if (associated(exp_items_list)) then
                list_item => exp_items_list
                do while(associated(list_item%next_item))
                   list_item => list_item%next_item
                enddo
                list_item%next_item => exp_item
             else
                exp_items_list => exp_item
             endif
          endif


       enddo grploop
    endif has_grps

    if ( present(nitems) ) then
       nitems = n_exp_items
    endif

    deallocate(sums_grps)

  end function shr_exp_parse

  ! -----------------------------------------------------------------
  ! deallocates memory occupied by linked list
  ! -----------------------------------------------------------------
  subroutine  shr_exp_list_destroy( list )
    type(shr_exp_item_t), pointer, intent(inout) :: list

    type(shr_exp_item_t), pointer :: item, next

    item => list
    do while(associated(item))
       next => item%next_item
       if (associated(item%vars)) then
          deallocate(item%vars)
          nullify(item%vars)
          deallocate(item%coeffs)
          nullify(item%coeffs)
       endif
       deallocate(item)
       nullify(item)
       item => next
    enddo

  end subroutine  shr_exp_list_destroy

end module shr_expr_parser_mod
