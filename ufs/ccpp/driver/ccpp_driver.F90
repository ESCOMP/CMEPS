module ccpp_driver

  use ccpp_api,           only: ccpp_t

  implicit none
  private

  public ccpp_step

  type(ccpp_t), pointer :: cdata => null()
  integer :: nthrds

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine ccpp_step(step, nblks, ierr)

    ! input/output variables
    character(len=*), intent(in)  :: step
    integer,          intent(in)  :: nblks
    integer,          intent(out) :: ierr

    ! local variables
    integer :: nb, nt
    character(len=*), parameter :: subname='(ccpp_step)'
    !-----------------------------------------------------------

    ierr = 0

    if (trim(step)=="init") then
       ! set number of threads
       ! TODO: also support OpenMP threading
       nthrds = 1

       ! allocate cdata structures for blocks and threads
       if (.not. allocated(cdata_block)) allocate(cdata_block(1:nblks,1:nthrds))

       ! loop over all blocks and threads
       do nt=1, nthrds
          do nb=1, nblks
             ! assign the correct block and thread numbers
             cdata_block(nb,nt)%blk_no = nb
             cdata_block(nb,nt)%thrd_no = nt
          end do
       end do
    end if

  end subroutine ccpp_step

end module ccpp_driver
