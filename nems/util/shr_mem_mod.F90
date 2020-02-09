module shr_mem_mod

  use shr_kind_mod, only : shr_kind_r8

  implicit none
  public

contains

  subroutine shr_mem_getusage(r_msize, r_mrss, prt)
    real(shr_kind_r8) :: r_msize,r_mrss
    logical, optional :: prt
    ! For now does nothing

  end subroutine shr_mem_getusage

end module shr_mem_mod
