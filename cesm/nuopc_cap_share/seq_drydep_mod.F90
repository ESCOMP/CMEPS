module seq_drydep_mod

  use shr_drydep_mod, only: seq_drydep_setHCoeff=>shr_drydep_setHCoeff
  implicit none

  ! method specification
  character(len=*), parameter :: DD_XLND = 'xactive_lnd' ! dry-dep land
  character(len=*), parameter :: drydep_method = DD_XLND ! XLND is the only option now
  logical, protected :: lnd_drydep

contains

  subroutine seq_drydep_readnl(NLFilename, drydep_nflds)

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer, intent(out)          :: drydep_nflds

    call shr_drydep_readnl(NLFilename, drydep_nflds)

    lnd_drydep = drydep_nflds>0

  end subroutine seq_drydep_readnl

end module seq_drydep_mod
