module perf_mod

  implicit none
  public

  public t_startf
  public t_stopf

contains

   subroutine t_startf(event, handle)
     character(len=*), intent(in) :: event
     integer,  optional :: handle
     ! Do nothing for nems right now
   end subroutine t_startf

   subroutine t_stopf(event, handle)
     character(len=*), intent(in) :: event
     integer,  optional :: handle
     ! Do nothing for nems right now
   end subroutine t_stopf

end module perf_mod
