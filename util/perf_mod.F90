module perf_mod

  use ESMF

  implicit none
  public

  public t_startf
  public t_stopf

contains

   subroutine t_startf(event, handle)
     character(len=*), intent(in) :: event
     integer,  optional :: handle
     call ESMF_TraceRegionEnter(event)
   end subroutine t_startf

   subroutine t_stopf(event, handle)
     character(len=*), intent(in) :: event
     integer,  optional :: handle
     call ESMF_TraceRegionExit(event)
   end subroutine t_stopf

end module perf_mod
