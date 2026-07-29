module med_ufs_trace_wrapper_mod

#ifdef UFS_TRACING
  use ufs_trace_mod, only: ufs_trace_init, ufs_trace, ufs_trace_finalize
#endif

  implicit none

  private

  public ufs_trace_init_wrapper
  public ufs_trace_wrapper
  public ufs_trace_finalize_wrapper

contains

  subroutine ufs_trace_init_wrapper()
#ifdef UFS_TRACING
    call ufs_trace_init
#endif
    return
  end subroutine ufs_trace_init_wrapper

  subroutine ufs_trace_wrapper(component, routine, ph)
    character(len=*), intent(in) :: component, routine, ph
#ifdef UFS_TRACING
    call ufs_trace(component, routine, ph)
#endif
    return
  end subroutine ufs_trace_wrapper

  subroutine ufs_trace_finalize_wrapper()
#ifdef UFS_TRACING
    call ufs_trace_finalize
#endif
  end subroutine ufs_trace_finalize_wrapper

end module med_ufs_trace_wrapper_mod
