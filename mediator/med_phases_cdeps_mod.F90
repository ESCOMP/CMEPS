module med_phases_cdeps_mod

  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_LogWrite
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGMSG_INFO

  use dshr_strdata_mod, only: shr_strdata_type
  use dshr_strdata_mod, only: shr_strdata_init_from_inline
  use perf_mod        , only: t_startf, t_stopf

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interf aces
  !--------------------------------------------------------------------------

  public med_phases_cdeps_run

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------



  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------



  character(*),parameter :: u_FILE_u = __FILE__

!============================================================================
contains
!============================================================================

  subroutine med_phases_cdeps_run(gcomp, rc)

    !------------------------------------------------------------------------
    ! Use CDEPS inline capability to read in data
    !------------------------------------------------------------------------

    use ESMF, only : ESMF_GridComp

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*)  , parameter :: subname='(med_phases_cdeps_run)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    !if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    !endif


    !if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    !endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_cdeps_run

end module med_phases_cdeps_mod
