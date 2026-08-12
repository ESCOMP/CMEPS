module med_utils_mod

  use med_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  public :: med_memcheck
  public :: med_utils_ChkErr
  public :: med_log_clock_advance
  public :: med_global_sums

  integer     , parameter :: memdebug_level=1
  character(*), parameter :: u_FILE_u = &
       __FILE__

  ! med_global_sums is generic over the number of fields being summed: the 1fld version takes
  ! summands for a single field and returns a scalar sum; the nflds version takes summands for
  ! an arbitrary number of fields and returns an array of sums.
  interface med_global_sums
     module procedure med_global_sums_1fld
     module procedure med_global_sums_nflds
  end interface med_global_sums

!===============================================================================
contains
!===============================================================================

  subroutine med_memcheck(string, level, maintask)
    character(len=*), intent(in) :: string
    integer, intent(in) :: level
    logical, intent(in) :: maintask
#ifdef CESMCOUPLED
    integer :: ierr
    integer, external :: GPTLprint_memusage
    if((maintask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif
#endif
  end subroutine med_memcheck

!===============================================================================

  logical function med_utils_ChkErr(rc, line, file, mpierr)
#ifndef NO_MPI2
    use mpi , only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
#else
    use mpi, only : MPI_SUCCESS
#endif
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_FAILURE, ESMF_LogWrite

    integer, intent(in) :: rc
    integer, intent(in) :: line

    character(len=*), intent(in) :: file
    logical, optional, intent(in) :: mpierr
#ifdef NO_MPI2
    integer, parameter :: MPI_MAX_ERROR_STRING=80
#else
    integer :: ierr, len
#endif
    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: lrc

    med_utils_ChkErr = .false.
    lrc = rc
    if (present(mpierr)) then
       if(mpierr) then
          if (rc == MPI_SUCCESS) return
#ifdef NO_MPI2
          write(lstring,*) "ERROR in mct mpi-serial library rc=",rc
#else
          call MPI_ERROR_STRING(rc, lstring, len, ierr)
#endif
          call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file)
          lrc = ESMF_FAILURE
       endif
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
      med_utils_ChkErr = .true.
    endif

  end function med_utils_ChkErr

!===============================================================================

  subroutine med_log_clock_advance(clock, component, logunit)
    use ESMF, only : ESMF_Clock, ESMF_ClockPrint

    type(ESMF_Clock) :: clock
    character(len=*), intent(in) :: component
    integer, intent(in) :: logunit

    character(len=CL) :: cvalue, prestring
    integer :: rc

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, &
         preString=trim(prestring), rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine med_log_clock_advance

!===============================================================================

  subroutine med_global_sums_nflds(gcomp, local_summands, global_sums, rc)

    ! Compute global sums of the given local summands.
    !
    ! local_summands has dimensions [number of local grid cells, number of fields];
    ! global_sums has dimensions [number of fields].
    !
    ! If the bfbflag attribute is true, the sums are computed in a manner that is
    ! independent of processor count; otherwise a cheaper, processor-count-dependent
    ! algorithm is used. Note that, for the sums to truly be independent of processor
    ! count, the caller needs to provide the individual summands here rather than
    ! pre-summing them locally.

    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF  , only : ESMF_VMAllreduce, ESMF_REDUCE_SUM, ESMF_SUCCESS
    use NUOPC , only : NUOPC_CompAttributeGet
    use shr_log_mod, only : shr_log_error
#ifdef CESMCOUPLED
    use shr_reprosum_mod, only : shr_reprosum_calc
#endif

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    real(r8), intent(in)  :: local_summands(:,:)
    real(r8), intent(out) :: global_sums(:)
    integer , intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer           :: nsummands  ! number of local grid cells
    integer           :: nflds      ! number of fields
    integer           :: n, nf
    integer           :: mpicom     ! MPI communicator of this component
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    logical           :: bfbflag    ! value of the bfbflag attribute
    real(r8) :: local_sums(size(local_summands, 2))  ! local sums over grid cells, for each field

    character(len=*), parameter :: subname = '(med_global_sums_nflds)'
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='bfbflag', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    bfbflag = .false.
    if (isPresent .and. isSet) then
       read(cvalue,*) bfbflag
    end if

    nsummands = size(local_summands, 1)
    nflds = size(local_summands, 2)

    if (size(global_sums) /= nflds) then
       call shr_log_error(subname//' ERROR: size of global_sums must agree with the '// &
            'second dimension of local_summands', line=__LINE__, file=u_FILE_u, rc=rc)
       return
    end if

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (bfbflag) then
#ifdef CESMCOUPLED
       call ESMF_VMGet(vm, mpiCommunicator=mpicom, rc=rc)
       if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_reprosum_calc(arr=local_summands, arr_gsum=global_sums, &
            nsummands=nsummands, dsummands=nsummands, nflds=nflds, commid=mpicom)
#else
       call shr_log_error(subname//' ERROR: bfbflag is true, but this build does not have '// &
            'access to the shr_reprosum module needed for bfbflag to be true', &
            line=__LINE__, file=u_FILE_u, rc=rc)
       return
#endif
    else
       do nf = 1, nflds
          local_sums(nf) = 0.0_r8
          do n = 1, nsummands
             local_sums(nf) = local_sums(nf) + local_summands(n,nf)
          end do
       end do
       call ESMF_VMAllreduce(vm, senddata=local_sums, recvdata=global_sums, count=nflds, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_global_sums_nflds

!===============================================================================

  subroutine med_global_sums_1fld(gcomp, local_summands, global_sum, rc)

    ! Compute the global sum of the given local summands, for a single field.
    !
    ! This is a convenience wrapper around med_global_sums_nflds for the common case of a
    ! single field; see that routine for details.

    use ESMF , only : ESMF_GridComp, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    real(r8)           , intent(in)  :: local_summands(:)
    real(r8)           , intent(out) :: global_sum
    integer            , intent(out) :: rc

    ! local variables
    ! Note that summands_2d is deliberately allocatable, rather than being an automatic array or
    ! an inline reshape of local_summands: those would give an array temporary, which some
    ! compilers place on the stack, and local_summands can be large (e.g., of order the number of
    ! local grid cells).
    real(r8), allocatable :: summands_2d(:,:)
    real(r8)              :: global_sums(1)
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    allocate(summands_2d(size(local_summands), 1))
    summands_2d(:,1) = local_summands(:)
    call med_global_sums_nflds(gcomp, summands_2d, global_sums, rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    global_sum = global_sums(1)

  end subroutine med_global_sums_1fld

end module med_utils_mod
