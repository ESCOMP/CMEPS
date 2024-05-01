module shr_dust_emis_mod

  !========================================================================
  ! Module for handling dust emissions.
  ! This module is shared by land and atmosphere models for the computation of
  ! dust emissions.
  !========================================================================

  use shr_sys_mod    , only : shr_sys_abort
  use shr_kind_mod   , only : CS => SHR_KIND_CS
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_log_mod    , only : shr_log_getLogUnit, errMsg => shr_log_errMsg

  implicit none
  private

  ! public member functions
  public :: shr_dust_emis_readnl           ! Read namelist
  public :: is_dust_emis_zender            ! If Zender_2003 dust emission method is being used
  public :: is_dust_emis_leung             ! If Leungr_2023 dust emission method is being used
  public :: is_zender_soil_erod_from_land  ! If Zender_2003 is being used and soil eroditability is in land
  public :: is_zender_soil_erod_from_atm   ! If Zender/_2003 is being used and soil eroditability is in atmosphere

  ! The following is only public for the sake of unit testing; it should not be called
  ! directly outside this module
  public :: dust_emis_set_options       ! Set the namelist options directory not through the namelist
  public :: is_NOT_initialized          ! Check if dust emission has NOT been initialized

  ! private data members:
  private :: check_options_finish_init  ! Check that the options are correct and finish initialization

  ! PRIVATE DATA:
  character(len=CS) :: dust_emis_method = 'Zender_2003'  ! Dust emisison method to use: Zender_2003 or Leung_2023
  character(len=CS) :: zender_soil_erod_source = 'none'  ! if calculated in lnd or atm (only when Zender_2003 is used)
  logical           :: dust_emis_initialized=.false.     ! If dust emissions have been initiatlized yet or not

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
CONTAINS
!===============================================================================

  subroutine shr_dust_emis_readnl(mpicom, NLFilename)

    !========================================================================
    ! reads dust_emis_inparm namelist to determine how dust emissions will
    ! be handled between the land and atmosphere models
    !========================================================================
    use shr_mpi_mod, only : shr_mpi_bcast, shr_mpi_commrank

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer         , intent(in)  :: mpicom     ! MPI communicator for broadcasting all all tasks

    !----- local -----
    integer       :: unitn            ! namelist unit number
    integer       :: ierr             ! error code
    logical       :: exists           ! if file exists or not
    integer       :: localPet         ! Local processor rank
    integer       :: s_logunit        ! Output log unit
    character(*),parameter :: F00   = "('(shr_dust_emis_read) ',8a)"
    character(*),parameter :: subName = '(shr_dust_emis_read) '
    !-----------------------------------------------------------------------------

    namelist /dust_emis_inparm/ dust_emis_method, zender_soil_erod_source

    !-----------------------------------------------------------------------------
    ! Read namelist, check if namelist file exists first
    !-----------------------------------------------------------------------------

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0  )then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if

    call shr_mpi_commrank( mpicom, localPet )

    call shr_log_getLogUnit(s_logunit)
    if (localPet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,F00) 'Read in dust_emis_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'dust_emis_inparm', ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0, no namelist is present.
             read(unitn, dust_emis_inparm, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort( subName//'ERROR:: problem on read of dust_emis_inparm ' &
                                    // 'namelist in shr_dust_emis_readnl')
             end if
          endif
          close( unitn )
       end if
    end if
    call shr_mpi_bcast(dust_emis_method, mpicom)
    call shr_mpi_bcast(zender_soil_erod_source, mpicom)

    call check_options_finish_init()

  end subroutine shr_dust_emis_readnl

!====================================================================================

  subroutine check_options_finish_init()
   ! Some error checking and mark initialization as finished
   integer :: s_logunit        ! Output log unit
   character(*),parameter :: subName = '(check_options_finish_init) '

   call shr_log_getLogUnit(s_logunit)
   if (trim(dust_emis_method) == 'Leung_2023') then
      if ( trim(zender_soil_erod_source) /= 'none' )then
         write(s_logunit,*) 'ERROR: '//errMsg(u_FILE_u, __LINE__)
         call shr_sys_abort(subName//"ERROR: zender_soil_erod_source should NOT be set, when dust_emis_method=Leung_2023" )
         return
      end if
   else if (trim(dust_emis_method) == 'Zender_2003') then
      if ( (trim(zender_soil_erod_source) /= 'lnd') .and. (trim(zender_soil_erod_source) /= 'atm') )then
         write(s_logunit,*) 'zender_soil_erod_source is NOT valid = ', trim(zender_soil_erod_source)
         write(s_logunit,*) 'ERROR: '//errMsg(u_FILE_u, __LINE__)
         call shr_sys_abort(subName//"ERROR: zender_soil_erod_source can only be lnd or atm" )
         return
      end if
   else
      write(s_logunit,*) 'dust_emis_method not recognized = ', trim(dust_emis_method)
      write(s_logunit,*) 'ERROR: '//errMsg(u_FILE_u, __LINE__)
      call shr_sys_abort(subName//"ERROR: dust_emis_method namelist item is not valid" )
      return
   end if

   dust_emis_initialized = .true.

  end subroutine check_options_finish_init

!====================================================================================

  logical function is_dust_emis_zender()
     ! is_dust_emis_zender – Logical function, true if the Zender 2003 scheme is being used
     if ( is_NOT_initialized() ) return
     if (trim(dust_emis_method) == 'Zender_2003') then
        is_dust_emis_zender = .true.
     else
        is_dust_emis_zender = .false.
     end if
  end function is_dust_emis_zender

!===============================================================================

  logical function is_dust_emis_leung()
     ! is_dust_emis_leung – Logical function, true if the Leung 2023 scheme is being used
     if ( is_NOT_initialized() ) return
     if (trim(dust_emis_method) == 'Leung_2023') then
        is_dust_emis_leung = .true.
     else
        is_dust_emis_leung = .false.
     end if
  end function is_dust_emis_leung

!===============================================================================

  logical function is_zender_soil_erod_from_land()
     ! is_zender_soil_erod_from_land – Logical function, true if the Zender method is being used and soil erodibility is in CTSM
     if ( is_NOT_initialized() ) return
     if ( is_dust_emis_zender() )then
        if (trim(zender_soil_erod_source) == 'lnd') then
           is_zender_soil_erod_from_land = .true.
        else
           is_zender_soil_erod_from_land = .false.
        end if
     else
        is_zender_soil_erod_from_land = .false.
     end if
  end function is_zender_soil_erod_from_land

!===============================================================================

  logical function is_zender_soil_erod_from_atm()
     !is_zender_soil_erod_from_atm – Logical function, true if the Zender method is being used and soil erodibility is in CAM
     if ( is_NOT_initialized() ) return
     if ( is_dust_emis_zender() )then
        if ( trim(zender_soil_erod_source) == 'atm') then
           is_zender_soil_erod_from_atm = .true.
        else
           is_zender_soil_erod_from_atm = .false.
        end if
     else
        is_zender_soil_erod_from_atm = .false.
     end if
  end function is_zender_soil_erod_from_atm

!===============================================================================

  logical function is_NOT_initialized()
     ! Check if this is NOT initialized and return true if so (false if initialized)
     ! Will abort with an error when using in the model
     ! For unit testing will return the logical state
     integer       :: s_logunit        ! Output log unit

     if ( dust_emis_initialized )then
        is_NOT_initialized = .false.
        return
     else
        is_NOT_initialized = .true.
        call shr_log_getLogUnit(s_logunit)
        write(s_logunit,*) 'ERROR: '//errMsg(u_FILE_u, __LINE__)
        call shr_sys_abort( 'ERROR: dust emission namelist has NOT been read in yet,' // &
                            ' shr_dust_emis_mod is NOT initialized ' )
     end if
   end function is_NOT_initialized

  subroutine dust_emis_set_options( dust_emis_method_in, zender_soil_erod_source_in)
    character(len=*), intent(IN) :: dust_emis_method_in         ! Dust emisison method to use: Zender_2003 or Leung_2023
    character(len=*), intent(IN) :: zender_soil_erod_source_in  ! if calculed in lnd or atm (only when Zender_2003 is used)

    dust_emis_method = dust_emis_method_in
    zender_soil_erod_source = zender_soil_erod_source_in
    call check_options_finish_init()
  end subroutine dust_emis_set_options

!===============================================================================

end module shr_dust_emis_mod
