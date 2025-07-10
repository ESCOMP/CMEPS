module med_phases_ocnalb_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : InternalState, logunit
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr          => med_utils_chkerr
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_internalstate_mod , only : mapconsf, mapnames, compatm, compocn, maintask
  use perf_mod              , only : t_startf, t_stopf
  use shr_orb_mod           , only : shr_orb_cosz, shr_orb_decl
  use shr_orb_mod           , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use shr_log_mod           , only : shr_log_unit, shr_log_error
  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_phases_ocnalb_run

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private med_phases_ocnalb_init
  private med_phases_ocnalb_orbital_update
  private med_phases_ocnalb_orbital_init

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type ocnalb_type
     real(r8) , pointer :: lats  (:) => null() ! latitudes  (degrees)
     real(r8) , pointer :: lons  (:) => null() ! longitudes (degrees)
     integer  , pointer :: mask  (:) => null() ! ocn domain mask: 0 <=> inactive cell
     real(r8) , pointer :: anidr (:) => null() ! albedo: near infrared, direct
     real(r8) , pointer :: avsdr (:) => null() ! albedo: visible      , direct
     real(r8) , pointer :: anidf (:) => null() ! albedo: near infrared, diffuse
     real(r8) , pointer :: avsdf (:) => null() ! albedo: visible      , diffuse
     logical            :: created   ! has memory been allocated here
  end type ocnalb_type

  character(*),parameter :: u_FILE_u = &
       __FILE__
  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  ! used, reused in module
  logical  :: flux_albav      ! use average dif and dir albedos
  logical  :: use_nextswcday  ! use the scalar field for next time (otherwise, will be set using clock)
  logical  :: use_min_albedo  ! apply minimum value of albedo for direct vis, nir
  real(R8) :: min_albedo      ! minimum value of albedo for direct vis, nir
  real(R8) :: albdif          ! 60 deg reference albedo, diffuse
  real(R8) :: albdir          ! 60 deg reference albedo, direct
!===============================================================================
contains
!===============================================================================

  subroutine med_phases_ocnalb_init(gcomp, ocnalb, rc)

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables and then use the module
    ! variables in the med_ocnalb phase
    ! All input field bundles are ASSUMED to be on the ocean grid
    !-----------------------------------------------------------------------

    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF  , only : ESMF_VM, ESMF_VMGet, ESMF_Mesh, ESMF_MeshGet
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet
    use NUOPC , only : NUOPC_CompAttributeGet
    use ESMF  , only : operator(==)

    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(ocnalb_type) , intent(inout) :: ocnalb
    integer           , intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_VM)            :: vm
    integer                  :: iam
    type(ESMF_Field)         :: lfield
    type(ESMF_Mesh)          :: lmesh
    integer                  :: n
    integer                  :: lsize
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    type(InternalState)      :: is_local
    real(R8), pointer        :: ownedElemCoords(:)
    character(len=CL)        :: tempc1,tempc2
    character(len=CS)        :: cvalue
    logical                  :: use_min_ocnalb
    logical                  :: isPresent, isSet
    integer                  :: fieldCount
    character(CL)            :: msg
    type(ESMF_Field), pointer :: fieldlist(:)
    character(*), parameter  :: subname = '(med_phases_ocnalb_init) '
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Set pointers to fields needed for albedo calculations
    !----------------------------------

    ! These must must be on the ocean grid since the ocean albedo computation is on the ocean grid
    ! The following sets pointers to the module arrays

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldname='So_avsdr', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=ocnalb%avsdr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldname='So_avsdf', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=ocnalb%avsdf, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldname='So_anidr', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=ocnalb%anidr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldname='So_anidf', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=ocnalb%anidf, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Get lat, lon, which are time-invariant
    !----------------------------------

    ! The following assumes that all fields in FBMed_ocnalb_o have the same grid - so
    ! only need to query field 1
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnalb_o, fieldlist=fieldlist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(fieldlist(1), mesh=lmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(fieldlist)
    call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(ocnalb%anidr)
    if (numOwnedElements /= lsize) then
       write(tempc1,'(i10)') numOwnedElements
       write(tempc2,'(i10)') lsize
       call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements "// trim(tempc1) // &
            " not equal to local size "// trim(tempc2), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(ocnalb%lons(numOwnedElements))
    allocate(ocnalb%lats(numOwnedElements))
    call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,lsize
       ocnalb%lons(n) = ownedElemCoords(2*n-1)
       ocnalb%lats(n) = ownedElemCoords(2*n)
    end do

    ! Initialize orbital values
    call  med_phases_ocnalb_orbital_init(gcomp, logunit, iam==0, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine if reference albedos are used
    flux_albav = .false.
    call NUOPC_CompAttributeGet(gcomp, name='flux_albav', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flux_albav
    end if
    ! Set reference albedo values
    call NUOPC_CompAttributeGet(gcomp, name="albdif", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) albdif
    else
       albdif = 0.06_r8
    end if
    call NUOPC_CompAttributeGet(gcomp, name="albdir", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) albdir
    else
       albdir = 0.07_r8
    end if
    ! Determine if direct albedo should have a minimum value
    call NUOPC_CompAttributeGet(gcomp, name="ocean_albedo_limit", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) min_albedo
       use_min_albedo = .true.
    else
       min_albedo = 0.0_R8
       use_min_ocnalb = .false.
    endif
    ! Allow setting of albedo timestep using the clock instead of the atm's next timestep
    use_nextswcday = .true.
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent ) then
       use_nextswcday = .false.
    endif

    if (flux_albav) then
       write(msg,'(2(A,f8.2))') trim(subname)//': mean albedos set: albdif = ',albdif,', albdir = ',albdir
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    else
       if (use_min_albedo) then
          write(msg,'(A,f8.2)') trim(subname)//': min_albedo setting = ',min_albedo
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end if
    end if
    write(msg,'(A,l1)') trim(subname)//': use_nextswcday setting is ',use_nextswcday
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnalb_init

  !===============================================================================

  subroutine med_phases_ocnalb_run(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Compute ocean albedos (on the ocean grid)
    !-----------------------------------------------------------------------

    use NUOPC_Mediator, only : NUOPC_MediatorGet
    use ESMF          , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_TimeInterval
    use ESMF          , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF          , only : ESMF_ClockIsCreated, ESMF_ClockGetNextTime
    use ESMF          , only : ESMF_VM, ESMF_VMGet
    use ESMF          , only : ESMF_LogWrite, ESMF_LogFoundError
    use ESMF          , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO
    use ESMF          , only : ESMF_Field, ESMF_FieldGet
    use ESMF          , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use ESMF          , only : operator(+)
    use NUOPC         , only : NUOPC_CompAttributeGet
    use med_constants_mod , only : shr_const_pi
    use med_phases_history_mod, only : med_phases_history_write_med

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    ! local variables
    type(ocnalb_type), save :: ocnalb
    type(ESMF_VM)           :: vm
    type(ESMF_Field)        :: lfield
    integer                 :: iam
    logical                 :: update_alb
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: clock
    type(ESMF_Clock)        :: dclock
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    character(CL)           :: cvalue
    character(CS)           :: starttype        ! config start type
    character(CL)           :: runtype          ! initial, continue, hybrid, branch
    real(R8)                :: nextsw_cday      ! calendar day of next atm shortwave
    real(R8), pointer       :: ofrac(:)
    real(R8), pointer       :: ofrad(:)
    real(R8), pointer       :: ifrac(:)
    real(R8), pointer       :: ifrad(:)
    integer                 :: lsize            ! local size
    integer                 :: n                ! indices
    real(R8)                :: rlat             ! gridcell latitude in radians
    real(R8)                :: rlon             ! gridcell longitude in radians
    real(R8)                :: cosz             ! Cosine of solar zenith angle
    real(R8)                :: eccen            ! Earth orbit eccentricity
    real(R8)                :: mvelpp           ! Earth orbit
    real(R8)                :: lambm0           ! Earth orbit
    real(R8)                :: obliqr           ! Earth orbit
    real(R8)                :: delta            ! Solar declination angle  in radians
    real(R8)                :: eccf             ! Earth orbit eccentricity factor
    real(R8), parameter     :: const_deg2rad = shr_const_pi/180.0_R8  ! deg to rads
    character(CL)           :: msg
    logical                 :: first_call = .true.
    character(len=*)  , parameter :: subname='(med_phases_ocnalb_run)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Determine main task
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine if ocnalb data type will be initialized - and if not return
    if (first_call) then
       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o, rc=rc)) then
          ocnalb%created = .true.
       else
          ocnalb%created = .false.
       end if
    end if
    if (.not. ocnalb%created) then
       return
    end if

    ! Note that in the mct version the atm was initialized first so
    ! that nextsw_cday could be passed to the other components - this
    ! assumed that atmosphere component was ALWAYS initialized first.
    ! In the nuopc version it will be easier to assume that on startup
    ! - nextsw_cday is just what cam was setting it as the current calendar day

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    call t_startf('MED:'//subname)

    ! get clock
    call ESMF_GridCompGet(gcomp, clock=clock)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then

       ! Initialize ocean albedo calculation
       call med_phases_ocnalb_init(gcomp, ocnalb, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) starttype

       if (trim(starttype) == trim('startup')) then
          runtype = "initial"
       else if (trim(starttype) == trim('continue') ) then
          runtype = "continue"
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
       else
          call ESMF_LogWrite( subname//' ERROR: unknown starttype', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (trim(runtype) == 'initial') then
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          ! obtain nextsw_cday from atm if it is in the import state
          if (use_nextswcday) then
             call State_GetScalar(&
                  state=is_local%wrap%NstateImp(compatm), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, &
                  scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
                  scalar_value=nextsw_cday, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       first_call = .false.

    else
       ! Note that med_methods_State_GetScalar includes a broadcast to all other pets
       if (use_nextswcday) then
          call State_GetScalar(&
               state=is_local%wrap%NstateImp(compatm), &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               flds_scalar_num=is_local%wrap%flds_scalar_num, &
               scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
               scalar_value=nextsw_cday, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Get orbital values
    call med_phases_ocnalb_orbital_update(clock, logunit, iam==0, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Calculate ocean albedos on the ocean grid
    update_alb = .false.
    lsize = size(ocnalb%anidr)

    if (flux_albav) then
       do n = 1,lsize
          ocnalb%anidr(n) = albdir
          ocnalb%avsdr(n) = albdir
          ocnalb%anidf(n) = albdif
          ocnalb%avsdf(n) = albdif
       end do
       update_alb = .true.
    else
       ! Solar declination
       ! Will only do albedo calculation if nextsw_cday is not -1.
       write(msg,*)trim(subname)//' nextsw_cday = ',nextsw_cday
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

       if (nextsw_cday >= -0.5_r8) then
          call shr_orb_decl(nextsw_cday, eccen, mvelpp,lambm0, obliqr, delta, eccf)

          ! Compute albedos
          do n = 1,lsize
             rlat = const_deg2rad * ocnalb%lats(n)
             rlon = const_deg2rad * ocnalb%lons(n)
             cosz = shr_orb_cosz( nextsw_cday, rlat, rlon, delta )
             if (cosz  >  0.0_r8) then !--- sun hit --
                ocnalb%anidr(n) = (.026_r8/(cosz**1.7_r8 + 0.065_r8)) +   &
                                  (.150_r8*(cosz         - 0.100_r8 ) *   &
                                  (cosz - 0.500_r8 ) * (cosz - 1.000_r8 )  )
                if (use_min_albedo) then
                   ocnalb%anidr(n) = max (ocnalb%anidr(n), min_albedo)
                end if
                ocnalb%avsdr(n) = ocnalb%anidr(n)
                ocnalb%anidf(n) = albdif
                ocnalb%avsdf(n) = albdif
             else !--- dark side of earth ---
                ocnalb%anidr(n) = 1.0_r8
                ocnalb%avsdr(n) = 1.0_r8
                ocnalb%anidf(n) = 1.0_r8
                ocnalb%avsdf(n) = 1.0_r8
             end if
          end do
          update_alb = .true.

       endif    ! nextsw_cday
    end if   ! flux_albav

    ! Update current ifrad/ofrad values if albedo was updated in field bundle
    if (update_alb) then
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ifrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=ifrac, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ifrad', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=ifrad, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ofrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=ofrac, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ofrad', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=ofrad, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ifrad(:) = ifrac(:)
       ofrad(:) = ofrac(:)
    endif

    if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o, rc=rc)) then
       call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_ClockIsCreated(dclock)) then
          call med_phases_history_write_med(gcomp, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBMed_ocnalb_o, string=trim(subname)//' FBMed_ocnalb_o', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnalb_run

!===============================================================================

  subroutine med_phases_ocnalb_orbital_init(gcomp, logunit, maintask, rc)

    !----------------------------------------------------------
    ! Obtain orbital related values
    !----------------------------------------------------------

    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_LogWrite, ESMF_LogFoundError, ESMF_LogSetError
    use ESMF  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    integer             , intent(in)    :: logunit         ! output logunit
    logical             , intent(in)    :: maintask
    integer             , intent(out)   :: rc              ! output error

    ! local variables

    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(med_phases_ocnalb_orbital_init)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
          write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call shr_log_error(msgstr, line=__LINE__, file=__FILE__, rc=rc)
          return
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call shr_log_error(msgstr, line=__LINE__, file=__FILE__, rc=rc)
          return
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
          write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
          write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call shr_log_error(msgstr, line=__LINE__, file=__FILE__, rc=rc)
          return
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call shr_log_error(msgstr, line=__LINE__, file=__FILE__, rc=rc)
       return
    endif
  end subroutine med_phases_ocnalb_orbital_init

  !===============================================================================

  subroutine med_phases_ocnalb_orbital_update(clock, logunit,  maintask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    use ESMF, only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF, only : ESMF_LogSetError, ESMF_RC_NOT_VALID, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: maintask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox long of perihelion plus pi (rad)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    integer, save     :: prev_orb_year=0
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    character(len=*) , parameter :: subname = "(med_phases_ocnalb_orbital_update)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    lprint = .false.
    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
    else
       orb_year = orb_iyear
    end if
    first_time = .false.

    eccen = orb_eccen
    shr_log_unit = logunit

    if(orb_year .ne. prev_orb_year) then
       prev_orb_year = orb_year
       lprint = maintask
    end if

    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call shr_log_error(msgstr, line=__LINE__, file=__FILE__, rc=rc)
       return
    endif

  end subroutine med_phases_ocnalb_orbital_update

!===============================================================================

end module med_phases_ocnalb_mod
