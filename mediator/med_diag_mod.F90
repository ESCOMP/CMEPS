module med_diag_mod

  !----------------------------------------------------------------------------
  ! Compute spatial and time averages of fluxed quatities for water and
  ! energy balance
  !
  ! Sign convention for fluxes is positive downward with hierarchy being
  !    atm/glc/lnd/rof/ice/ocn
  ! Sign convention:
  !    positive value <=> the model is gaining water, heat, momentum, etc.
  ! Unit convention:
  !    heat flux     ~ W/m^2
  !    momentum flux ~ N/m^2
  !    water flux    ~ (kg/s)/m^2
  !    salt  flux    ~ (kg/s)/m^2
  !----------------------------------------------------------------------------

  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
  use NUOPC_Mediator        , only : NUOPC_MediatorGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR
  use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
  use ESMF                  , only : ESMF_VM, ESMF_VMReduce, ESMF_REDUCE_SUM
  use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockGetNextTime
  use ESMF                  , only : ESMF_Alarm, ESMF_ClockGetAlarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
  use ESMF                  , only : ESMF_FieldBundle, ESMF_Field, ESMF_FieldGet
  use med_constants_mod     , only : shr_const_rearth, shr_const_pi, shr_const_latice, shr_const_latvap
  use med_constants_mod     , only : shr_const_ice_ref_sal, shr_const_ocn_ref_sal, shr_const_isspval
  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : InternalState, logunit, maintask, diagunit, samegrid_atmlnd
  use med_methods_mod       , only : fldbun_getdata2d => med_methods_FB_getdata2d
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_fldChk    => med_methods_FB_FldChk
  use med_time_mod          , only : alarmInit        => med_time_alarmInit
  use med_utils_mod         , only : chkerr           => med_utils_ChkErr
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_diag_init
  public  :: med_diag_zero
  public  :: med_phases_diag_accum
  public  :: med_phases_diag_print
  public  :: med_phases_diag_atm
  public  :: med_phases_diag_lnd
  public  :: med_phases_diag_rof
  public  :: med_phases_diag_glc
  public  :: med_phases_diag_ocn
  public  :: med_phases_diag_ice_ice2med
  public  :: med_phases_diag_ice_med2ice

  private :: med_diag_sum_main
  private :: med_diag_print_atm
  private :: med_diag_print_lnd_ice_ocn
  private :: med_diag_print_summary

  type, public :: budget_diag_type
     character(CS) :: name
  end type budget_diag_type

  type, public :: budget_diag_indices
     type(budget_diag_type), pointer :: comps(:) => null()
     type(budget_diag_type), pointer :: fields(:) => null()
     type(budget_diag_type), pointer :: periods(:) => null()
  end type budget_diag_indices
  type(budget_diag_indices) :: budget_diags

  interface med_diag_zero
     module procedure med_diag_zero_mode
     module procedure med_diag_zero_select
  end interface

  ! ---------------------------------
  ! print options (obtained from mediator config input)
  ! ---------------------------------

  ! sets the diagnotics level of the annual budgets. [0,1,2,3],
  !  0 = none,
  !  1 = net summary budgets
  !  2 = 1 + detailed lnd/ocn/ice component budgets
  !  3 = 2 + detailed atm budgets

  integer :: budget_print_inst  ! default is 0
  integer :: budget_print_daily ! default is 0
  integer :: budget_print_month ! default is 1
  integer :: budget_print_ann   ! default is 1
  integer :: budget_print_ltann ! default is 1
  integer :: budget_print_ltend ! default is 0

  ! formats for output tables
  character(*), parameter :: F00  = "('(med_phases_diag_print) ',4a)"
  character(*), parameter :: FAH  = "(4a,i9,i6)"
  character(*), parameter :: FA0  = "('    ',12x,6(6x,a8,1x))"
  character(*), parameter :: FA1  = "('    ',a12,6f15.8)"
  character(*), parameter :: FA0r = "('    ',12x,8(6x,a8,1x))"
  character(*), parameter :: FA1r = "('    ',a12,8f15.8)"
  character(*), parameter :: FA0s = "('    ',12x,8(7x,a8,2x))"
  character(*), parameter :: FA1s = "('    ',a12,8g18.8)"

  ! ---------------------------------
  ! C for component
  ! ---------------------------------

  ! "r" is receive by the mediator from the component
  ! "s" is send from the mediator  to the component

  integer  :: c_atm_send  ! model index: atm
  integer  :: c_atm_recv  ! model index: atm
  integer  :: c_inh_send  ! model index: ice, northern
  integer  :: c_inh_recv  ! model index: ice, northern
  integer  :: c_ish_send  ! model index: ice, southern
  integer  :: c_ish_recv  ! model index: ice, southern
  integer  :: c_lnd_send  ! model index: lnd
  integer  :: c_lnd_recv  ! model index: lnd
  integer  :: c_ocn_send  ! model index: ocn
  integer  :: c_ocn_recv  ! model index: ocn
  integer  :: c_rof_send  ! model index: rof
  integer  :: c_rof_recv  ! model index: rof
  integer  :: c_glc_send  ! model index: glc
  integer  :: c_glc_recv  ! model index: glc

  ! The folowing is needed for detailing the atm budgets and breakdown into components
  integer  :: c_inh_asend  ! model index: ice, northern, on atm grid
  integer  :: c_inh_arecv  ! model index: ice, northern, on atm grid
  integer  :: c_ish_asend  ! model index: ice, southern, on atm grid
  integer  :: c_ish_arecv  ! model index: ice, southern, on atm grid
  integer  :: c_lnd_asend  ! model index: lnd, on atm grid
  integer  :: c_lnd_arecv  ! model index: lnd, on atm grid
  integer  :: c_ocn_asend  ! model index: ocn, on atm grid
  integer  :: c_ocn_arecv  ! model index: ocn, on atm grid

  ! ---------------------------------
  ! F for field
  ! ---------------------------------
  integer, parameter :: unset_index = -999
  integer :: f_area          = unset_index ! area (wrt to unit sphere)
  integer :: f_heat_frz      = unset_index ! heat : latent, freezing
  integer :: f_heat_melt     = unset_index ! heat : latent, melting
  integer :: f_heat_swnet    = unset_index ! heat : short wave, net
  integer :: f_heat_lwdn     = unset_index ! heat : longwave down
  integer :: f_heat_lwup     = unset_index ! heat : longwave up
  integer :: f_heat_latvap   = unset_index ! heat : latent, vaporization
  integer :: f_heat_latf     = unset_index ! heat : latent, fusion, snow
  integer :: f_heat_ioff     = unset_index ! heat : latent, fusion, frozen runoff
  integer :: f_heat_sen      = unset_index ! heat : sensible
  integer :: f_heat_rain     = unset_index ! heat : heat content of rain
  integer :: f_heat_snow     = unset_index ! heat : heat content of snow
  integer :: f_heat_evap     = unset_index ! heat : heat content of evaporation
  integer :: f_heat_cond     = unset_index ! heat : heat content of evaporation
  integer :: f_heat_rofl     = unset_index ! heat : heat content of liquid runoff
  integer :: f_heat_rofi     = unset_index ! heat : heat content of ice runoff

  integer :: f_watr_frz      = unset_index ! water: freezing
  integer :: f_watr_melt     = unset_index ! water: melting
  integer :: f_watr_rain     = unset_index ! water: precip, liquid
  integer :: f_watr_snow     = unset_index ! water: precip, frozen
  integer :: f_watr_evap     = unset_index ! water: evaporation
  integer :: f_watr_salt     = unset_index ! water: water equivalent of salt flux
  integer :: f_watr_roff     = unset_index ! water: runoff/flood
  integer :: f_watr_ioff     = unset_index ! water: frozen runoff
  integer :: f_watr_frz_16O  = unset_index ! water isotope: freezing
  integer :: f_watr_melt_16O = unset_index ! water isotope: melting
  integer :: f_watr_rain_16O = unset_index ! water isotope: precip, liquid
  integer :: f_watr_snow_16O = unset_index ! water isotope: prcip, frozen
  integer :: f_watr_evap_16O = unset_index ! water isotope: evaporation
  integer :: f_watr_roff_16O = unset_index ! water isotope: runoff/flood
  integer :: f_watr_ioff_16O = unset_index ! water isotope: frozen runoff
  integer :: f_watr_frz_18O  = unset_index ! water isotope: freezing
  integer :: f_watr_melt_18O = unset_index ! water isotope: melting
  integer :: f_watr_rain_18O = unset_index ! water isotope: precip, liquid
  integer :: f_watr_snow_18O = unset_index ! water isotope: precip, frozen
  integer :: f_watr_evap_18O = unset_index ! water isotope: evaporation
  integer :: f_watr_roff_18O = unset_index ! water isotope: runoff/flood
  integer :: f_watr_ioff_18O = unset_index ! water isotope: frozen runoff
  integer :: f_watr_frz_HDO  = unset_index ! water isotope: freezing
  integer :: f_watr_melt_HDO = unset_index ! water isotope: melting
  integer :: f_watr_rain_HDO = unset_index ! water isotope: precip, liquid
  integer :: f_watr_snow_HDO = unset_index ! water isotope: precip, frozen
  integer :: f_watr_evap_HDO = unset_index ! water isotope: evaporation
  integer :: f_watr_roff_HDO = unset_index ! water isotope: runoff/flood
  integer :: f_watr_ioff_HDO = unset_index ! water isotope: frozen runoff

  integer :: f_heat_beg      = unset_index ! 1st index  for heat
  integer :: f_heat_end      = unset_index ! Last index for heat
  integer :: f_watr_beg      = unset_index ! 1st index  for water
  integer :: f_watr_end      = unset_index ! Last index for water
  integer :: f_salt_beg      = unset_index ! 1st index  for salt
  integer :: f_salt_end      = unset_index ! Last index for salt

  integer :: f_16O_beg       = unset_index ! 1st index  for 16O water isotope
  integer :: f_16O_end       = unset_index ! Last index for 16O water isotope
  integer :: f_18O_beg       = unset_index ! 1st index  for 18O water isotope
  integer :: f_18O_end       = unset_index ! Last index for 18O water isotope
  integer :: f_HDO_beg       = unset_index ! 1st index  for HDO water isotope
  integer :: f_HDO_end       = unset_index ! Last index for HDO water isotope

  ! ---------------------------------
  ! water isotopes names and indices
  ! ---------------------------------

  logical :: flds_wiso  = .false.! If water isotope fields are active -
  ! TODO: for now set to .false. - but this needs to be set in an initialization phase

  integer, parameter :: nisotopes = 3
  integer :: iso0(nisotopes)
  integer :: isof(nisotopes)
  character(len=5) :: isoname(nisotopes)

  ! ---------------------------------
  ! P for period
  ! ---------------------------------

  integer :: period_inst=0
  integer :: period_day=0
  integer :: period_mon=0
  integer :: period_ann=0
  integer :: period_inf=0

  ! ---------------------------------
  ! local constants
  ! ---------------------------------

  real(r8), parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
       &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
       &    (shr_const_ocn_ref_sal*shr_const_latice)

  ! WFLX (kg/m^2s) = -SFLX (kg/m^2s) / ocn_ref_sal (psu) (34.7g/kg) / 1.e-3 kg/g
  real(r8), parameter :: SFLXtoWFLX = & ! water flux implied by salt flux (kg/m^2s)
       -1._r8/(shr_const_ocn_ref_sal*1.e-3_r8)

  ! ---------------------------------
  ! public data members
  ! ---------------------------------

  ! note: call med_diag_sum_main then save budget_global and budget_counter on restart from/to root pe ---

  real(r8), allocatable :: budget_local  (:,:,:) ! local sum, valid on all pes
  real(r8), allocatable :: budget_global (:,:,:) ! global sum, valid only on root pe
  real(r8), allocatable :: budget_counter(:,:,:) ! counter, valid only on root pe
  real(r8), allocatable :: budget_global_1d(:)   ! needed for ESMF_VMReduce call

  character(len=*), parameter :: modName   = "(med_diag) "
  character(len=*), parameter :: u_FILE_u  = &
      __FILE__

  character(len=CS) :: budget_table_version

!===============================================================================
contains
!===============================================================================

  subroutine med_diag_init(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Initialize module variables and allocate dynamic memory
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(out)   :: rc

    ! local variables

    integer           :: c_size   ! number of component send/recvs
    integer           :: f_size   ! number of fields
    integer           :: p_size   ! number of period types
    character(CS)     :: cvalue
    logical           :: isPresent, isSet
    character(*), parameter :: subName = '(med_phases_diag_init) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if(maintask) then
       write(logunit,'(a)') ' Creating budget_diags%comps '
    end if

    call NUOPC_CompAttributeGet(gcomp, name="budget_table_version", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (isPresent .and. isSet) then
       read(cvalue,*) budget_table_version
    else
       budget_table_version = 'v1'
    end if
    if (maintask) then
       write(logunit,'(a)') trim(subname) //' budget table version is '//trim(budget_table_version)
    end if

    call add_to_budget_diag(budget_diags%comps, c_atm_send , 'c2a_atm' ) ! comp index: atm
    call add_to_budget_diag(budget_diags%comps, c_atm_recv , 'a2c_atm' ) ! comp index: atm
    call add_to_budget_diag(budget_diags%comps, c_inh_send , 'c2i_inh' ) ! comp index: ice, northern
    call add_to_budget_diag(budget_diags%comps, c_inh_recv , 'i2c_inh' ) ! comp index: ice, northern
    call add_to_budget_diag(budget_diags%comps, c_ish_send , 'c2i_ish' ) ! comp index: ice, southern
    call add_to_budget_diag(budget_diags%comps, c_ish_recv , 'i2c_ish' ) ! comp index: ice, southern
    call add_to_budget_diag(budget_diags%comps, c_lnd_send , 'c2l_lnd' ) ! comp index: lnd
    call add_to_budget_diag(budget_diags%comps, c_lnd_recv , 'l2c_lnd' ) ! comp index: lnd
    call add_to_budget_diag(budget_diags%comps, c_ocn_send , 'c2o_ocn' ) ! comp index: ocn
    call add_to_budget_diag(budget_diags%comps, c_ocn_recv , 'o2c_ocn' ) ! comp index: ocn
    call add_to_budget_diag(budget_diags%comps, c_rof_send , 'c2r_rof' ) ! comp index: rof
    call add_to_budget_diag(budget_diags%comps, c_rof_recv , 'r2c_rof' ) ! comp index: rof
    call add_to_budget_diag(budget_diags%comps, c_glc_send , 'c2g_glc' ) ! comp index: glc
    call add_to_budget_diag(budget_diags%comps, c_glc_recv , 'g2c_glc' ) ! comp index: glc
    call add_to_budget_diag(budget_diags%comps, c_inh_asend, 'c2a_inh' ) ! comp index: ice, northern, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_inh_arecv, 'a2c_inh' ) ! comp index: ice, northern, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_ish_asend, 'c2a_ish' ) ! comp index: ice, southern, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_ish_arecv, 'a2c_ish' ) ! comp index: ice, southern, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_lnd_asend, 'c2a_lnd' ) ! comp index: lnd, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_lnd_arecv, 'a2c_lnd' ) ! comp index: lnd, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_ocn_asend, 'c2a_ocn' ) ! comp index: ocn, on atm grid
    call add_to_budget_diag(budget_diags%comps, c_ocn_arecv, 'a2c_ocn' ) ! comp index: ocn, on atm grid

    call add_to_budget_diag(budget_diags%fields, f_area          ,'area'        ) ! field  area (wrt to unit sphere)

    ! -----------------------------------------
    ! Heat fluxes budget terms
    ! -----------------------------------------

    ! Note that this order is important here to determine f_heat_beg and f_heat_end
    call add_to_budget_diag(budget_diags%fields, f_heat_frz      ,'hfreeze'     ) ! field  heat : latent, freezing
    call add_to_budget_diag(budget_diags%fields, f_heat_melt     ,'hmelt'       ) ! field  heat : latent, melting
    call add_to_budget_diag(budget_diags%fields, f_heat_swnet    ,'hnetsw'      ) ! field  heat : short wave, net
    call add_to_budget_diag(budget_diags%fields, f_heat_lwdn     ,'hlwdn'       ) ! field  heat : longwave down
    call add_to_budget_diag(budget_diags%fields, f_heat_lwup     ,'hlwup'       ) ! field  heat : longwave up
    call add_to_budget_diag(budget_diags%fields, f_heat_latvap   ,'hlatvap'     ) ! field  heat : latent, vaporization
    call add_to_budget_diag(budget_diags%fields, f_heat_latf     ,'hlatfus'     ) ! field  heat : latent, fusion, snow
    call add_to_budget_diag(budget_diags%fields, f_heat_ioff     ,'hiroff'      ) ! field  heat : latent, fusion, frozen runoff
    call add_to_budget_diag(budget_diags%fields, f_heat_sen      ,'hsen'        ) ! field  heat : sensible
    if (trim(budget_table_version) == 'v0') then
       f_heat_beg = f_heat_frz      ! field  first index for heat
       f_heat_end = f_heat_sen      ! field  last  index for heat
    else if (trim(budget_table_version) == 'v1') then
       call add_to_budget_diag(budget_diags%fields, f_heat_rain  ,'hrain'       ) ! field  heat : enthalpy of rain
       call add_to_budget_diag(budget_diags%fields, f_heat_snow  ,'hsnow'       ) ! field  heat : enthalpy of snow
       call add_to_budget_diag(budget_diags%fields, f_heat_evap  ,'hevap'       ) ! field  heat : enthalpy of evaporation
       call add_to_budget_diag(budget_diags%fields, f_heat_cond  ,'hcond'       ) ! field  heat : enthalpy of evaporation
       call add_to_budget_diag(budget_diags%fields, f_heat_rofl  ,'hrofl'       ) ! field  heat : enthalpy of liquid runoff
       call add_to_budget_diag(budget_diags%fields, f_heat_rofi  ,'hrofi'       ) ! field  heat : enthalpy of ice runoff
       f_heat_beg = f_heat_frz      ! field  first index for heat
       f_heat_end = f_heat_rofi     ! field  last  index for heat
    end if

    ! -----------------------------------------
    ! Water fluxes budget terms
    ! -----------------------------------------

    ! Note that this order is important here to determine f_watr_beg and f_watr_end
    if (trim(budget_table_version) == 'v0') then
       call add_to_budget_diag(budget_diags%fields, f_watr_frz   ,'wfreeze'     ) ! field  water: freezing
    end if
    call add_to_budget_diag(budget_diags%fields, f_watr_melt     ,'wmelt'       ) ! field  water: melting
    call add_to_budget_diag(budget_diags%fields, f_watr_rain     ,'wrain'       ) ! field  water: precip, liquid
    call add_to_budget_diag(budget_diags%fields, f_watr_snow     ,'wsnow'       ) ! field  water: precip, frozen
    call add_to_budget_diag(budget_diags%fields, f_watr_evap     ,'wevap'       ) ! field  water: evaporation
    if (trim(budget_table_version) == 'v0') then
       call add_to_budget_diag(budget_diags%fields, f_watr_salt  ,'weqsaltf'    ) ! field  water: water equivalent of salt flux
    endif
    call add_to_budget_diag(budget_diags%fields, f_watr_roff     ,'wrunoff'     ) ! field  water: runoff/flood
    call add_to_budget_diag(budget_diags%fields, f_watr_ioff     ,'wfrzrof'     ) ! field  water: frozen runoff
    if (trim(budget_table_version) == 'v0') then
       f_watr_beg = f_watr_frz  ! field  firs  index for water
    else
       f_watr_beg = f_watr_melt ! field  firs  index for water
    end if
    f_watr_end = f_watr_ioff    ! field  last  index for water

    if (flds_wiso) then
       call add_to_budget_diag(budget_diags%fields, f_watr_frz_16O  ,'wfreeze_16O' ) ! field  water isotope: freezing
       call add_to_budget_diag(budget_diags%fields, f_watr_melt_16O ,'wmelt_16O'   ) ! field  water isotope: melting
       call add_to_budget_diag(budget_diags%fields, f_watr_rain_16O ,'wrain_16O'   ) ! field  water isotope: precip, liquid
       call add_to_budget_diag(budget_diags%fields, f_watr_snow_16O ,'wsnow_16O'   ) ! field  water isotope: prcip, frozen
       call add_to_budget_diag(budget_diags%fields, f_watr_evap_16O ,'wevap_16O'   ) ! field  water isotope: evaporation
       call add_to_budget_diag(budget_diags%fields, f_watr_roff_16O ,'wrunoff_16O' ) ! field  water isotope: runoff/flood
       call add_to_budget_diag(budget_diags%fields, f_watr_ioff_16O ,'wfrzrof_16O' ) ! field  water isotope: frozen runoff
       f_16O_beg  = f_watr_frz_16O  ! field 1st  index for 16O water isotope
       f_16O_end  = f_watr_ioff_16O ! field Last index for 16O water isotope

       call add_to_budget_diag(budget_diags%fields, f_watr_frz_18O  ,'wfreeze_18O' ) ! field  water isotope: freezing
       call add_to_budget_diag(budget_diags%fields, f_watr_melt_18O ,'wmelt_18O'   ) ! field  water isotope: melting
       call add_to_budget_diag(budget_diags%fields, f_watr_rain_18O ,'wrain_18O'   ) ! field  water isotope: precip, liquid
       call add_to_budget_diag(budget_diags%fields, f_watr_snow_18O ,'wsnow_18O'   ) ! field  water isotope: precip, frozen
       call add_to_budget_diag(budget_diags%fields, f_watr_evap_18O ,'wevap_18O'   ) ! field  water isotope: evaporation
       call add_to_budget_diag(budget_diags%fields, f_watr_roff_18O ,'wrunoff_18O' ) ! field  water isotope: runoff/flood
       call add_to_budget_diag(budget_diags%fields, f_watr_ioff_18O ,'wfrzrof_18O' ) ! field  water isotope: frozen runoff
       f_18O_beg  = f_watr_frz_18O  ! field 1st  index for 18O water isotope
       f_18O_end  = f_watr_ioff_18O ! field Last index for 18O water isotope

       call add_to_budget_diag(budget_diags%fields, f_watr_frz_HDO  ,'wfreeze_HDO' ) ! field  water isotope: freezing
       call add_to_budget_diag(budget_diags%fields, f_watr_melt_HDO ,'wmelt_HDO'   ) ! field  water isotope: melting
       call add_to_budget_diag(budget_diags%fields, f_watr_rain_HDO ,'wrain_HDO'   ) ! field  water isotope: precip, liquid
       call add_to_budget_diag(budget_diags%fields, f_watr_snow_HDO ,'wsnow_HDO'   ) ! field  water isotope: precip, frozen
       call add_to_budget_diag(budget_diags%fields, f_watr_evap_HDO ,'wevap_HDO'   ) ! field  water isotope: evaporation
       call add_to_budget_diag(budget_diags%fields, f_watr_roff_HDO ,'wrunoff_HDO' ) ! field  water isotope: runoff/flood
       call add_to_budget_diag(budget_diags%fields, f_watr_ioff_HDO ,'wfrzrof_HDO' ) ! field  water isotope: frozen runoff
       f_HDO_beg  = f_watr_frz_HDO  ! field 1st  index for HDO water isotope
       f_HDO_end  = f_watr_ioff_HDO ! field Last index for HDO water isotope

       ! water isotopes
       iso0(:)    = (/ f_16O_beg, f_18O_beg, f_hdO_beg /)
       isof(:)    = (/ f_16O_end, f_18O_end, f_hdO_end /)
       isoname(:) = (/ 'H216O',   'H218O',   '  HDO'   /)
    end if

    ! -----------------------------------------
    ! Salt fluxes budget terms (for v1 only)
    ! -----------------------------------------

    if (trim(budget_table_version) == 'v1') then
       call add_to_budget_diag(budget_diags%fields, f_watr_salt  ,'saltf') ! field  water: salt flux
       f_salt_beg = f_watr_salt
       f_salt_end = f_watr_salt
    endif

    !-------------------------------------------------------------------------------
    ! Get config variables
    !-------------------------------------------------------------------------------

    budget_print_inst = get_diag_attribute(gcomp, 'budget_inst', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_print_daily = get_diag_attribute(gcomp, 'budget_daily', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_print_month = get_diag_attribute(gcomp, 'budget_month', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_print_ann = get_diag_attribute(gcomp, 'budget_ann', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_print_ltann = get_diag_attribute(gcomp, 'budget_ltann', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_print_ltend = get_diag_attribute(gcomp, 'budget_ltend', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! period types
    call add_to_budget_diag(budget_diags%periods, period_inst,'    inst')
    if(budget_print_daily > 0) call add_to_budget_diag(budget_diags%periods, period_day ,'   daily')
    if(budget_print_month > 0) call add_to_budget_diag(budget_diags%periods, period_mon ,' monthly')
    if(budget_print_ann   > 0) call add_to_budget_diag(budget_diags%periods, period_ann ,'  annual')
    call add_to_budget_diag(budget_diags%periods, period_inf ,'all_time')

    ! allocate module budget arrays
    c_size = size(budget_diags%comps)
    f_size = size(budget_diags%fields)
    p_size = size(budget_diags%periods)

    allocate(budget_local    (f_size , c_size , p_size)) ! local sum, valid on all pes
    allocate(budget_global   (f_size , c_size , p_size)) ! global sum, valid only on root pe
    allocate(budget_counter  (f_size , c_size , p_size)) ! counter, valid only on root pe
    allocate(budget_global_1d(f_size * c_size * p_size)) ! needed for ESMF_VMReduce call

  end subroutine med_diag_init

  integer function get_diag_attribute(gcomp, name, rc)
    type(ESMF_GridComp) , intent(inout) :: gcomp
    character(len=*), intent(in) :: name
    integer, intent(out) :: rc

    character(CS)     :: cvalue
    logical :: isPresent

    rc = ESMF_SUCCESS
    get_diag_attribute = 0
    call NUOPC_CompAttributeGet(gcomp, name=name, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name=name, value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) get_diag_attribute
    else
       call NUOPC_CompAttributeAdd(gcomp, (/name/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name=name, value='0', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
  end function get_diag_attribute

  !===============================================================================
  subroutine med_diag_zero_mode(mode, rc)

    ! ------------------------------------------------------------------
    ! Zero out global budget diagnostic data.
    ! ------------------------------------------------------------------

    ! input/output variables
    character(len=*) , intent(in)  :: mode
    integer          , intent(out) :: rc

    ! local variables
    character(*), parameter :: subName = '(med_diag_zero) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (trim(mode) == 'inst') then
       budget_local(:,:,period_inst) = 0.0_r8
       budget_global(:,:,period_inst) = 0.0_r8
       budget_counter(:,:,period_inst) = 0.0_r8
    elseif (trim(mode) == 'day') then
       budget_local(:,:,period_day) = 0.0_r8
       budget_global(:,:,period_day) = 0.0_r8
       budget_counter(:,:,period_day) = 0.0_r8
    elseif (trim(mode) == 'mon') then
       budget_local(:,:,period_mon) = 0.0_r8
       budget_global(:,:,period_mon) = 0.0_r8
       budget_counter(:,:,period_mon) = 0.0_r8
    elseif (trim(mode) == 'ann') then
       budget_local(:,:,period_ann) = 0.0_r8
       budget_global(:,:,period_ann) = 0.0_r8
       budget_counter(:,:,period_ann) = 0.0_r8
    elseif (trim(mode) == 'inf') then
       budget_local(:,:,period_inf) = 0.0_r8
       budget_global(:,:,period_inf) = 0.0_r8
       budget_counter(:,:,period_inf) = 0.0_r8
    elseif (trim(mode) == 'all') then
       budget_local(:,:,:) = 0.0_r8
       budget_global(:,:,:) = 0.0_r8
       budget_counter(:,:,period_inst) = 0.0_r8
       budget_counter(:,:,period_inst+1:) = 1.0_r8
    else
       call ESMF_LogWrite(trim(subname)//' mode '//trim(mode)//&
            ' not recognized', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    endif
  end subroutine med_diag_zero_mode

  !===============================================================================
  subroutine med_diag_zero_select(year, mon, day, tod)

    ! ------------------------------------------------------------------
    ! Zero out global budget diagnostic data.
    ! ------------------------------------------------------------------

    ! input/output variables
    integer, intent(in)  :: year
    integer, intent(in)  :: mon
    integer, intent(in)  :: day
    integer, intent(in)  :: tod

    ! local variables
    integer :: ip
    character(*), parameter :: subName = '(med_diag_zero_select) '
    ! ------------------------------------------------------------------

    do ip = 1,size(budget_diags%periods)
       if (ip == period_inst) then
          budget_local(:,:,ip) = 0.0_r8
          budget_global(:,:,ip) = 0.0_r8
          budget_counter(:,:,ip) = 0.0_r8
       endif
       if (ip==period_day .and. tod==0) then
          budget_local(:,:,ip) = 0.0_r8
          budget_global(:,:,ip) = 0.0_r8
          budget_counter(:,:,ip) = 0.0_r8
       endif
       if (ip==period_mon .and. day==1 .and. tod==0) then
          budget_local(:,:,ip) = 0.0_r8
          budget_global(:,:,ip) = 0.0_r8
          budget_counter(:,:,ip) = 0.0_r8
       endif
       if (ip==period_ann .and. mon==1 .and. day==1 .and. tod==0) then
          budget_local(:,:,ip) = 0.0_r8
          budget_global(:,:,ip) = 0.0_r8
          budget_counter(:,:,ip) = 0.0_r8
       endif
    enddo
  end subroutine med_diag_zero_select

  !===============================================================================
  subroutine med_phases_diag_accum(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Accumulate out global budget diagnostic data.
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer     :: ip
    character(*), parameter :: subName = '(med_diag_accum) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS
    do ip = period_inst+1,size(budget_diags%periods)
       budget_local(:,:,ip) = budget_local(:,:,ip) + budget_local(:,:,period_inst)
    enddo
    budget_counter(:,:,:) = budget_counter(:,:,:) + 1.0_r8

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_accum

  !===============================================================================
  subroutine med_diag_sum_main(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Sum local values to global on root
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: count
    integer       :: c_size  ! number of component send/recvs
    integer       :: f_size  ! number of fields
    integer       :: p_size  ! number of period types
    character(*), parameter :: subName = '(med_diag_sum_main) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    f_size = size(budget_diags%fields)
    c_size = size(budget_diags%comps)
    p_size = size(budget_diags%periods)

    count  = size(budget_global)
    budget_global_1d(:) = 0.0_r8

    call ESMF_VMReduce(vm, reshape(budget_local,(/count/)) , budget_global_1d, count, ESMF_REDUCE_SUM, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    budget_global = reshape(budget_global_1d,(/f_size,c_size,p_size/))

    budget_local(:,:,period_inst) = 0.0_r8

    call t_stopf('MED:'//subname)

  end subroutine med_diag_sum_main

  !===============================================================================
  subroutine med_phases_diag_atm(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global atm input/output flux diagnostics
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : compatm

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,nf,ip
    real(r8), pointer   :: afrac(:)
    real(r8), pointer   :: lfrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: areas(:)
    real(r8), pointer   :: lats(:)
    character(*), parameter :: subName = '(med_phases_diag_atm) '
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get fractions on atm mesh
    if (samegrid_atmlnd) then
      call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'lfrin', lfrac, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    areas => is_local%wrap%mesh_info(compatm)%areas
    lats  => is_local%wrap%mesh_info(compatm)%lats
    allocate(afrac(size(areas)))
    afrac = 1.0_R8

    !-------------------------------
    ! from atm to mediator (_recv suffix is what the mediator is receiving)
    !-------------------------------

    ip = period_inst

    do n = 1,size(afrac)
       nf = f_area
       budget_local(nf,c_atm_recv ,ip) = budget_local(nf,c_atm_recv ,ip) - areas(n)*afrac(n)
       budget_local(nf,c_lnd_arecv,ip) = budget_local(nf,c_lnd_arecv,ip) + areas(n)*lfrac(n)
       budget_local(nf,c_ocn_arecv,ip) = budget_local(nf,c_ocn_arecv,ip) + areas(n)*ofrac(n)
       if (is_local%wrap%mesh_info(compatm)%lats(n) > 0.0_r8) then
          budget_local(nf,c_inh_arecv,ip) = budget_local(nf,c_inh_arecv,ip) + areas(n)*ifrac(n)
       else
          budget_local(nf,c_ish_arecv,ip) = budget_local(nf,c_ish_arecv,ip) + areas(n)*ifrac(n)
       end if
    end do

    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', f_heat_swnet, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', f_heat_lwdn, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Note that passing f_watr_rain twice will just add up contributions from Faxa_rainc and Faxa_rainl
    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', f_watr_rain, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', f_watr_rain, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Note that passing f_watr_rain twice will just add up contributions from Faxa_snowc and Faxa_snowl
    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', f_watr_snow, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', f_watr_snow, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call diag_atm_wiso_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', &
            f_watr_rain_16O, f_watr_rain_18O, f_watr_rain_HDO, areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_atm_wiso_recv(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', &
            f_watr_rain_16O, f_watr_rain_18O, f_watr_rain_HDO, areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! heat implied by snow flux from atm to mediator
    budget_local(f_heat_latf,c_atm_recv ,ip) = -budget_local(f_watr_snow,c_atm_recv ,ip)*shr_const_latice
    budget_local(f_heat_latf,c_lnd_arecv,ip) = -budget_local(f_watr_snow,c_lnd_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_ocn_arecv,ip) = -budget_local(f_watr_snow,c_ocn_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_inh_arecv,ip) = -budget_local(f_watr_snow,c_inh_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_ish_arecv,ip) = -budget_local(f_watr_snow,c_ish_arecv,ip)*shr_const_latice

    !-------------------------------
    ! from mediator to atm
    !-------------------------------

    ip = period_inst

    do n = 1,size(afrac)
       budget_local(f_area,c_atm_send ,ip) = budget_local(f_area,c_atm_send ,ip) - areas(n)*afrac(n)
       budget_local(f_area,c_lnd_asend,ip) = budget_local(f_area,c_lnd_asend,ip) + areas(n)*lfrac(n)
       budget_local(f_area,c_ocn_asend,ip) = budget_local(f_area,c_ocn_asend,ip) + areas(n)*ofrac(n)
       if (is_local%wrap%mesh_info(compatm)%lats(n) > 0.0_r8) then
          budget_local(f_area,c_inh_asend,ip) = budget_local(f_area,c_inh_asend,ip) + areas(n)*ifrac(n)
       else
          budget_local(f_area,c_ish_asend,ip) = budget_local(f_area,c_ish_asend,ip) + areas(n)*ifrac(n)
       end if
    end do

    call diag_atm_send(is_local%wrap%FBExp(compatm), 'Faxx_lwup', f_heat_lwup, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_send(is_local%wrap%FBExp(compatm), 'Faxx_lat', f_heat_latvap, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_send(is_local%wrap%FBExp(compatm), 'Faxx_sen', f_heat_sen, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_atm_send(is_local%wrap%FBExp(compatm), 'Faxx_evap', f_watr_evap, &
         areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! water isotopes
    if (flds_wiso) then
       call diag_atm_wiso_send(is_local%wrap%FBImp(compatm,compatm), 'Faxa_evap_wiso', &
            f_watr_evap_16O, f_watr_evap_18O, f_watr_evap_HDO, &
            areas, lats, afrac, lfrac, ofrac, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    deallocate(afrac)
    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_atm

  subroutine diag_atm_recv(FB, fldname, nf, areas, lats, afrac, lfrac, ofrac, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: afrac(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(in)    :: ofrac(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1,size(data)
          budget(nf,c_atm_recv,ip)  = budget(nf,c_atm_recv,ip)  - areas(n)*data(n)*afrac(n)
          budget(nf,c_lnd_arecv,ip) = budget(nf,c_lnd_arecv,ip) + areas(n)*data(n)*lfrac(n)
          budget(nf,c_ocn_arecv,ip) = budget(nf,c_ocn_arecv,ip) + areas(n)*data(n)*ofrac(n)
          if (lats(n) > 0.0_r8) then
             budget(nf,c_inh_arecv,ip) = budget(nf,c_inh_arecv,ip) + areas(n)*data(n)*ifrac(n)
          else
             budget(nf,c_ish_arecv,ip) = budget(nf,c_ish_arecv,ip) + areas(n)*data(n)*ifrac(n)
          end if
       end do
    end if
  end subroutine diag_atm_recv

  subroutine diag_atm_send(FB, fldname, nf, areas, lats, afrac, lfrac, ofrac, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: afrac(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(in)    :: ofrac(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1,size(data)
          budget(nf,c_atm_send,ip)  = budget(nf,c_atm_send,ip)  - areas(n)*data(n)*afrac(n)
          budget(nf,c_lnd_asend,ip) = budget(nf,c_lnd_asend,ip) + areas(n)*data(n)*lfrac(n)
          budget(nf,c_ocn_asend,ip) = budget(nf,c_ocn_asend,ip) + areas(n)*data(n)*ofrac(n)
          if (lats(n) > 0.0_r8) then
             budget(nf,c_inh_asend,ip) = budget(nf,c_inh_asend,ip) + areas(n)*data(n)*ifrac(n)
          else
             budget(nf,c_ish_asend,ip) = budget(nf,c_ish_asend,ip) + areas(n)*data(n)*ifrac(n)
          end if
       end do
    end if
  end subroutine diag_atm_send

  subroutine diag_atm_wiso_recv(FB, fldname, nf_16O, nf_18O, nf_HDO, areas, lats, &
       afrac, lfrac, ofrac, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: afrac(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(in)    :: ofrac(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1,size(data, dim=2)
          budget(nf_16O,c_atm_recv,ip) = budget(nf_16O,c_atm_recv,ip) - areas(n)*afrac(n)*data(1,n)
          budget(nf_16O,c_lnd_arecv,ip) = budget(nf_16O,c_lnd_arecv,ip) + areas(n)*lfrac(n)*data(1,n)
          budget(nf_16O,c_ocn_arecv,ip) = budget(nf_16O,c_ocn_arecv,ip) + areas(n)*ofrac(n)*data(1,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_16O,c_inh_arecv,ip) = budget(nf_16O,c_inh_arecv,ip) + areas(n)*ifrac(n)*data(1,n)
          else
             budget(nf_16O,c_ish_arecv,ip) = budget(nf_16O,c_ish_arecv,ip) + areas(n)*ifrac(n)*data(1,n)
          end if

          budget(nf_18O,c_atm_recv,ip) = budget(nf_18O,c_atm_recv,ip) - areas(n)*afrac(n)*data(2,n)
          budget(nf_18O,c_lnd_arecv,ip) = budget(nf_18O,c_lnd_arecv,ip) + areas(n)*lfrac(n)*data(2,n)
          budget(nf_18O,c_ocn_arecv,ip) = budget(nf_18O,c_ocn_arecv,ip) + areas(n)*ofrac(n)*data(2,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_18O,c_inh_arecv,ip) = budget(nf_18O,c_inh_arecv,ip) + areas(n)*ifrac(n)*data(2,n)
          else
             budget(nf_18O,c_ish_arecv,ip) = budget(nf_18O,c_ish_arecv,ip) + areas(n)*ifrac(n)*data(2,n)
          end if

          budget(nf_HDO,c_atm_recv,ip) = budget(nf_HDO,c_atm_recv,ip) - areas(n)*afrac(n)*data(3,n)
          budget(nf_HDO,c_lnd_arecv,ip) = budget(nf_HDO,c_lnd_arecv,ip) + areas(n)*lfrac(n)*data(3,n)
          budget(nf_HDO,c_ocn_arecv,ip) = budget(nf_HDO,c_ocn_arecv,ip) + areas(n)*ofrac(n)*data(3,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_HDO,c_inh_arecv,ip) = budget(nf_HDO,c_inh_arecv,ip) + areas(n)*ifrac(n)*data(3,n)
          else
             budget(nf_HDO,c_ish_arecv,ip) = budget(nf_HDO,c_ish_arecv,ip) + areas(n)*ifrac(n)*data(3,n)
          end if
       end do
    end if
  end subroutine diag_atm_wiso_recv

  subroutine diag_atm_wiso_send(FB, fldname, nf_16O, nf_18O, nf_HDO, areas, lats, &
       afrac, lfrac, ofrac, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: afrac(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(in)    :: ofrac(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1,size(data, dim=2)
          budget(nf_16O,c_atm_send,ip) = budget(nf_16O,c_atm_send,ip) - areas(n)*afrac(n)*data(1,n)
          budget(nf_16O,c_lnd_asend,ip) = budget(nf_16O,c_lnd_asend,ip) + areas(n)*lfrac(n)*data(1,n)
          budget(nf_16O,c_ocn_asend,ip) = budget(nf_16O,c_ocn_asend,ip) + areas(n)*ofrac(n)*data(1,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_16O,c_inh_asend,ip) = budget(nf_16O,c_inh_asend,ip) + areas(n)*ifrac(n)*data(1,n)
          else
             budget(nf_16O,c_ish_asend,ip) = budget(nf_16O,c_ish_asend,ip) + areas(n)*ifrac(n)*data(1,n)
          end if

          budget(nf_18O,c_atm_send,ip) = budget(nf_18O,c_atm_send,ip) - areas(n)*afrac(n)*data(2,n)
          budget(nf_18O,c_lnd_asend,ip) = budget(nf_18O,c_lnd_asend,ip) + areas(n)*lfrac(n)*data(2,n)
          budget(nf_18O,c_ocn_asend,ip) = budget(nf_18O,c_ocn_asend,ip) + areas(n)*ofrac(n)*data(2,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_18O,c_inh_asend,ip) = budget(nf_18O,c_inh_asend,ip) + areas(n)*ifrac(n)*data(2,n)
          else
             budget(nf_18O,c_ish_asend,ip) = budget(nf_18O,c_ish_asend,ip) + areas(n)*ifrac(n)*data(2,n)
          end if

          budget(nf_HDO,c_atm_send,ip) = budget(nf_HDO,c_atm_send,ip) - areas(n)*afrac(n)*data(3,n)
          budget(nf_HDO,c_lnd_asend,ip) = budget(nf_HDO,c_lnd_asend,ip) + areas(n)*lfrac(n)*data(3,n)
          budget(nf_HDO,c_ocn_asend,ip) = budget(nf_HDO,c_ocn_asend,ip) + areas(n)*ofrac(n)*data(3,n)
          if (lats(n) > 0.0_r8) then
             budget(nf_HDO,c_inh_asend,ip) = budget(nf_HDO,c_inh_asend,ip) + areas(n)*ifrac(n)*data(3,n)
          else
             budget(nf_HDO,c_ish_asend,ip) = budget(nf_HDO,c_ish_asend,ip) + areas(n)*ifrac(n)*data(3,n)
          end if
       end do
    end if
  end subroutine diag_atm_wiso_send

  !===============================================================================
  subroutine med_phases_diag_lnd( gcomp, rc)

   ! ------------------------------------------------------------------
    ! Compute global lnd input/output flux diagnostics
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : complnd

    ! intput/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    real(r8), pointer   :: lfrac(:)
    integer             :: n,ip, ic
    real(r8), pointer   :: areas(:)
    character(*), parameter :: subName = '(med_phases_diag_lnd) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get fractions on lnd mesh
    call fldbun_getdata1d(is_local%wrap%FBfrac(complnd), 'lfrin', lfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    areas => is_local%wrap%mesh_info(complnd)%areas

    !-------------------------------
    ! from land to mediator
    !-------------------------------

    ic = c_lnd_recv
    ip = period_inst
    do n = 1, size(lfrac)
       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + areas(n)*lfrac(n)
    end do

    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Fall_swnet', f_heat_swnet  , ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup' , f_heat_lwup   , ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat'  , f_heat_latvap , ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen'  , f_heat_sen    , ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap' , f_watr_evap   , ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsur', f_watr_roff, ic,&
         areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofgwl', f_watr_roff, ic,&
         areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsub', f_watr_roff, ic,&
         areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Flrl_irrig' , f_watr_roff, ic,&
         areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi'  , f_watr_ioff, ic,&
         areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call diag_lnd_wiso(is_local%wrap%FBImp(complnd,complnd), 'Flrl_evap_wiso', &
            f_watr_evap_16O, f_watr_evap_18O, f_watr_evap_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofl_wiso', &
            f_watr_roff_16O, f_watr_roff_18O, f_watr_roff_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi_wiso', &
            f_watr_ioff_16O, f_watr_ioff_18O, f_watr_ioff_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! to land from mediator
    !-------------------------------

    ic = c_lnd_send
    ip = period_inst

    do n = 1,size(lfrac)
       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + areas(n)*lfrac(n)
    end do
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Faxa_lwdn' , f_heat_lwdn, ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Faxa_rainc', f_watr_rain, ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Faxa_rainl', f_watr_rain, ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Faxa_snowc', f_watr_snow, ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Faxa_snowl', f_watr_snow, ic, areas, lfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_lnd(is_local%wrap%FBExp(complnd), 'Flrl_flood', f_watr_roff, ic, areas, lfrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call diag_lnd_wiso(is_local%wrap%FBExp(complnd), 'Faxa_rainc_wiso', &
            f_watr_rain_16O, f_watr_rain_18O, f_watr_rain_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBExp(complnd), 'Faxa_rainl_wiso', &
            f_watr_rain_16O, f_watr_rain_18O, f_watr_rain_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBExp(complnd), 'Faxa_snowc_wiso', &
            f_watr_snow_16O, f_watr_snow_18O, f_watr_snow_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBExp(complnd), 'Faxa_snowl_wiso', &
            f_watr_snow_16O, f_watr_snow_18O, f_watr_snow_HDO, ic, areas, lfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_lnd_wiso(is_local%wrap%FBExp(complnd), 'Flrl_flood_wiso', &
            f_watr_roff_16O, f_watr_roff_18O, f_watr_roff_HDO, ic, areas, lfrac, budget_local, minus=.true., rc=rc)
    end if

    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_lnd

  subroutine diag_lnd(FB, fldname, nf, ic, areas, lfrac, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data)
          if (present(minus)) then
             budget(nf,ic,ip) = budget(nf,ic,ip) - areas(n)*lfrac(n)*data(n)
          else
             budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*lfrac(n)*data(n)
          end if
       end do
    end if
  end subroutine diag_lnd

  subroutine diag_lnd_wiso(FB, fldname, nf_16O, nf_18O, nf_HDO, ic, areas, lfrac, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lfrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data, dim=2)
          if (present(minus)) then
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) - areas(n)*lfrac(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) - areas(n)*lfrac(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) - areas(n)*lfrac(n)*data(3,n)
          else
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) + areas(n)*lfrac(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) + areas(n)*lfrac(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) + areas(n)*lfrac(n)*data(3,n)
          end if
       end do
    end if
  end subroutine diag_lnd_wiso

  !===============================================================================
  subroutine med_phases_diag_rof( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global river input/output
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : comprof

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: ic, ip
    real(r8), pointer   :: areas(:)
    character(*), parameter :: subName = '(med_phases_diag_rof) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    areas => is_local%wrap%mesh_info(comprof)%areas

    !-------------------------------
    ! from river to mediator
    !-------------------------------

    ic = c_rof_recv
    ip = period_inst

    call diag_rof(is_local%wrap%FBImp(comprof,comprof), 'Flrr_flood', f_watr_roff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofl' , f_watr_roff, ic, areas, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofi' , f_watr_ioff, ic, areas, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldbun_fldchk(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofl_glc', rc=rc)) then
      call diag_rof(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofi_glc' , f_watr_roff, ic, areas, budget_local, minus=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if ( fldbun_fldchk(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofi_glc', rc=rc)) then
      call diag_rof(is_local%wrap%FBImp(comprof,comprof), 'Forr_rofi_glc' , f_watr_ioff, ic, areas, budget_local, minus=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call diag_rof_wiso(is_local%wrap%FBExp(comprof), 'Forr_flood_wiso', &
            f_watr_ioff_16O, f_watr_ioff_18O, f_watr_ioff_HDO, ic, areas, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_rof_wiso(is_local%wrap%FBExp(comprof), 'Forr_rofl_wiso', &
            f_watr_roff_16O, f_watr_roff_18O, f_watr_roff_HDO, ic, areas, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_rof_wiso(is_local%wrap%FBExp(comprof), 'Forr_rofi_wiso', &
            f_watr_ioff_16O, f_watr_ioff_18O, f_watr_ioff_HDO, ic, areas, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! to river from mediator
    !-------------------------------

    ic = c_rof_send
    ip = period_inst

    call diag_rof(is_local%wrap%FBExp(comprof), 'Flrl_rofsur', f_watr_roff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBExp(comprof), 'Flrl_rofgwl', f_watr_roff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBExp(comprof), 'Flrl_rofsub', f_watr_roff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBExp(comprof), 'Flrl_irrig' , f_watr_roff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_rof(is_local%wrap%FBExp(comprof), 'Flrl_rofi'  , f_watr_ioff, ic, areas, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (fldbun_fldchk(is_local%wrap%FBExp(comprof), 'Fgrg_rofl', rc=rc)) then
      call diag_rof(is_local%wrap%FBExp(comprof), 'Fgrg_rofl'  , f_watr_roff, ic, areas, budget_local, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (fldbun_fldchk(is_local%wrap%FBExp(comprof), 'Fgrg_rofi', rc=rc)) then
      call diag_rof(is_local%wrap%FBExp(comprof), 'Fgrg_rofi'  , f_watr_ioff, ic, areas, budget_local, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call diag_rof_wiso(is_local%wrap%FBExp(comprof), 'Flrl_rofl_wiso', &
            f_watr_roff_16O, f_watr_roff_18O, f_watr_roff_HDO, ic, areas, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_rof_wiso(is_local%wrap%FBExp(comprof), 'Flrl_rofi_wiso', &
         f_watr_ioff_16O, f_watr_ioff_18O, f_watr_ioff_HDO, ic, areas, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_rof

  subroutine diag_rof(FB, fldname, nf, ic, areas, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data)
          if (present(minus)) then
             budget(nf,ic,ip) = budget(nf,ic,ip) - areas(n)*data(n)
          else
             budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*data(n)
          end if
       end do
    end if
  end subroutine diag_rof

  subroutine diag_rof_wiso(FB, fldname, nf_16O, nf_18O, nf_HDO, ic, areas, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data, dim=2)
          if (present(minus)) then
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) - areas(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) - areas(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) - areas(n)*data(3,n)
          else
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) + areas(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) + areas(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) + areas(n)*data(3,n)
          end if
       end do
    end if
  end subroutine diag_rof_wiso

  !===============================================================================
  subroutine med_phases_diag_glc( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global glc output
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : compglc

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: ic, ip, ns
    real(r8), pointer   :: areas(:)
    character(*), parameter :: subName = '(med_phases_diag_glc) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    !-------------------------------
    ! from glc to mediator
    !-------------------------------

    ! TODO: this will not be correct if there is more than 1 ice sheet
    ic = c_glc_recv
    ip = period_inst

    do ns = 1,is_local%wrap%num_icesheets
       areas => is_local%wrap%mesh_info(compglc(ns))%areas
       call diag_glc(is_local%wrap%FBImp(compglc(ns),compglc(ns)), 'Fgrg_rofl', f_watr_roff, ic, areas, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_glc(is_local%wrap%FBImp(compglc(ns),compglc(ns)), 'Fgrg_rofi', f_watr_ioff, ic, areas, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_glc(is_local%wrap%FBImp(compglc(ns),compglc(ns)), 'Figg_rofi', f_watr_ioff, ic, areas, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_glc

  subroutine diag_glc(FB, fldname, nf, ic, areas, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data)
          if (present(minus)) then
             budget(nf,ic,ip) = budget(nf,ic,ip) - areas(n)*data(n)
          else
             budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*data(n)
          end if
       end do
    end if
  end subroutine diag_glc

  !===============================================================================
  subroutine med_phases_diag_ocn( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ocn input from mediator
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : compocn, compatm

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8)            :: wgt_i,wgt_o
    real(r8), pointer   :: ifrac(:)  ! ice fraction in ocean grid cell
    real(r8), pointer   :: ofrac(:)  ! non-ice fraction nin ocean grid cell
    real(r8), pointer   :: sfrac(:)  ! sum of ifrac and ofrac
    real(r8), pointer   :: areas(:)
    real(r8), pointer   :: data(:)
    character(*), parameter :: subName = '(med_phases_diag_ocn) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbun_getdata1d(is_local%wrap%FBfrac(compocn), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(sfrac(size(ofrac)))
    sfrac(:) = 1._r8

    !areas => is_local%wrap%mesh_info(compocn)%areas
    call fldbun_getdata1d(is_local%wrap%FBarea(compocn), 'area', areas, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from ocn to mediator
    !-------------------------------

    ip = period_inst
    ic = c_ocn_recv
    do n = 1,size(ofrac)
       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + areas(n)*ofrac(n)
    end do

    if ( fldbun_fldchk(is_local%wrap%FBImp(compocn,compocn), 'Fioo_q', rc=rc)) then
       call fldbun_getdata1d(is_local%wrap%FBImp(compocn,compocn), 'Fioo_q', data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(ifrac)
          wgt_o = areas(n) * ofrac(n)
          wgt_i = areas(n) * ifrac(n)
          budget_local(f_heat_frz,ic,ip) = budget_local(f_heat_frz,ic,ip) + (wgt_o + wgt_i)*max(0.0_r8,data(n))
       end do
    end if

    if (f_watr_frz /= unset_index) then
       budget_local(f_watr_frz,ic,ip) = budget_local(f_heat_frz,ic,ip) * HFLXtoWFLX
    end if

    !-------------------------------
    ! from mediator to ocn
    !-------------------------------

    ic = c_ocn_send
    ip = period_inst

    do n = 1,size(ofrac)
       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + areas(n)*ofrac(n)
    end do

    if (fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lwnet', rc=rc)) then ! MOM6
       call diag_ocn(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup', f_heat_lwup, ic, areas, ofrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn(is_local%wrap%FBImp(compatm,compocn), 'Faxa_lwdn', f_heat_lwdn, ic, areas, ofrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else ! POP
       call diag_ocn(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup' , f_heat_lwup   , ic, areas, ofrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Faxa_lwdn' , f_heat_lwdn   , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call diag_ocn(is_local%wrap%FBMed_aoflux_o, 'Faox_sen' , f_heat_sen    , ic, areas, ofrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', f_watr_evap   , ic, areas, ofrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lat', rc=rc)) then ! POP
       call diag_ocn(is_local%wrap%FBMed_aoflux_o, 'Faox_lat'  , f_heat_latvap , ic, areas, ofrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else ! MOM6
       call diag_ocn(is_local%wrap%FBMed_aoflux_o, 'Faox_evap' , f_heat_latvap , ic, areas, ofrac, budget_local, &
            scale=shr_const_latvap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if


    call diag_ocn(is_local%wrap%FBExp(compocn), 'Fioi_meltw', f_watr_melt   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Fioi_bergw', f_watr_melt   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Fioi_melth', f_heat_melt   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Fioi_bergh', f_heat_melt   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call diag_ocn(is_local%wrap%FBExp(compocn), 'Fioi_salt' , f_watr_salt   , ic, areas, sfrac, budget_local, &
         scale=SFLXtoWFLX, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_swnet', f_heat_swnet  , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc) .and. &
             fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', rc=rc) .and. &
             fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', rc=rc) .and. &
             fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', rc=rc)) then
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', f_heat_swnet  , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', f_heat_swnet  , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', f_heat_swnet  , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', f_heat_swnet  , ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call diag_ocn(is_local%wrap%FBExp(compocn), 'Faxa_rain' , f_watr_rain   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Faxa_snow' , f_watr_snow   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , f_watr_roff   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , f_watr_ioff   , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc' , rc=rc)) then
      call diag_ocn(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc' , f_watr_roff   , ic, areas, sfrac, budget_local, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if ( fldbun_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc' , rc=rc)) then
      call diag_ocn(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc' , f_watr_ioff   , ic, areas, sfrac, budget_local, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call diag_ocn_wiso(is_local%wrap%FBMed_aoflux_o, 'Faox_evap_wiso', &
            f_watr_evap_16O, f_watr_evap_18O, f_watr_evap_HDO, ic, areas, ofrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn_wiso(is_local%wrap%FBExp(compocn), 'Fioi_meltw_wiso', &
            f_watr_melt_16O, f_watr_melt_HDO, f_watr_melt_HDO, ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn_wiso(is_local%wrap%FBExp(compocn), 'Fioi_rain_wiso' , &
            f_watr_rain_16O, f_watr_rain_HDO, f_watr_rain_HDO, ic, areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn_wiso(is_local%wrap%FBExp(compocn), 'Fioi_snow_wiso' , &
            f_watr_snow_16O, f_watr_snow_HDO, f_watr_snow_HDO, ic,  areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn_wiso(is_local%wrap%FBExp(compocn), 'Foxx_rofl_wiso' , &
            f_watr_roff_16O, f_watr_roff_HDO, f_watr_roff_HDO, ic,  areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ocn_wiso(is_local%wrap%FBExp(compocn), 'Foxx_rofi_wiso' , &
            f_watr_ioff_16O, f_watr_ioff_HDO, f_watr_ioff_HDO, ic,  areas, sfrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hrain', f_heat_rain , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', f_heat_snow , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hevap', f_heat_evap , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hcond', f_heat_cond , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', f_heat_rofl , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ocn(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', f_heat_rofi , ic, areas, sfrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    deallocate(sfrac)
    call t_stopf('MED:'//subname)

  end subroutine med_phases_diag_ocn

  subroutine diag_ocn(FB, fldname, nf, ic, areas, frac, budget, scale, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: frac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    real(r8), optional     , intent(in)    :: scale
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data)
          if (present(scale)) then
             budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*frac(n)*data(n)*scale
          else
             budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*frac(n)*data(n)
          end if
       end do
    end if
  end subroutine diag_ocn

  subroutine diag_ocn_wiso(FB, fldname, nf_16O, nf_18O, nf_HDO, ic, areas, frac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    integer                , intent(in)    :: ic
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: frac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data, dim=2)
          budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) + areas(n)*frac(n)*data(1,n)
          budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) + areas(n)*frac(n)*data(2,n)
          budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) + areas(n)*frac(n)*data(3,n)
       end do
    end if
  end subroutine diag_ocn_wiso

  !===============================================================================
  subroutine med_phases_diag_ice_ice2med( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ice input/output flux diagnostics
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : compice

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: areas(:)
    real(r8), pointer   :: lats(:)
    character(*), parameter :: subName = '(med_phases_diag_ice_ice2med) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbun_getdata1d(is_local%wrap%FBFrac(compice), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(is_local%wrap%FBFrac(compice), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    areas => is_local%wrap%mesh_info(compice)%areas
    lats  => is_local%wrap%mesh_info(compice)%lats

    ip = period_inst

    do n = 1,size(ifrac)
       if (lats(n) > 0.0_r8) then
          ic = c_inh_recv
       else
          ic = c_ish_recv
       endif
       budget_local(f_area ,ic,ip) = budget_local(f_area ,ic,ip) + areas(n)*ifrac(n)
    end do

    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_melth', f_heat_melt, &
         areas, lats, ifrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_meltw', f_watr_melt, &
         areas, lats, ifrac, budget_local, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_salt', f_watr_salt, &
         areas, lats, ifrac, budget_local, minus=.true., scale=SFLXtoWFLX, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldbun_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc) .and. &
         fldbun_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdf', rc=rc) .and. &
         fldbun_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idr', rc=rc) .and. &
         fldbun_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idf', rc=rc)) then
       call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', f_heat_swnet, &
            areas, lats, ifrac, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdf', f_heat_swnet, &
            areas, lats, ifrac, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idr', f_heat_swnet, &
            areas, lats, ifrac, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idf', f_heat_swnet, &
            areas, lats, ifrac, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen', f_heat_swnet, &
            areas, lats, ifrac, budget_local, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Faii_swnet', f_heat_swnet, &
         areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', f_heat_lwup, &
         areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Faii_lat', f_heat_latvap, &
         areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Faii_sen', f_heat_sen, &
         areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_recv(is_local%wrap%FBImp(compice,compice), 'Faii_evap', f_watr_evap, &
         areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call diag_ice_recv_wiso(is_local%wrap%FBImp(compice,compice), 'Fioi_meltw_wiso', &
            f_watr_melt_16O, f_watr_melt_18O, f_watr_melt_HDO, areas, lats, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ice_recv_wiso(is_local%wrap%FBImp(compice,compice), 'Faii_evap_wiso', &
            f_watr_evap_16O, f_watr_evap_18O, f_watr_evap_HDO, areas, lats, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_ice_ice2med

  subroutine diag_ice_recv(FB, fldname, nf, areas, lats, ifrac, budget, minus, scale, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical,  optional     , intent(in)    :: minus
    real(r8), optional     , intent(in)    :: scale
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ic, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1,size(data)
          if (lats(n) > 0.0_r8) then
             ic = c_inh_recv
          else
             ic = c_ish_recv
          endif
          if (present(minus)) then
             if (present(scale)) then
                budget(nf ,ic,ip) = budget(nf ,ic,ip) - areas(n)*ifrac(n)*data(n)*scale
             else
                budget(nf ,ic,ip) = budget(nf ,ic,ip) - areas(n)*ifrac(n)*data(n)
             end if
          else
             if (present(scale)) then
                budget(nf ,ic,ip) = budget(nf ,ic,ip) + areas(n)*ifrac(n)*data(n)*scale
             else
                budget(nf ,ic,ip) = budget(nf ,ic,ip) + areas(n)*ifrac(n)*data(n)
             end if
          end if
       end do
    end if
  end subroutine diag_ice_recv

  subroutine diag_ice_recv_wiso(FB, fldname, nf_16O, nf_18O, nf_HDO, areas, lats, ifrac, budget, minus, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    logical, optional      , intent(in)    :: minus
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ic, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data, dim=2)
          if (lats(n) > 0.0_r8) then
             ic = c_inh_recv
          else
             ic = c_ish_recv
          endif
          if (present(minus)) then
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) - areas(n)*ifrac(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) - areas(n)*ifrac(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) - areas(n)*ifrac(n)*data(3,n)
          else
             budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) + areas(n)*ifrac(n)*data(1,n)
             budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) + areas(n)*ifrac(n)*data(2,n)
             budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) + areas(n)*ifrac(n)*data(3,n)
          end if
       end do
    end if
  end subroutine diag_ice_recv_wiso

  !===============================================================================
  subroutine med_phases_diag_ice_med2ice( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ice input/output flux diagnostics
    ! ------------------------------------------------------------------

    use med_internalstate_mod, only : compice

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8)            :: wgt_i, wgt_o
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: data(:)
    real(r8), pointer   :: areas(:)
    real(r8), pointer   :: lats(:)
    character(*), parameter :: subName = '(med_phases_diag_ice_med2ice) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)
    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbun_getdata1d(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getdata1d(is_local%wrap%FBfrac(compice), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    areas => is_local%wrap%mesh_info(compice)%areas
    lats  => is_local%wrap%mesh_info(compice)%lats

    ip = period_inst

    do n = 1,size(ifrac)
       if (lats(n) > 0.0_r8) then
          ic = c_inh_send
       else
          ic = c_ish_send
       endif
       budget_local(f_area ,ic,ip) = budget_local(f_area ,ic,ip) + areas(n)*ifrac(n)
    end do

    call diag_ice_send(is_local%wrap%FBExp(compice), 'Faxa_lwdn', f_heat_lwdn, areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_send(is_local%wrap%FBExp(compice), 'Faxa_rain', f_watr_rain, areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call diag_ice_send(is_local%wrap%FBExp(compice), 'Faxa_snow', f_watr_snow, areas, lats, ifrac, budget_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldbun_fldchk(is_local%wrap%FBExp(compice), 'Fioo_q', rc=rc)) then
       call fldbun_getdata1d(is_local%wrap%FBExp(compice), 'Fioo_q', data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(data)
          wgt_o = areas(n) * ofrac(n)
          wgt_i = areas(n) * ifrac(n)
          if (lats(n) > 0.0_r8) then
             ic = c_inh_send
          else
             ic = c_ish_send
          endif
          budget_local(f_heat_frz,ic,ip) = budget_local(f_heat_frz,ic,ip) - (wgt_o + wgt_i)*max(0.0_r8,data(n))
       end do
    end if

    ic = c_inh_send
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    if (trim(budget_table_version) == 'v0') then
       budget_local(f_watr_frz ,ic,ip) =  budget_local(f_heat_frz ,ic,ip)*HFLXtoWFLX
    end if

    ic = c_ish_send
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    if (trim(budget_table_version) == 'v0') then
       budget_local(f_watr_frz ,ic,ip) =  budget_local(f_heat_frz ,ic,ip)*HFLXtoWFLX
    end if

    if (flds_wiso) then
       call diag_ice_send_wiso(is_local%wrap%FBExp(compice), 'Faxa_rain_wiso', &
            f_watr_rain_16O, f_watr_rain_18O, f_watr_rain_HDO, areas, lats, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call diag_ice_send_wiso(is_local%wrap%FBExp(compice), 'Faxa_snow_wiso', &
            f_watr_snow_16O, f_watr_snow_18O, f_watr_snow_HDO, areas, lats, ifrac, budget_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)
  end subroutine med_phases_diag_ice_med2ice

  subroutine diag_ice_send(FB, fldname, nf, areas, lats, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc
    ! local variables
    integer           :: n, ic, ip
    real(r8), pointer :: data(:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    ip = period_inst
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata1d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(data)
          if (lats(n) > 0.0_r8) then
             ic = c_inh_send
          else
             ic = c_ish_send
          endif
          budget(nf,ic,ip) = budget(nf,ic,ip) + areas(n)*ifrac(n)*data(n)
       end do
    end if
  end subroutine diag_ice_send

  subroutine diag_ice_send_wiso(FB, fldname, nf_16O, nf_18O, nf_HDO, areas, lats, ifrac, budget, rc)
    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: fldname
    integer                , intent(in)    :: nf_16O
    integer                , intent(in)    :: nf_18O
    integer                , intent(in)    :: nf_HDO
    real(r8)               , intent(in)    :: areas(:)
    real(r8)               , intent(in)    :: lats(:)
    real(r8)               , intent(in)    :: ifrac(:)
    real(r8)               , intent(inout) :: budget(:,:,:)
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n, ic, ip
    real(r8), pointer :: data(:,:)
    ! ------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if ( fldbun_fldchk(FB, trim(fldname), rc=rc)) then
       call fldbun_getdata2d(FB, trim(fldname), data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ip = period_inst
       do n = 1, size(data, dim=2)
          if (lats(n) > 0.0_r8) then
             ic = c_inh_send
          else
             ic = c_ish_send
          endif
          budget(nf_16O,ic,ip) = budget(nf_16O,ic,ip) + areas(n)*ifrac(n)*data(1,n)
          budget(nf_18O,ic,ip) = budget(nf_18O,ic,ip) + areas(n)*ifrac(n)*data(2,n)
          budget(nf_HDO,ic,ip) = budget(nf_HDO,ic,ip) + areas(n)*ifrac(n)*data(3,n)
       end do
    end if
  end subroutine diag_ice_send_wiso

  !===============================================================================
  subroutine med_phases_diag_print(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Print global budget diagnostics.
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    ! local variables
    type(ESMF_Clock)      :: clock
    type(ESMF_Alarm)      :: stop_alarm
    type(ESMF_Time)       :: nextTime
    integer               :: date        ! coded date, seconds
    integer               :: year
    integer               :: mon
    integer               :: day
    integer               :: tod
    integer               :: output_level ! print level
    logical               :: sumdone      ! has a sum been computed yet
    integer               :: ip
    integer               :: c_size       ! number of component send/recvs
    integer               :: f_size       ! number of fields
    integer               :: p_size       ! number of period types
    real(r8), allocatable :: datagpr(:,:,:)
    logical, save         :: firstcall = .true.
#ifdef DEBUG
    character(len=CL)     :: timestr
#endif
    character(*), parameter :: subName = '(med_phases_diag_print) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------------
    ! Print budget data if appropriate
    !-------------------------------------------------------------------------------

    ! Get clock and alarm info
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! NOTE - we are using the next time to ensure that budgets are
    ! written at the end of the run correctly This duplicates the
    ! behavior in the restart and history file output in that the time
    ! stamp is the next time and not the actual current time
    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( nextTime, yy=year, mm=mon, dd=day, s=tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    date = year*10000 + mon*100 + day

#ifdef DEBUG
    if(maintask) then
       write(timestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') year,'-',mon,'-',day,'-',tod
       write(logunit,' (a)') trim(subname)//": time = "//trim(timestr)
    endif
#endif

    if (firstcall) then
       firstcall = .false.
       return
    endif

    sumdone = .false.
    do ip = 1,size(budget_diags%periods)

       ! Determine output level for this period type
       output_level = 0
       if (ip == period_inst) then
          output_level = max(output_level, budget_print_inst)
       end if
       if (ip == period_day .and. tod == 0) then
          output_level = max(output_level, budget_print_daily)
       end if
       if (ip == period_mon .and. day == 1 .and. tod == 0) then
          output_level = max(output_level, budget_print_month)
       end if
       if (ip == period_ann .and. mon == 1 .and. day == 1 .and. tod == 0) then
          output_level = max(output_level, budget_print_ann)
       end if
       if (ip == period_inf .and. mon == 1 .and. day == 1 .and. tod == 0) then
          output_level = max(output_level, budget_print_ltann)
       end if
       if (ip == period_inf) then
          call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=stop_alarm, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (ESMF_AlarmIsRinging(stop_alarm, rc=rc)) then
             output_level = max(output_level, budget_print_ltend)
             call ESMF_AlarmRingerOff( stop_alarm, rc=rc )
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif

       ! Currently output_level is limited to levels of 0,1,2, 3
       ! (see comment for print options at top)

       if (output_level > 0) then
          if (.not. sumdone) then
             ! Some budgets will be printed for this period type
             ! Determine sums if not already done
             call med_diag_sum_main(gcomp, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             sumdone = .true.
          end if

          if (maintask) then
             c_size = size(budget_diags%comps)
             f_size = size(budget_diags%fields)
             p_size = size(budget_diags%periods)
             allocate(datagpr(f_size, c_size, p_size))
             datagpr(:,:,:) = budget_global(:,:,:)

             ! budget normalizations (global area and 1e6 for water)
             datagpr = datagpr/(4.0_r8*shr_const_pi)
             datagpr(f_watr_beg:f_watr_end,:,:) = datagpr(f_watr_beg:f_watr_end,:,:) * 1.0e6_r8
             if ( flds_wiso ) then
                datagpr(iso0(1):isof(nisotopes),:,:) = datagpr(iso0(1):isof(nisotopes),:,:) * 1.0e6_r8
             end if
             datagpr(:,:,:) = datagpr(:,:,:)/budget_counter(:,:,:)

             ! Write diagnostic tables to logunit (maintask only)
             if (output_level >= 3) then
                ! detail atm budgets and breakdown into components ---
                call med_diag_print_atm(datagpr, ip, date, tod)
             end if
             if (output_level >= 2) then
                ! detail lnd/ocn/ice component budgets ----
                call med_diag_print_lnd_ice_ocn(datagpr, ip, date, tod)
             end if
             if (output_level >= 1) then
                ! net summary budgets
                call med_diag_print_summary(datagpr, ip, date, tod)
             endif
             write(diagunit,*) ' '

             deallocate(datagpr)

          endif ! output_level > 0 and maintask
       end if ! if maintask
    enddo  ! ip = 1, period_types

    !-------------------------------------------------------------------------------
    ! Zero budget data
    !-------------------------------------------------------------------------------

    call med_diag_zero(year, mon, day, tod)

  end subroutine med_phases_diag_print

  !===============================================================================
  subroutine med_diag_print_atm(data, ip, date, tod)

    ! ---------------------------------------------------------
    ! detail atm budgets and breakdown into components
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such
    integer , intent(in) :: ip          ! period index
    integer , intent(in) :: date
    integer , intent(in) :: tod

    ! local variables
    integer           :: ic,nf,is ! data array indicies
    integer           :: ica,icl
    integer           :: icn,ics,ico
    character(len=40) :: str         ! string
    character(*), parameter:: subName = '(med_phases_diag_print_atm) '
    ! ------------------------------------------------------------------

    ica = 0
    icl = 0
    icn = 0
    ics = 0
    ico = 0
    str = ""
    do ic = 1,2
       if (ic == 1) then    ! from atm to mediator
          ica = c_atm_recv ! total from atm
          icl = c_lnd_arecv ! from land   to med on atm grid
          icn = c_inh_arecv ! from ice-nh to med on atm grid
          ics = c_ish_arecv ! from ice-sh to med on atm grid
          ico = c_ocn_arecv ! from ocn to to med on atm grid
          str = "ATM_to_CPL"
       elseif (ic == 2) then ! from mediator to atm
          ica = c_atm_send  ! merged to atm
          icl = c_lnd_asend  ! from land   to atm
          icn = c_inh_asend  ! from ice-nh to atm
          ics = c_ish_asend  ! from ice-sh to atm
          ico = c_ocn_asend  ! from ocn    to atm
          str = "CPL_TO_ATM"
       endif

       write(diagunit,*) ' '
       write(diagunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ', &
            trim(budget_diags%periods(ip)%name), ': date = ', date, tod
       write(diagunit,FA0) &
            budget_diags%comps(ica)%name,&
            budget_diags%comps(icl)%name,&
            budget_diags%comps(icn)%name,&
            budget_diags%comps(ics)%name,&
            budget_diags%comps(ico)%name,' *SUM*  '
       write(diagunit,FA1) budget_diags%fields(f_area)%name,&
            data(f_area,ica,ip), &
            data(f_area,icl,ip), &
            data(f_area,icn,ip), &
            data(f_area,ics,ip), &
            data(f_area,ico,ip), &
            data(f_area,ica,ip) + data(f_area,icl,ip) + &
            data(f_area,icn,ip) + data(f_area,ics,ip) + data(f_area,ico,ip)

       write(diagunit,*) ' '
       write(diagunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',&
            trim(budget_diags%periods(ip)%name),': date = ',date,tod
       write(diagunit,FA0) &
            budget_diags%comps(ica)%name,&
            budget_diags%comps(icl)%name,&
            budget_diags%comps(icn)%name,&
            budget_diags%comps(ics)%name,&
            budget_diags%comps(ico)%name,' *SUM*  '
       do nf = f_heat_beg, f_heat_end
          write(diagunit,FA1) budget_diags%fields(nf)%name,&
               data(nf,ica,ip), &
               data(nf,icl,ip), &
               data(nf,icn,ip), &
               data(nf,ics,ip), &
               data(nf,ico,ip), &
               data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
       enddo
       write(diagunit,FA1)    '   *SUM*'   ,&
            sum(data(f_heat_beg:f_heat_end,ica,ip)), &
            sum(data(f_heat_beg:f_heat_end,icl,ip)), &
            sum(data(f_heat_beg:f_heat_end,icn,ip)), &
            sum(data(f_heat_beg:f_heat_end,ics,ip)), &
            sum(data(f_heat_beg:f_heat_end,ico,ip)), &
            sum(data(f_heat_beg:f_heat_end,ica,ip)) + sum(data(f_heat_beg:f_heat_end,icl,ip)) + &
            sum(data(f_heat_beg:f_heat_end,icn,ip)) + sum(data(f_heat_beg:f_heat_end,ics,ip)) + &
            sum(data(f_heat_beg:f_heat_end,ico,ip))

       write(diagunit,*) ' '
       write(diagunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',&
            trim(budget_diags%periods(ip)%name),': date = ',date,tod
       write(diagunit,FA0) &
            budget_diags%comps(ica)%name,&
            budget_diags%comps(icl)%name,&
            budget_diags%comps(icn)%name,&
            budget_diags%comps(ics)%name,&
            budget_diags%comps(ico)%name,' *SUM*  '
       do nf = f_watr_beg, f_watr_end
          write(diagunit,FA1) budget_diags%fields(nf)%name,&
               data(nf,ica,ip), &
               data(nf,icl,ip), &
               data(nf,icn,ip), &
               data(nf,ics,ip), &
               data(nf,ico,ip), &
               data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
       enddo
       write(diagunit,FA1)    '   *SUM*'   ,&
            sum(data(f_watr_beg:f_watr_end,ica,ip)), &
            sum(data(f_watr_beg:f_watr_end,icl,ip)), &
            sum(data(f_watr_beg:f_watr_end,icn,ip)), &
            sum(data(f_watr_beg:f_watr_end,ics,ip)), &
            sum(data(f_watr_beg:f_watr_end,ico,ip)), &
            sum(data(f_watr_beg:f_watr_end,ica,ip)) + sum(data(f_watr_beg:f_watr_end,icl,ip)) + &
            sum(data(f_watr_beg:f_watr_end,icn,ip)) + sum(data(f_watr_beg:f_watr_end,ics,ip)) + &
            sum(data(f_watr_beg:f_watr_end,ico,ip))

       if ( flds_wiso ) then
          do is = 1, nisotopes
             write(diagunit,*) ' '
             write(diagunit,FAH) subname,trim(str)//' '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                  trim(budget_diags%periods(ip)%name),': date = ',date,tod
             write(diagunit,FA0) &
                  budget_diags%comps(ica)%name,&
                  budget_diags%comps(icl)%name,&
                  budget_diags%comps(icn)%name,&
                  budget_diags%comps(ics)%name,&
                  budget_diags%comps(ico)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(diagunit,FA1) budget_diags%fields(nf)%name,&
                     data(nf,ica,ip), &
                     data(nf,icl,ip), &
                     data(nf,icn,ip), &
                     data(nf,ics,ip), &
                     data(nf,ico,ip), &
                     data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
             enddo
             write(diagunit,FA1)    '   *SUM*', &
                  sum(data(iso0(is):isof(is),ica,ip)), &
                  sum(data(iso0(is):isof(is),icl,ip)), &
                  sum(data(iso0(is):isof(is),icn,ip)), &
                  sum(data(iso0(is):isof(is),ics,ip)), &
                  sum(data(iso0(is):isof(is),ico,ip)), &
                  sum(data(iso0(is):isof(is),ica,ip)) + sum(data(iso0(is):isof(is),icl,ip)) + &
                  sum(data(iso0(is):isof(is),icn,ip)) + sum(data(iso0(is):isof(is),ics,ip)) + &
                  sum(data(iso0(is):isof(is),ico,ip))
          end do
       end if

    enddo

  end subroutine med_diag_print_atm

  !===============================================================================
  subroutine med_diag_print_lnd_ice_ocn(data, ip, date, tod)

    ! ---------------------------------------------------------
    ! detail lnd/ocn/ice component budgets
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such
    integer , intent(in) :: ip
    integer , intent(in) :: date
    integer , intent(in) :: tod

    ! local variables
    integer           :: ic,nf,is ! data array indicies
    integer           :: icar,icas
    integer           :: icxs,icxr
    character(len=40) :: str      ! string
    character(*), parameter :: subName = '(med_diag_print_lnd_ice_ocn) '
    ! ------------------------------------------------------------------
    icar = 0
    icxs = 0
    icxr = 0
    icas = 0
    str = ""
    do ic = 1,4

       if (ic == 1) then
          icar = c_lnd_arecv
          icxs = c_lnd_send
          icxr = c_lnd_recv
          icas = c_lnd_asend
          str = "LND"
       elseif (ic == 2) then
          icar = c_ocn_arecv
          icxs = c_ocn_send
          icxr = c_ocn_recv
          icas = c_ocn_asend
          str = "OCN"
       elseif (ic == 3) then
          icar = c_inh_arecv
          icxs = c_inh_send
          icxr = c_inh_recv
          icas = c_inh_asend
          str = "ICE_NH"
       elseif (ic == 4) then
          icar = c_ish_arecv
          icxs = c_ish_send
          icxr = c_ish_recv
          icas = c_ish_asend
          str = "ICE_SH"
       endif

       ! heat budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh,

       write(diagunit,*) ' '
       write(diagunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',&
            trim(budget_diags%periods(ip)%name),': date = ',date,tod
       write(diagunit,FA0) budget_diags%comps(icar)%name,&
            budget_diags%comps(icxs)%name,&
            budget_diags%comps(icxr)%name,&
            budget_diags%comps(icas)%name,' *SUM*  '
       do nf = f_heat_beg, f_heat_end
          write(diagunit,FA1) budget_diags%fields(nf)%name,&
               -data(nf,icar,ip), &
                data(nf,icxs,ip), &
                data(nf,icxr,ip), &
               -data(nf,icas,ip), &
               -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
       enddo
       write(diagunit,FA1)'   *SUM*',&
            -sum(data(f_heat_beg:f_heat_end,icar,ip)), &
             sum(data(f_heat_beg:f_heat_end,icxs,ip)), &
             sum(data(f_heat_beg:f_heat_end,icxr,ip)), &
            -sum(data(f_heat_beg:f_heat_end,icas,ip)), &
            -sum(data(f_heat_beg:f_heat_end,icar,ip)) + sum(data(f_heat_beg:f_heat_end,icxs,ip)) + &
             sum(data(f_heat_beg:f_heat_end,icxr,ip)) - sum(data(f_heat_beg:f_heat_end,icas,ip))

       ! water budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh,

       write(diagunit,*) ' '
       write(diagunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',&
            trim(budget_diags%periods(ip)%name),': date = ',date,tod
       write(diagunit,FA0) &
            budget_diags%comps(icar)%name,&
            budget_diags%comps(icxs)%name,&
            budget_diags%comps(icxr)%name,&
            budget_diags%comps(icas)%name,' *SUM*  '
       do nf = f_watr_beg, f_watr_end
          write(diagunit,FA1) budget_diags%fields(nf)%name,&
               -data(nf,icar,ip),&
                data(nf,icxs,ip), &
                data(nf,icxr,ip),&
               -data(nf,icas,ip), &
               -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
       enddo
       write(diagunit,FA1)    '   *SUM*',&
            -sum(data(f_watr_beg:f_watr_end,icar,ip)), &
             sum(data(f_watr_beg:f_watr_end,icxs,ip)), &
             sum(data(f_watr_beg:f_watr_end,icxr,ip)), &
            -sum(data(f_watr_beg:f_watr_end,icas,ip)), &
            -sum(data(f_watr_beg:f_watr_end,icar,ip)) + sum(data(f_watr_beg:f_watr_end,icxs,ip)) + &
             sum(data(f_watr_beg:f_watr_end,icxr,ip)) - sum(data(f_watr_beg:f_watr_end,icas,ip))

       if ( flds_wiso ) then
          do is = 1, nisotopes

            ! heat budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh for water isotopes

             write(diagunit,*) ' '
             write(diagunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',&
                  trim(budget_diags%periods(ip)%name), &
                  ': date = ',date,tod
             write(diagunit,FA0) &
                  budget_diags%comps(icar)%name,&
                  budget_diags%comps(icxs)%name,&
                  budget_diags%comps(icxr)%name,&
                  budget_diags%comps(icas)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(diagunit,FA1) budget_diags%fields(nf)%name,&
                     -data(nf,icar,ip), &
                      data(nf,icxs,ip), &
                      data(nf,icxr,ip), &
                     -data(nf,icas,ip), &
                     -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
             enddo
             write(diagunit,FA1)    '   *SUM*',&
                  -sum(data(iso0(is):isof(is),icar,ip)),&
                   sum(data(iso0(is):isof(is),icxs,ip)), &
                   sum(data(iso0(is):isof(is),icxr,ip)), &
                  -sum(data(iso0(is):isof(is),icas,ip)), &
                  -sum(data(iso0(is):isof(is),icar,ip)) + sum(data(iso0(is):isof(is),icxs,ip)) + &
                   sum(data(iso0(is):isof(is),icxr,ip)) - sum(data(iso0(is):isof(is),icas,ip))

             ! water budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh for water isotopes

             write(diagunit,*) ' '
             write(diagunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',&
                  trim(budget_diags%periods(ip)%name),&
                  ': date = ',date,tod
             write(diagunit,FA0) &
                  budget_diags%comps(icar)%name,&
                  budget_diags%comps(icxs)%name,&
                  budget_diags%comps(icxr)%name,&
                  budget_diags%comps(icas)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(diagunit,FA1) budget_diags%fields(nf)%name,&
                     -data(nf,icar,ip), &
                      data(nf,icxs,ip), &
                      data(nf,icxr,ip), &
                     -data(nf,icas,ip), &
                     -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
             enddo
             write(diagunit,FA1)    '   *SUM*',                &
                  -sum(data(iso0(is):isof(is), icar, ip)), &
                   sum(data(iso0(is):isof(is), icxs, ip)), &
                   sum(data(iso0(is):isof(is), icxr, ip)), &
                  -sum(data(iso0(is):isof(is), icas, ip)), &
                  -sum(data(iso0(is):isof(is), icar, ip)) + sum(data(iso0(is):isof(is), icxs, ip)) + &
                   sum(data(iso0(is):isof(is), icxr, ip)) - sum(data(iso0(is):isof(is), icas, ip))
          end do
       end if
    enddo

  end subroutine med_diag_print_lnd_ice_ocn

  !===============================================================================
  subroutine med_diag_print_summary(data, ip, date, tod)

    ! ---------------------------------------------------------
    ! net summary budgets
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such
    integer , intent(in) :: ip
    integer , intent(in) :: date
    integer , intent(in) :: tod

    ! local variables
    integer  :: nf,is ! data array indicies
    real(r8) :: atm_area, lnd_area, ocn_area
    real(r8) :: ice_area_nh, ice_area_sh
    real(r8) :: sum_area
    real(r8) :: net_water_atm    , sum_net_water_atm
    real(r8) :: net_water_lnd    , sum_net_water_lnd
    real(r8) :: net_water_rof    , sum_net_water_rof
    real(r8) :: net_water_ocn    , sum_net_water_ocn
    real(r8) :: net_water_glc    , sum_net_water_glc
    real(r8) :: net_water_ice_nh , sum_net_water_ice_nh
    real(r8) :: net_water_ice_sh , sum_net_water_ice_sh
    real(r8) :: net_water_tot    , sum_net_water_tot
    real(r8) :: net_heat_atm     , sum_net_heat_atm
    real(r8) :: net_heat_lnd     , sum_net_heat_lnd
    real(r8) :: net_heat_rof     , sum_net_heat_rof
    real(r8) :: net_heat_ocn     , sum_net_heat_ocn
    real(r8) :: net_heat_glc     , sum_net_heat_glc
    real(r8) :: net_heat_ice_nh  , sum_net_heat_ice_nh
    real(r8) :: net_heat_ice_sh  , sum_net_heat_ice_sh
    real(r8) :: net_heat_tot     , sum_net_heat_tot
    real(r8) :: net_salt_atm     , sum_net_salt_atm
    real(r8) :: net_salt_lnd     , sum_net_salt_lnd
    real(r8) :: net_salt_rof     , sum_net_salt_rof
    real(r8) :: net_salt_ocn     , sum_net_salt_ocn
    real(r8) :: net_salt_glc     , sum_net_salt_glc
    real(r8) :: net_salt_ice_nh  , sum_net_salt_ice_nh
    real(r8) :: net_salt_ice_sh  , sum_net_salt_ice_sh
    real(r8) :: net_salt_tot     , sum_net_salt_tot
    character(*), parameter:: subName = '(med_diag_print_summary) '
    ! ------------------------------------------------------------------

    call t_startf('MED:'//subname)

    ! write out areas
    write(diagunit,*) ' '
    write(diagunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',&
         trim(budget_diags%periods(ip)%name),&
         ': date = ',date,tod
    write(diagunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
    atm_area    = data(f_area,c_atm_recv,ip)
    lnd_area    = data(f_area,c_lnd_recv,ip)
    ocn_area    = data(f_area,c_ocn_recv,ip)
    ice_area_nh = data(f_area,c_inh_recv,ip)
    ice_area_sh = data(f_area,c_ish_recv,ip)
    sum_area    = atm_area + lnd_area + ocn_area + ice_area_nh + ice_area_sh
    write(diagunit,FA1) budget_diags%fields(f_area)%name, atm_area, lnd_area, ocn_area, ice_area_nh, ice_area_sh, sum_area

    ! -----------------------------
    ! write out net heat budgets
    ! -----------------------------

    write(diagunit,*) ' '
    write(diagunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',&
         trim(budget_diags%periods(ip)%name), ': date = ',date,tod
    write(diagunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
    do nf = f_heat_beg, f_heat_end
       net_heat_atm    = data(nf, c_atm_recv, ip) + data(nf, c_atm_send, ip)
       net_heat_lnd    = data(nf, c_lnd_recv, ip) + data(nf, c_lnd_send, ip)
       net_heat_rof    = data(nf, c_rof_recv, ip) + data(nf, c_rof_send, ip)
       net_heat_ocn    = data(nf, c_ocn_recv, ip) + data(nf, c_ocn_send, ip)
       net_heat_ice_nh = data(nf, c_inh_recv, ip) + data(nf, c_inh_send, ip)
       net_heat_ice_sh = data(nf, c_ish_recv, ip) + data(nf, c_ish_send, ip)
       net_heat_glc    = data(nf, c_glc_recv, ip) + data(nf, c_glc_send, ip)
       net_heat_tot    = net_heat_atm + net_heat_lnd + net_heat_rof + net_heat_ocn + &
                         net_heat_ice_nh + net_heat_ice_sh + net_heat_glc

       write(diagunit,FA1r) budget_diags%fields(nf)%name,&
            net_heat_atm, net_heat_lnd, net_heat_rof, net_heat_ocn, &
            net_heat_ice_nh, net_heat_ice_sh, net_heat_glc, net_heat_tot
    end do

    ! Write out sum over all net heat budgets (sum over f_heat_beg -> f_heat_end)
    sum_net_heat_atm    = sum(data(f_heat_beg:f_heat_end, c_atm_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_atm_send, ip))
    sum_net_heat_lnd    = sum(data(f_heat_beg:f_heat_end, c_lnd_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_lnd_send, ip))
    sum_net_heat_rof    = sum(data(f_heat_beg:f_heat_end, c_rof_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_rof_send, ip))
    sum_net_heat_ocn    = sum(data(f_heat_beg:f_heat_end, c_ocn_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_ocn_send, ip))
    sum_net_heat_ice_nh = sum(data(f_heat_beg:f_heat_end, c_inh_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_inh_send, ip))
    sum_net_heat_ice_sh = sum(data(f_heat_beg:f_heat_end, c_ish_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_ish_send, ip))
    sum_net_heat_glc    = sum(data(f_heat_beg:f_heat_end, c_glc_recv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_glc_send, ip))
    sum_net_heat_tot    = sum_net_heat_atm + sum_net_heat_lnd + sum_net_heat_rof + sum_net_heat_ocn + &
                          sum_net_heat_ice_nh + sum_net_heat_ice_sh + sum_net_heat_glc

    write(diagunit,FA1r)'   *SUM*',&
         sum_net_heat_atm, sum_net_heat_lnd, sum_net_heat_rof, sum_net_heat_ocn, &
         sum_net_heat_ice_nh, sum_net_heat_ice_sh, sum_net_heat_glc, sum_net_heat_tot

    ! -----------------------------
    ! write out net water budgets
    ! -----------------------------

    write(diagunit,*) ' '
    write(diagunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',&
         trim(budget_diags%periods(ip)%name), ': date = ',date,tod
    write(diagunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
    do nf = f_watr_beg, f_watr_end
       net_water_atm    = data(nf, c_atm_recv, ip) + data(nf, c_atm_send, ip)
       net_water_lnd    = data(nf, c_lnd_recv, ip) + data(nf, c_lnd_send, ip)
       net_water_rof    = data(nf, c_rof_recv, ip) + data(nf, c_rof_send, ip)
       net_water_ocn    = data(nf, c_ocn_recv, ip) + data(nf, c_ocn_send, ip)
       net_water_ice_nh = data(nf, c_inh_recv, ip) + data(nf, c_inh_send, ip)
       net_water_ice_sh = data(nf, c_ish_recv, ip) + data(nf, c_ish_send, ip)
       net_water_glc    = data(nf, c_glc_recv, ip) + data(nf, c_glc_send, ip)
       net_water_tot    = net_water_atm + net_water_lnd + net_water_rof + net_water_ocn + &
                          net_water_ice_nh + net_water_ice_sh + net_water_glc

       write(diagunit,FA1r) budget_diags%fields(nf)%name,&
            net_water_atm, net_water_lnd, net_water_rof, net_water_ocn, &
            net_water_ice_nh, net_water_ice_sh, net_water_glc, net_water_tot
    enddo

    ! Write out sum over all net water budgets (sum over f_watr_beg -> f_watr_end)
    sum_net_water_atm    = sum(data(f_watr_beg:f_watr_end, c_atm_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_atm_send, ip))
    sum_net_water_lnd    = sum(data(f_watr_beg:f_watr_end, c_lnd_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_lnd_send, ip))
    sum_net_water_rof    = sum(data(f_watr_beg:f_watr_end, c_rof_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_rof_send, ip))
    sum_net_water_ocn    = sum(data(f_watr_beg:f_watr_end, c_ocn_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_ocn_send, ip))
    sum_net_water_ice_nh = sum(data(f_watr_beg:f_watr_end, c_inh_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_inh_send, ip))
    sum_net_water_ice_sh = sum(data(f_watr_beg:f_watr_end, c_ish_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_ish_send, ip))
    sum_net_water_glc    = sum(data(f_watr_beg:f_watr_end, c_glc_recv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_glc_send, ip))
    sum_net_water_tot    = sum_net_water_atm + sum_net_water_lnd + sum_net_water_rof + sum_net_water_ocn + &
                           sum_net_water_ice_nh + sum_net_water_ice_sh + sum_net_water_glc

    write(diagunit,FA1r)'   *SUM*',&
         sum_net_water_atm, sum_net_water_lnd, sum_net_water_rof, sum_net_water_ocn, &
         sum_net_water_ice_nh, sum_net_water_ice_sh, sum_net_water_glc, sum_net_water_tot

    ! write out net water water-isoptope budgets

    if ( flds_wiso ) then

       do is = 1, nisotopes
          write(diagunit,*) ' '
          write(diagunit,FAH) subname,'NET '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
               trim(budget_diags%periods(ip)%name),': date = ',date,tod
          write(diagunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
          do nf = iso0(is), isof(is)
             net_water_atm    = data(nf, c_atm_recv, ip) + data(nf, c_atm_send, ip)
             net_water_lnd    = data(nf, c_lnd_recv, ip) + data(nf, c_lnd_send, ip)
             net_water_rof    = data(nf, c_rof_recv, ip) + data(nf, c_rof_send, ip)
             net_water_ocn    = data(nf, c_ocn_recv, ip) + data(nf, c_ocn_send, ip)
             net_water_ice_nh = data(nf, c_inh_recv, ip) + data(nf, c_inh_send, ip)
             net_water_ice_sh = data(nf, c_ish_recv, ip) + data(nf, c_ish_send, ip)
             net_water_glc    = data(nf, c_glc_recv, ip) + data(nf, c_glc_send, ip)
             net_water_tot    = net_water_atm + net_water_lnd + net_water_rof + net_water_ocn + &
                                net_water_ice_nh + net_water_ice_sh + net_water_glc

             write(diagunit,FA1r) budget_diags%fields(nf)%name,&
                  net_water_atm, net_water_lnd, net_water_rof, net_water_ocn, &
                  net_water_ice_nh, net_water_ice_sh, net_water_glc, net_water_tot
          enddo

          sum_net_water_atm    = sum(data(iso0(is):isof(is), c_atm_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_atm_send, ip))
          sum_net_water_lnd    = sum(data(iso0(is):isof(is), c_lnd_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_lnd_send, ip))
          sum_net_water_rof    = sum(data(iso0(is):isof(is), c_rof_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_rof_send, ip))
          sum_net_water_ocn    = sum(data(iso0(is):isof(is), c_ocn_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_ocn_send, ip))
          sum_net_water_ice_nh = sum(data(iso0(is):isof(is), c_inh_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_inh_send, ip))
          sum_net_water_ice_sh = sum(data(iso0(is):isof(is), c_ish_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_ish_send, ip))
          sum_net_water_glc    = sum(data(iso0(is):isof(is), c_glc_recv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_glc_send, ip))
          sum_net_water_tot    = sum_net_water_atm + sum_net_water_lnd + sum_net_water_rof + &
                                 sum_net_water_ocn + sum_net_water_ice_nh + sum_net_water_ice_sh + &
                                 sum_net_water_glc

          write(diagunit,FA1r)'   *SUM*',&
               sum_net_water_atm, sum_net_water_lnd, sum_net_water_rof, sum_net_water_ocn, &
               sum_net_water_ice_nh, sum_net_water_ice_sh, sum_net_water_glc, sum_net_water_tot
       end do
    end if

    ! -----------------------------
    ! write out net salt budgets
    ! -----------------------------

    if (trim(budget_table_version) == 'v1') then
       write(diagunit,*) ' '
       write(diagunit,FAH) subname,'NET SALT BUDGET (kg/m2s): period = ',&
            trim(budget_diags%periods(ip)%name), ': date = ',date,tod
       write(diagunit,FA0s) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
       do nf = f_salt_beg, f_salt_end
          net_salt_atm    = data(nf, c_atm_recv, ip) + data(nf, c_atm_send, ip)
          net_salt_lnd    = data(nf, c_lnd_recv, ip) + data(nf, c_lnd_send, ip)
          net_salt_rof    = data(nf, c_rof_recv, ip) + data(nf, c_rof_send, ip)
          net_salt_ocn    = data(nf, c_ocn_recv, ip) + data(nf, c_ocn_send, ip)
          net_salt_ice_nh = data(nf, c_inh_recv, ip) + data(nf, c_inh_send, ip)
          net_salt_ice_sh = data(nf, c_ish_recv, ip) + data(nf, c_ish_send, ip)
          net_salt_glc    = data(nf, c_glc_recv, ip) + data(nf, c_glc_send, ip)
          net_salt_tot    = net_salt_atm + net_salt_lnd + net_salt_rof + net_salt_ocn + &
               net_salt_ice_nh + net_salt_ice_sh + net_salt_glc

          write(diagunit,FA1s) budget_diags%fields(nf)%name,&
               net_salt_atm, net_salt_lnd, net_salt_rof, net_salt_ocn, &
               net_salt_ice_nh, net_salt_ice_sh, net_salt_glc, net_salt_tot
       enddo

       ! Write out sum over all net heat budgets (sum over f_salt_beg -> f_salt_end)
       sum_net_salt_atm    = sum(data(f_salt_beg:f_salt_end, c_atm_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_atm_send, ip))
       sum_net_salt_lnd    = sum(data(f_salt_beg:f_salt_end, c_lnd_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_lnd_send, ip))
       sum_net_salt_rof    = sum(data(f_salt_beg:f_salt_end, c_rof_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_rof_send, ip))
       sum_net_salt_ocn    = sum(data(f_salt_beg:f_salt_end, c_ocn_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_ocn_send, ip))
       sum_net_salt_ice_nh = sum(data(f_salt_beg:f_salt_end, c_inh_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_inh_send, ip))
       sum_net_salt_ice_sh = sum(data(f_salt_beg:f_salt_end, c_ish_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_ish_send, ip))
       sum_net_salt_glc    = sum(data(f_salt_beg:f_salt_end, c_glc_recv, ip)) + &
            sum(data(f_salt_beg:f_salt_end, c_glc_send, ip))
       sum_net_salt_tot    = sum_net_salt_atm + sum_net_salt_lnd + sum_net_salt_rof + sum_net_salt_ocn + &
            sum_net_salt_ice_nh + sum_net_salt_ice_sh + sum_net_salt_glc

       write(diagunit,FA1s)'   *SUM*',&
            sum_net_salt_atm, sum_net_salt_lnd, sum_net_salt_rof, sum_net_salt_ocn, &
            sum_net_salt_ice_nh, sum_net_salt_ice_sh, sum_net_salt_glc, sum_net_salt_tot
    end if

    call t_stopf('MED:'//subname)
  end subroutine med_diag_print_summary

  !===============================================================================
  subroutine add_to_budget_diag(entries, index, name)

    ! input/output variablesn
    type(budget_diag_type) , pointer       :: entries(:)
    integer                , intent(out)   :: index
    character(len=*)       , intent(in)    :: name

    ! local variables
    integer :: n
    integer :: oldsize
    logical :: found
    type(budget_diag_type), pointer :: new_entries(:)
    character(len=*), parameter :: subname='(add_to_budget_diag)'
    !----------------------------------------------------------------------

    if (associated(entries)) then
       oldsize = size(entries)
       found = .false.
       do n= 1,oldsize
          if (trim(name) == trim(entries(n)%name)) then
             found = .true.
             exit
          end if
       end do
    else
       oldsize = 0
       found = .false.
    end if
    index = oldsize + 1

    ! create new entry if fldname is not in original list

    if (.not. found) then
       if(maintask) write(logunit,*) ' Add ',trim(name),' to budgets with index ',index
       ! 1) allocate newfld to be size (one element larger than input flds)
       allocate(new_entries(index))

       ! 2) copy entries into first N-1 elements of new_entries
       do n = 1,oldsize
          new_entries(n)%name = entries(n)%name
       end do

       ! 3) deallocate / nullify entries
       if (oldsize >  0) then
          deallocate(entries)
          nullify(entries)
       end if
       entries => new_entries

       ! 4) point entries => new_entries
       entries => new_entries

       ! 5) now update entries information for new entry
       entries(index)%name  = trim(name)
    end if

  end subroutine add_to_budget_diag

end module med_diag_mod
