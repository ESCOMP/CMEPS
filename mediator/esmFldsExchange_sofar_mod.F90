module esmFldsExchange_sofar_mod

  use ESMF
  use NUOPC
  use med_utils_mod         , only : chkerr => med_utils_chkerr
  use med_kind_mod          , only : CX=>SHR_KIND_CX
  use med_kind_mod          , only : CS=>SHR_KIND_CS
  use med_kind_mod          , only : CL=>SHR_KIND_CL
  use med_kind_mod          , only : R8=>SHR_KIND_R8
  use med_internalstate_mod , only : compmed
  use med_internalstate_mod , only : compatm
  use med_internalstate_mod , only : compocn
  use med_internalstate_mod , only : compwav
  use med_internalstate_mod , only : compice    ! Sofar added: for eventual coupling of cice6
  use med_internalstate_mod , only : ncomps
  use med_internalstate_mod , only : coupling_mode
  use esmFlds               , only : addfld_ocnalb => med_fldList_addfld_ocnalb

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_sofar

  character(*), parameter :: u_FILE_u = &
       __FILE__

  type gcomp_attr
    character(len=CX) :: atm2ocn_fmap = 'unset'
    character(len=CX) :: atm2ocn_smap = 'unset'
    character(len=CX) :: atm2ocn_vmap = 'unset'
    character(len=CX) :: atm2wav_smap = 'unset'
    character(len=CX) :: ocn2atm_fmap = 'unset'
    character(len=CX) :: ocn2atm_smap = 'unset'
    character(len=CX) :: ocn2wav_smap = 'unset'
    character(len=CX) :: wav2ocn_smap = 'unset'
    character(len=CX) :: wav2atm_smap = 'unset'
    character(len=CS) :: mapnorm      = 'one'
    logical           :: atm_present  = .false.
    logical           :: ocn_present  = .false.
    logical           :: wav_present  = .false.
  end type

!===============================================================================
contains
!===============================================================================

  subroutine esmFldsExchange_sofar(gcomp, phase, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    character(len=*) , parameter   :: subname='(esmFldsExchange_sofar)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    if (phase == 'advertise') then
      call esmFldsExchange_sofar_advt(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'fieldcheck') then
      call esmFldsExchange_sofar_fchk(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'initialize') then
      call esmFldsExchange_sofar_init(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogSetError(ESMF_FAILURE, &
         msg=trim(subname)//": Phase is set to "//trim(phase), &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_sofar

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_sofar_advt(gcomp, phase, rc)

    use esmFlds, only : addfld_to => med_fldList_addfld_to
    use esmFlds, only : addfld_from => med_fldList_addfld_from

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    integer             :: n
    logical             :: isPresent
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: sofar_attr
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_sofar_advt)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !=====================================================================
    ! scalar information
    !=====================================================================

    call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld_from(n, trim(cvalue))
          call addfld_to(n, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_sofar_attr(gcomp, sofar_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    if (sofar_attr%ocn_present) then                 ! Sofar system: added
        call addfld_from(compocn, 'So_omask')
    endif                                            ! Sofar system: added

    !----------------------------------------------------------
    ! to med: frac from components
    !----------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then    ! Sofar system: added
!       call addfld_to(compatm, 'So_ofrac')                          ! Sofar system: added
        call addfld_from(compatm , 'Sa_oceanfrac')                   
    endif                                                            ! Sofar system: added 

    !----------------------------------------------------------
    ! from med: ocean albedos (not sent to the ATM in UFS).
    !----------------------------------------------------------
    if (trim(coupling_mode(1:5)) == 'sofar') then
       if (phase == 'advertise') then
          call addfld_ocnalb('So_avsdr')
          call addfld_ocnalb('So_avsdf')
          call addfld_ocnalb('So_anidr')
          call addfld_ocnalb('So_anidf')
       end if
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to atm: surface temperatures from ocn
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
       if (trim(coupling_mode(1:5)) == 'sofar') then
          allocate(S_flds(1))
          S_flds = (/'So_t'/)  ! sea_surface_temperature
       else
           allocate(S_flds(0))
       end if
       do n = 1,size(S_flds)
           fldname = trim(S_flds(n))
           call addfld_from(compocn, trim(fldname))
           call addfld_to(compatm, trim(fldname))
       end do
       if (allocated(S_flds)) deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Charnock parameter and surface roughness length
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%wav_present) then
       allocate(S_flds(1))
       !allocate(S_flds(2))                                          ! Sofar system !ISSUE: add option to change the export vars at runtime 
       S_flds = (/'Sw_charno'/)    ! Charnock parameter
       !S_flds = (/'Sw_z0rlen', &  ! wave_z0_roughness_length        ! Sofar system
       !           'Sw_charno', &  ! Charnock parameter              ! Sofar system
       !          /)                                                 ! Sofar system
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          call addfld_from(compwav, trim(fldname))
          call addfld_to(compatm, trim(fldname))
       end do
       deallocate(S_flds)
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
        if (trim(coupling_mode(1:5)) == 'sofar') then
            allocate(S_flds(1))
            S_flds = (/'Sa_pslv'/) ! inst_pres_height_surface
        else
            allocate(S_flds(0))
        endif
        do n = 1,size(S_flds)
            fldname = trim(S_flds(n))
            call addfld_from(compatm, trim(fldname))
            call addfld_to(compocn, trim(fldname))
        enddo
        if (allocated(S_flds)) deallocate(S_flds)
    endif

    ! ---------------------------------------------------------------------
    ! to ocn: flux fields
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
        if (trim(coupling_mode) == 'sofar.nofluxes') then
            print *, "esmFldsExchange_sofar_mod.F90:: coupling_mode==sofar.test, Skip fluxes..."
        elseif (trim(coupling_mode(1:5)) == 'sofar') then
            allocate(F_flds(10,2))
            F_flds(1 ,:) = (/'Faxa_taux     ','Foxx_taux     '/) ! inst_zonal_moment_flx_atm
            F_flds(2 ,:) = (/'Faxa_tauy     ','Foxx_tauy     '/) ! inst_merid_moment_flx_atm
            F_flds(3 ,:) = (/'Faxa_rain     ','Faxa_rain     '/) ! inst_prec_rate
            F_flds(4 ,:) = (/'Faxa_lwnet    ','Foxx_lwnet    '/) ! inst_net_lw_flx
            F_flds(5 ,:) = (/'Faxa_sen      ','Foxx_sen      '/) ! inst_sensi_heat_flx
            F_flds(6 ,:) = (/'Faxa_evap     ','Foxx_evap     '/) ! inst_evap_rate
            F_flds(7 ,:) = (/'Faxa_swndr    ','Foxx_swnet_idr'/) ! inst_down_sw_ir_dir_flx
            F_flds(8 ,:) = (/'Faxa_swndf    ','Foxx_swnet_idf'/) ! inst_down_sw_ir_dif_flx
            F_flds(9 ,:) = (/'Faxa_swvdr    ','Foxx_swnet_vdr'/) ! inst_down_sw_vis_dir_flx
            F_flds(10,:) = (/'Faxa_swvdf    ','Foxx_swnet_vdf'/) ! inst_down_sw_vis_dif_flx
        else
            allocate(F_flds(0,1))
        endif
        do n = 1,size(F_flds,1)
            fldname1 = trim(F_flds(n,1))
            fldname2 = trim(F_flds(n,2))
            call addfld_from(compatm, trim(fldname1))
            call addfld_to(compocn, trim(fldname2))
        end do
        if (allocated(F_flds)) deallocate(F_flds)
    endif

    ! ---------------------------------------------------------------------
    ! to ocn: wave parameters
    ! ---------------------------------------------------------------------
    if (sofar_attr%wav_present .and. sofar_attr%ocn_present) then
        ! See here for fields that the ocean model can actually accept:
        ! https://github.com/NOAA-GFDL/MOM6/blob/2f2b7905c08e95a729d3dd3f8b02e0a0bed10602/config_src/drivers/nuopc_cap/mom_cap.F90#L787
        if (trim(coupling_mode) == 'sofar.wav2ocn') then
            S_flds = (/'Sw_uscurr', &      ! Stokes Drift 3D
                       'Sw_vscurr', &      ! 
                       'Sw_x1pstk', &      ! Partitioned Stokes Drift 3 2D fields
                       'Sw_y1pstk', &      ! 
                       'Sw_x2pstk', &      ! 
                       'Sw_y2pstk', &      ! 
                       'Sw_x3pstk', &      ! 
                       'Sw_y3pstk', &      ! 
                       'Sw_wbcuru', &      ! Bottom Currents
                       'Sw_wbcurv', &      ! 
                       'Sw_wbcurp', &      ! 
                       'Sw_wavsuu', &      ! Radiation stresses 2D
                       'Sw_wavsuv', &      ! 
                       'Sw_wavsvv' &       ! 
                      /)                   
        else
            allocate(S_flds(0))
        endif
        do n = 1,size(S_flds)
            fldname = trim(S_flds(n))
            call addfld_from(compwav, trim(fldname))
            call addfld_to(compocn, trim(fldname))
        enddo
        deallocate(S_flds)
        if (allocated(F_flds)) deallocate(F_flds)
    endif

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to wav: 10-m wind components, air surface density, and air-sea temp difference
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%wav_present) then
        allocate(S_flds(4))          ! Sofar system !ISSUE: add option to change the export vars at runtime
        S_flds = (/'Sa_u10n', &      ! zonal diagnosed 10m neutral wind component
                   'Sa_v10n', &      ! meridional diagnosed 10m neutral wind component
                   'Sa_rhoa', &      ! atmospheric surface density
                   'Sa_astdiff' &    ! air minus sea surface temperature difference
                  /)                 ! Sofar system with diagnosed 10m winds
        do n = 1,size(S_flds)
            fldname = trim(S_flds(n))
            call addfld_from(compatm, trim(fldname))
            call addfld_to(compwav, trim(fldname))
        enddo
        deallocate(S_flds)
    endif

    ! ---------------------------------------------------------------------
    ! to wav: ocean surface components
    ! ---------------------------------------------------------------------
    if (sofar_attr%ocn_present .and. sofar_attr%wav_present) then
        if (trim(coupling_mode) == 'sofar.ocn2wav')
            allocate(S_flds(3))          ! Sofar system !ISSUE: add option to change the export vars at runtime
            S_flds = (/'So_u', &         ! zonal ocean surface current
                       'So_v', &         ! meridional ocean surface current
                       'So_ssh'  &       ! Sea surface height
                   !   'So_rhoo', &      ! ocean surface density
                   !   'So_t' &          ! ocean surface temperature (eventually pass gustiness from atm to wav)
                      /)                 
        else
            allocate(S_flds(0))
        endif
        do n = 1,size(S_flds)
            fldname = trim(S_flds(n))
            call addfld_from(compocn, trim(fldname))
            call addfld_to(compwav, trim(fldname))
        enddo
        deallocate(S_flds)
    endif

    ! ---------------------------------------------------------------------
    ! to wav: sea ice (not yet supported)
    ! ---------------------------------------------------------------------
    if (sofar_attr%ocn_present .and. sofar_attr%wav_present) then
        if (trim(coupling_mode) == 'sofar.ice2wav')
            allocate(S_flds(1))          ! Sofar system !ISSUE: add option to change the export vars at runtime
            S_flds = (/'Si_seaice'       ! Sea ice fraction / concentration
                      /)                 
        else
            allocate(S_flds(0))
        endif
        do n = 1,size(S_flds)
            fldname = trim(S_flds(n))
            call addfld_from(compice, trim(fldname))
            call addfld_to(compwav, trim(fldname))
        enddo
        deallocate(S_flds)
    endif


    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_sofar_advt

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_sofar_fchk(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    character(len=*) , parameter   :: subname='(esmFldsExchange_sofar_fchk)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fldchk(is_local%wrap%FBImp(compocn,compocn),'So_omask',rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": Field connected "//"So_omask", &
          ESMF_LOGMSG_INFO)
    else
       call ESMF_LogSetError(ESMF_FAILURE, &
          msg=trim(subname)//": Field is not connected "//"So_omask", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_sofar_fchk

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_sofar_init(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd
    use med_internalstate_mod , only : mapfillv_bilnr
    use med_internalstate_mod , only : mapnstod_consf
    use esmFlds               , only : addmap_from => med_fldList_addmap_from
    use esmFlds               , only : addmrg_to   => med_fldList_addmrg_to

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: n
    character(len=CS)   :: fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: sofar_attr
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_sofar_init)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------
    ! Merging arguments:
    ! mrg_fromN = source component index that for the field to be merged
    ! mrg_fldN  = souce field name to be merged
    ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum',
    !                         'sum_with_weights', 'merge')
    ! NOTE:
    ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn
    ! fluxes or for ocn albedos
    !
    ! NOTE:
    ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in
    ! the atm/ocn flux calculations. Input fields required from either the atm
    ! or the ocn for these computation will use the logical 'use_med_aoflux'
    ! below. This is used to determine mappings between the atm and ocn needed
    ! for these computations.
    !--------------------------------------

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_sofar_attr(gcomp, sofar_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to atm: sea surface temperature
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
       if (trim(coupling_mode) == 'sofar') then
          allocate(S_flds(1))
          S_flds = (/'So_t'/)  ! sea_surface_temperature
                 !   'So_u', & ! surface zonal current
                 !   'So_v'/)  ! surface meridional current
       elseif (trim(coupling_mode) == 'sofar.test') then
          allocate(S_flds(1))
          S_flds = (/'So_t'/) ! sea_surface_temperature
       elseif (trim(coupling_mode) == 'sofar.hycom') then
          allocate(S_flds(1))
          S_flds = (/'So_t'/) ! sea_surface_temperature
       else
           allocate(S_flds(0))
       endif
       do n = 1,size(S_flds)
           fldname = trim(S_flds(n))
           if (fldchk(is_local%wrap%FBExp(compatm),trim(fldname),rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
               ) then
               call addmap_from(compocn, trim(fldname), compatm, &
                     mapfillv_bilnr, sofar_attr%mapnorm, sofar_attr%ocn2atm_smap)
               call addmrg_to(compatm, trim(fldname), &
                     mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
           endif
       enddo
       deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Charnock parameter and surface roughness length
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%wav_present) then
       allocate(S_flds(1))
       !allocate(S_flds(2))                                          ! Sofar system !ISSUE: add option to change the export vars at runtime
!      S_flds = (/'Sw_z0rlen'/)    ! wave_z0_roughness_length
       S_flds = (/'Sw_charno'/)    ! Charnock parameter
       !S_flds = (/'Sw_z0rlen', &  ! wave_z0_roughness_length        ! Sofar system
       !           'Sw_charno', &  ! Charnock parameter              ! Sofar system
       !          /)
       do n = 1,size(S_flds)
          fldname = trim(S_flds(n))
          if (fldchk(is_local%wrap%FBExp(compatm),trim(fldname),rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compwav,compwav),trim(fldname),rc=rc) &
             ) then
             call addmap_from(compwav, trim(fldname), compatm, &
                  mapfillv_bilnr, sofar_attr%mapnorm, sofar_attr%wav2atm_smap)
             call addmrg_to(compatm, trim(fldname), &
                  mrg_from=compwav, mrg_fld=trim(fldname), mrg_type='copy')
          endif
       enddo
       deallocate(S_flds)
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
       if (trim(coupling_mode) == 'sofar') then
          allocate(S_flds(1))
          S_flds = (/'Sa_pslv'/) ! inst_pres_height_surface
       elseif (trim(coupling_mode) == 'sofar.test') then
          allocate(S_flds(1))
          S_flds = (/'Sa_t2m ' /) ! inst_temp_height2m
       elseif (trim(coupling_mode) == 'sofar.hycom') then
          allocate(S_flds(6))
          S_flds = (/'Sa_u10n', & ! inst_zonal_wind_height10m
                     'Sa_v10n', & ! inst_merid_wind_height10m
                     'Sa_t2m ', & ! inst_temp_height2m
                     'Sa_q2m ', & ! inst_spec_humid_height2m
                     'Sa_pslv', & ! inst_pres_height_surface
                     'Sa_tskn' /) ! inst_temp_height_surface
       else
           allocate(S_flds(0))
       endif
       do n = 1,size(S_flds)
           fldname = trim(S_flds(n))
           if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
               ) then
               call addmap_from(compatm, trim(fldname), compocn, &
                     mapfillv_bilnr, sofar_attr%mapnorm, sofar_attr%atm2ocn_smap)
               call addmrg_to(compocn, trim(fldname), &
                     mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
           endif
       enddo
       if (allocated(S_flds)) deallocate(S_flds)
    endif

    ! ---------------------------------------------------------------------
    ! to ocn: flux fields
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%ocn_present) then
       if (trim(coupling_mode) == 'sofar') then
          allocate(F_flds(10,2))
          F_flds(1 ,:) = (/'Faxa_taux     ','Foxx_taux     '/) ! inst_zonal_moment_flx_atm
          F_flds(2 ,:) = (/'Faxa_tauy     ','Foxx_tauy     '/) ! inst_merid_moment_flx_atm
          F_flds(3 ,:) = (/'Faxa_rain     ','Faxa_rain     '/) ! inst_prec_rate
          F_flds(4 ,:) = (/'Faxa_lwnet    ','Foxx_lwnet    '/) ! inst_net_lw_flx
          F_flds(5 ,:) = (/'Faxa_sen      ','Foxx_sen      '/) ! inst_sensi_heat_flx
          F_flds(6 ,:) = (/'Faxa_evap     ','Foxx_evap     '/) ! inst_evap_rate
          F_flds(7 ,:) = (/'Faxa_swndr    ','Foxx_swnet_idr'/) ! inst_down_sw_ir_dir_flx
          F_flds(8 ,:) = (/'Faxa_swndf    ','Foxx_swnet_idf'/) ! inst_down_sw_ir_dif_flx
          F_flds(9 ,:) = (/'Faxa_swvdr    ','Foxx_swnet_vdr'/) ! inst_down_sw_vis_dir_flx
          F_flds(10,:) = (/'Faxa_swvdf    ','Foxx_swnet_vdf'/) ! inst_down_sw_vis_dif_flx
       elseif (trim(coupling_mode) == 'sofar.test') then
          print *, "esmFldsExchange_sofar_mod.F90:: coupling_mode==sofar.test, Skip fluxes..."
       elseif (trim(coupling_mode) == 'sofar.hycom') then
          allocate(F_flds(7,2))
          F_flds(1,:) = (/'Faxa_taux ','Faxa_taux '/) ! inst_zonal_moment_flx_atm
          F_flds(2,:) = (/'Faxa_tauy ','Faxa_tauy '/) ! inst_merid_moment_flx_atm
          F_flds(3,:) = (/'Faxa_rain ','Faxa_rain '/) ! inst_prec_rate
          F_flds(4,:) = (/'Faxa_swnet','Faxa_swnet'/) ! inst_net_sw_flx
          F_flds(5,:) = (/'Faxa_lwnet','Faxa_lwnet'/) ! inst_net_lw_flx
          F_flds(6,:) = (/'Faxa_sen  ','Faxa_sen  '/) ! inst_sensi_heat_flx
          F_flds(7,:) = (/'Faxa_lat  ','Faxa_lat  '/) ! inst_laten_heat_flx
       else
           allocate(F_flds(0,1))
       endif
       do n = 1,size(F_flds,1)
           fldname1 = trim(F_flds(n,1))
           fldname2 = trim(F_flds(n,2))
           if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname2),rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname1),rc=rc) &
               ) then
               call addmap_from(compatm, trim(fldname1), compocn, &
                     mapfillv_bilnr, sofar_attr%mapnorm, sofar_attr%atm2ocn_smap)
               call addmrg_to(compocn, trim(fldname2), &
                     mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')
           endif
       enddo
       deallocate(F_flds)
    end if

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to wav: 10-m wind components
    ! ---------------------------------------------------------------------
    if (sofar_attr%atm_present .and. sofar_attr%wav_present) then
!     allocate(S_flds(2))
!     S_flds = (/'Sa_u10m', 'Sa_v10m'/)
      allocate(S_flds(4))          ! Sofar system !ISSUE: add option to change the export vars at runtime
      S_flds = (/'Sa_u10n', &      ! zonal diagnosed 10m neutral wind component
                 'Sa_v10n', &      ! meridional diagnosed 10m neutral wind component
                 'Sa_rhoa', &      ! atmospheric surface density
                 'Sa_astdiff' &    ! air minus sea surface temperature difference
                /)                 ! Sofar system with diagnosed 10m winds
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        if (fldchk(is_local%wrap%FBexp(compwav),trim(fldname),rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname),rc=rc) &
           ) then
           call addmap_from(compatm, trim(fldname), compwav, &
                mapfillv_bilnr, sofar_attr%mapnorm, sofar_attr%atm2wav_smap)
           call addmrg_to(compwav, trim(fldname), &
                mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
        end if
      end do
      deallocate(S_flds)
    end if

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_sofar_init

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_sofar_attr(gcomp, sofar_attr, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    type(gcomp_attr) , intent(inout) :: sofar_attr
    integer          , intent(inout) :: rc

    ! local variables:
    character(32)       :: cname
    integer             :: verbosity, diagnostic
    character(len=CL)   :: cvalue
    logical             :: isPresent, isSet
    character(len=*) , parameter   :: subname='(esmFldsExchange_sofar_attr)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call NUOPC_CompGet(gcomp, name=cname, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Component active or not?
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'satm') sofar_attr%atm_present = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='OCN_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'socn') sofar_attr%ocn_present = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='WAV_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'swav') sofar_attr%wav_present = .true.
    end if

    !----------------------------------------------------------
    ! Normalization type
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='normalization', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='normalization', &
          value=sofar_attr%mapnorm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
          value=sofar_attr%ocn2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
          value=sofar_attr%ocn2atm_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ocn
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
          value=sofar_attr%atm2ocn_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       value=sofar_attr%atm2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
          value=sofar_attr%atm2ocn_vmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to wav
    call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', &
          value=sofar_attr%atm2wav_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', &
          value=sofar_attr%ocn2wav_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! from wav
    call NUOPC_CompAttributeGet(gcomp, name='wav2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='wav2atm_smapname', &
          value=sofar_attr%wav2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smapname', &
          value=sofar_attr%wav2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Log Attribute Settings
    if (btest(verbosity,16)) then
       write(cvalue,"(I0)") verbosity
       call ESMF_LogWrite(trim(subname)//': Verbosity        = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       write(cvalue,"(I0)") diagnostic
       call ESMF_LogWrite(trim(subname)//': Diagnostic       = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': normalization    = '// &
          trim(sofar_attr%mapnorm), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_smapname = '// &
          trim(sofar_attr%ocn2atm_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_fmapname = '// &
          trim(sofar_attr%ocn2atm_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_fmapname = '// &
          trim(sofar_attr%atm2ocn_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_smapname = '// &
          trim(sofar_attr%atm2ocn_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_vmapname = '// &
          trim(sofar_attr%atm2ocn_vmap), ESMF_LOGMSG_INFO)
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_sofar_attr

end module esmFldsExchange_sofar_mod
