module esmFldsExchange_ufs_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_ufs

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange_ufs(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use med_internalstate_mod , only : compmed, compatm, compocn, compice, complnd, compwav, ncomps
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use med_internalstate_mod , only : mapconsf_aofrac, mapbilnr_nstod
    use med_internalstate_mod , only : coupling_mode, mapnames
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld_to => med_fldList_addfld_to
    use esmFlds               , only : addmrg_to => med_fldList_addmrg_to
    use esmFlds               , only : addfld_from => med_fldList_addfld_from
    use esmFlds               , only : addmap_from => med_fldList_addmap_from
    use esmFlds               , only : addfld_aoflux => med_fldList_addfld_aoflux
    use esmFlds               , only : addmap_aoflux => med_fldList_addmap_aoflux
    use esmFlds               , only : addfld_ocnalb => med_fldList_addfld_ocnalb
    use esmFlds               , only : addmap_ocnalb => med_fldList_addmap_ocnalb

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: n, maptype
    logical             :: med_aoflux_to_ocn
    character(len=CX)   :: msgString
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS), allocatable :: flds(:), oflds(:), aflds(:), iflds(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_ufs)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set maptype according to coupling_mode
    if (trim(coupling_mode) == 'ufs.nfrac' .or. trim(coupling_mode) == 'ufs.nfrac.aoflux') then
      maptype = mapnstod_consf
    else
      maptype = mapconsf
    end if
    write(msgString,'(A,i6,A)') trim(subname)//': maptype is ',maptype,', '//mapnames(maptype)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    if (trim(coupling_mode) == 'ufs.nfrac.aoflux' .or. trim(coupling_mode) == 'ufs.frac.aoflux') then
       med_aoflux_to_ocn = .true.
    else
       med_aoflux_to_ocn = .false.
    end if

    !=====================================================================
    ! scalar information
    !=====================================================================

    if (phase == 'advertise') then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld_to(n   , trim(cvalue))
          call addfld_from(n , trim(cvalue))
       end do
    end if

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    ! masks from components
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice)) call addfld_from(compice, 'Si_imask')
       if (is_local%wrap%comp_present(compocn)) call addfld_from(compocn, 'So_omask')
       if (is_local%wrap%comp_present(complnd)) call addfld_from(complnd, 'Sl_lfrin')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_omask', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_omask', rc=rc)) then
          call addmap_from(compocn, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
       end if
    end if

    ! fields required for atm/ocn flux calculation
    if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
       ! from atm: states for fluxes
       allocate(flds(13))
       flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', 'Sa_pslv', &
                'Sa_shum', 'Sa_ptem', 'Sa_dens', 'Sa_u10m', 'Sa_v10m', 'Sa_t2m ', &
                'Sa_q2m '/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             call addfld_from(compatm , fldname)
          else
             if ( fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
                call addmap_from(compatm, fldname, compocn, maptype, 'one', 'unset')
             end if
          end if
       end do
       deallocate(flds)

       ! from med: fields returned by the atm/ocn flux computation, otherwise unadvertised
       allocate(flds(12))
       flds = (/'So_tref       ', 'So_qref       ', 'So_ustar      ', 'So_re         ', 'So_ssq        ', &
                'So_u10        ', 'So_duu10n     ', 'Faox_lat      ', 'So_ugustOut   ', 'So_u10withGust', &
                'So_u10res     ', 'Faxa_rainc    '/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             call addfld_aoflux(fldname)
          end if
       end do
       deallocate(flds)
    end if

    ! from med: ocean albedos (not sent to the ATM in UFS).
    if (phase == 'advertise') then
       call addfld_ocnalb('So_avsdr')
       call addfld_ocnalb('So_avsdf')
       call addfld_ocnalb('So_anidr')
       call addfld_ocnalb('So_anidf')
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE (compatm)
    !=====================================================================

    ! to atm: fractions (computed in med_phases_prep_atm)
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compice , 'Si_ifrac')
          call addfld_to(compatm   , 'Si_ifrac')
       end if
       ! ofrac used by atm
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compatm , 'Sa_ofrac')
       end if
       ! lfrac used by atm
       if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_to(compatm  , 'Sl_lfrac')
       end if
    end if

    ! to atm: unmerged from ice
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux,
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux
    ! - evaporation water flux from water
    ! - mean ice volume per unit area
    ! - mean snow volume per unit area
    ! - surface temperatures
    allocate(flds(9))
    flds = (/'Faii_taux', 'Faii_tauy', 'Faii_lat ', 'Faii_sen ', 'Faii_lwup', 'Faii_evap', &
             'Si_vice  ', 'Si_vsno  ', 'Si_t     '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
             call addfld_from(compice , fldname)
             call addfld_to(compatm   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), fldname, rc=rc)) then
             call addmap_from(compice, fldname, compatm, maptype, 'ifrac', 'unset')
             call addmrg_to(compatm, fldname, mrg_from=compice, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to atm: unmerged sea ice albedo, 4 bands
    allocate(flds(4))
    flds = (/'Si_avsdr', 'Si_avsdf', 'Si_anidr', 'Si_anidf'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
             call addfld_from(compice , fldname)
             call addfld_to(compatm   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), fldname, rc=rc)) then
             call addmap_from(compice, fldname, compatm, maptype, 'ifrac', 'unset')
             call addmrg_to(compatm, fldname, mrg_from=compice, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to atm: unmerged surface temperatures from ocn
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compocn , 'So_t')
          call addfld_to(compatm   , 'So_t')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap_from(compocn, 'So_t', compatm, maptype, 'ofrac', 'unset')
          call addmrg_to(compatm, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! to atm: unmerged flux components from lnd
    if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compatm)) then
       allocate(flds(6))
       flds = (/ 'lat ', 'sen ', 'evap', 'gflx', 'roff', 'soff' /)
       if (phase == 'advertise') then
          do n = 1,size(flds)
             call addfld_from(complnd, 'Fall_'//trim(flds(n)))
             call addfld_to(compatm, 'Fall_'//trim(flds(n)))
          end do
       else
          do n = 1,size(flds)
             if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Fall_'//trim(flds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_'//trim(flds(n)), rc=rc)) then
                call addmap_from(complnd, 'Fall_'//trim(flds(n)), compatm, maptype, 'lfrac', 'unset')
                call addmrg_to(compatm, 'Fall_'//trim(flds(n)), mrg_from=complnd, mrg_fld='Fall_'//trim(flds(n)), mrg_type='copy')
             end if
          end do
       end if
       deallocate(flds)
    end if

    ! to atm: unmerged state variables from lnd
    if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compatm)) then
       allocate(flds(7))
       flds = (/ 'sfrac', 'tref ', 'qref ', 'q    ', 'cmm  ', 'chh  ', 'zvfun' /)
       if (phase == 'advertise') then
          do n = 1,size(flds)
             call addfld_from(complnd, 'Sl_'//trim(flds(n)))
             call addfld_to(compatm, 'Sl_'//trim(flds(n)))
          end do
       else
          do n = 1,size(flds)
             if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sl_'//trim(flds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_'//trim(flds(n)), rc=rc)) then
                call addmap_from(complnd, 'Sl_'//trim(flds(n)), compatm, maptype, 'lfrac', 'unset')
                call addmrg_to(compatm, 'Sl_'//trim(flds(n)), mrg_from=complnd, mrg_fld='Sl_'//trim(flds(n)), mrg_type='copy')
             end if
          end do
       end if
       deallocate(flds)
    end if

    ! to atm: unmerged from mediator, merge will be done under FV3/CCPP composite step
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux,
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux
    allocate(flds(5))
    flds = (/ 'Faox_lat ', 'Faox_sen ', 'Faox_lwup', 'Faox_taux', 'Faox_tauy' /)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_aoflux(fldname)
             call addfld_to(compatm , fldname)
          end if
       else
          if (fldchk(is_local%wrap%FBMed_aoflux_o, fldname, rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux(fldname, compatm, maptype, 'ofrac', 'unset')
             end if
             call addmrg_to(compatm, fldname, mrg_from=compmed, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to atm: surface roughness length from wav
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compwav , 'Sw_z0')
          call addfld_to(compatm   , 'Sw_z0')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sw_z0', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_z0', rc=rc)) then
          call addmap_from(compwav, 'Sw_z0', compatm, mapbilnr_nstod, 'one', 'unset')
          call addmrg_to(compatm, 'Sw_z0', mrg_from=compwav, mrg_fld='Sw_z0', mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! to ocn: sea level pressure from atm
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
          call addfld_from(compatm , 'Sa_pslv')
          call addfld_to(compocn   , 'Sa_pslv')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Sa_pslv', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_pslv', rc=rc)) then
          call addmap_from(compatm, 'Sa_pslv', compocn, maptype, 'one', 'unset')
          call addmrg_to(compocn, 'Sa_pslv', mrg_from=compatm, mrg_fld='Sa_pslv', mrg_type='copy')
       end if
    end if

    ! to ocn: swpen thru ice w/o bands
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
          call addfld_from(compice , 'Fioi_swpen')
       end if
    else
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen', rc=rc)) then
          call addmap_from(compice, 'Fioi_swpen', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if
    ! to ocn: from sw from atm and sw net from ice (custom merge in med_phases_prep_ocn)
    ! - downward direct  near-infrared ("n" or "i") incident solar radiation
    ! - downward diffuse near-infrared ("n" or "i") incident solar radiation
    ! - downward direct visible ("v") incident solar radiation
    ! - downward diffuse visible ("v") incident solar radiation
    allocate(oflds(4))
    allocate(aflds(4))
    allocate(iflds(4))
    oflds = (/'Foxx_swnet_idr', 'Foxx_swnet_idf', 'Foxx_swnet_vdr', 'Foxx_swnet_vdf'/)
    aflds = (/'Faxa_swndr'    , 'Faxa_swndf'    , 'Faxa_swvdr'    , 'Faxa_swvdf'/)
    iflds = (/'Fioi_swpen_idr', 'Fioi_swpen_idf', 'Fioi_swpen_vdr', 'Fioi_swpen_vdf'/)
    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compatm , trim(aflds(n)))
             call addfld_to(compocn   , trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap_from(compatm, trim(aflds(n)), compocn, maptype, 'one', 'unset')
          end if
       end if
    end do

    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compice , trim(iflds(n)))
             call addfld_to(compocn   , trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(iflds(n)), rc=rc)) then
             call addmap_from(compice, trim(iflds(n)), compocn, mapfcopy, 'unset', 'unset')
          end if
       end if
    end do
    deallocate(oflds)
    deallocate(aflds)
    deallocate(iflds)

    ! to ocn: rain and snow via auto merge
    allocate(flds(2))
    flds = (/'Faxa_rain', 'Faxa_snow'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compatm , fldname)
             call addfld_to(compocn   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
             call addmap_from(compatm, fldname, compocn, maptype, 'one', 'unset')
             call addmrg_to(compocn, fldname, &
                  mrg_from=compatm, mrg_fld=fldname, mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end do

    !to ocn: surface stress from mediator or atm and ice stress via auto merge
    flds = (/'taux', 'tauy'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_aoflux('Faox_'//fldname)
             call addfld_from(compatm , 'Faxa_'//fldname)
             call addfld_from(compice , 'Fioi_'//fldname)
             call addfld_to(compocn   , 'Foxx_'//fldname)
          end if
       else
          if (med_aoflux_to_ocn) then
             if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_'//fldname, rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_'//fldname, rc=rc) .and. &
                  fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_'//fldname, rc=rc)) then
                call addmap_from(compice,  'Fioi_'//fldname, compocn, mapfcopy, 'unset', 'unset')
                call addmrg_to(compocn,  'Foxx_'//fldname, &
                     mrg_from=compmed, mrg_fld='Faox_'//fldname, mrg_type='merge', mrg_fracname='ofrac')
                call addmrg_to(compocn,  'Foxx_'//fldname, &
                     mrg_from=compice, mrg_fld='Fioi_'//fldname, mrg_type='merge', mrg_fracname='ifrac')
             end if
          else
             if ( fldchk(is_local%wrap%FBexp(compocn)   , 'Foxx_'//fldname, rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_'//fldname, rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_'//fldname, rc=rc)) then
                call addmap_from(compice, 'Fioi_'//fldname, compocn, mapfcopy, 'unset', 'unset')
                call addmap_from(compatm, 'Faxa_'//fldname, compocn, mapconsf_aofrac, 'aofrac', 'unset')
                call addmrg_to(compocn, 'Foxx_'//fldname, &
                     mrg_from=compice, mrg_fld='Fioi_'//fldname, mrg_type='merge', mrg_fracname='ifrac')
                call addmrg_to(compocn, 'Foxx_'//fldname, &
                     mrg_from=compatm, mrg_fld='Faxa_'//fldname, mrg_type='merge', mrg_fracname='ofrac')
             end if
          end if
       end if
       end do
    deallocate(flds)

    ! to ocn: net long wave via auto merge
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
          call addfld_aoflux('Faox_lwup')
          call addfld_from(compatm , 'Faxa_lwnet')
          call addfld_from(compatm , 'Faxa_lwdn')
          call addfld_to(compocn   , 'Foxx_lwnet')
       end if
    else
       if (med_aoflux_to_ocn) then
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc)) then
             call addmap_from(compatm, 'Faxa_lwdn', compocn, maptype, 'one', 'unset')
             call addmrg_to(compocn, 'Foxx_lwnet', &
                  mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
             call addmrg_to(compocn, 'Foxx_lwnet', &
                  mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='merge', mrg_fracname='ofrac')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwnet', rc=rc)) then
             call addmap_from(compatm, 'Faxa_lwnet', compocn, mapconsf_aofrac, 'aofrac', 'unset')
             call addmrg_to(compocn, 'Foxx_lwnet', &
                  mrg_from=compatm, mrg_fld='Faxa_lwnet', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! to ocn: sensible heat flux
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
          call addfld_aoflux('Faox_sen')
          call addfld_from(compatm , 'Faxa_sen')
          call addfld_to(compocn   , 'Foxx_sen')
       end if
    else
       if (med_aoflux_to_ocn) then
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_sen', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_sen' , rc=rc)) then
             call addmrg_to(compocn, 'Foxx_sen', &
                  mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_sen', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_sen', rc=rc)) then
             call addmap_from(compatm, 'Faxa_sen', compocn, mapconsf_aofrac, 'aofrac', 'unset')
             call addmrg_to(compocn, 'Foxx_sen', &
                  mrg_from=compatm, mrg_fld='Faxa_sen', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! to ocn: evaporation water flux
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
          call addfld_aoflux('Faox_evap')
          call addfld_from(compatm , 'Faxa_evap')
          call addfld_to(compocn   , 'Foxx_evap')
       end if
    else
       if (med_aoflux_to_ocn) then
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap' , rc=rc)) then
             call addmrg_to(compocn, 'Foxx_evap', &
                  mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_evap' , rc=rc)) then
             call addmap_from(compatm, 'Faxa_evap', compocn, mapconsf_aofrac, 'aofrac', 'unset')
             call addmrg_to(compocn, 'Foxx_evap', &
                  mrg_from=compatm, mrg_fld='Faxa_evap', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! to ocn: unmerged fluxes from ice
    ! - water flux due to melting ice from ice
    ! - heat flux from melting ice from ice
    ! - salt flux from ice
    allocate(flds(3))
    flds = (/'Fioi_meltw', 'Fioi_melth', 'Fioi_salt '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compice , fldname)
             call addfld_to(compocn   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), fldname, rc=rc)) then
             call addmap_from(compice, fldname, compocn,  mapfcopy, 'unset', 'unset')
             call addmrg_to(compocn, fldname, &
                  mrg_from=compice, mrg_fld=fldname, mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end do
    deallocate(flds)

    ! to ocn: partitioned stokes drift from wav
    allocate(flds(2))
    flds = (/'Sw_pstokes_x', 'Sw_pstokes_y'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compwav , fldname)
             call addfld_to(compocn   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav,compwav), fldname, rc=rc)) then
             call addmap_from(compwav, fldname, compocn, mapbilnr_nstod, 'one', 'unset')
             call addmrg_to(compocn, fldname, mrg_from=compwav, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! to ice: fluxes from atm
    ! - downward longwave heat flux from atm
    ! - downward direct near-infrared incident solar radiation  from atm
    ! - downward direct visible incident solar radiation        from atm
    ! - downward diffuse near-infrared incident solar radiation from atm
    ! - downward Diffuse visible incident solar radiation       from atm
    ! - rain from atm
    ! - snow from atm

    allocate(flds(7))
    flds = (/'Faxa_lwdn ', 'Faxa_swndr', 'Faxa_swvdr', 'Faxa_swndf', 'Faxa_swvdf', &
             'Faxa_rain ', 'Faxa_snow '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compice)) then
             call addfld_from(compatm , fldname)
             call addfld_to(compice   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
             call addmap_from(compatm, fldname, compice, maptype, 'one', 'unset')
             call addmrg_to(compice, fldname, mrg_from=compatm, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to ice: states from atm
    ! - height at the lowest model level from atm
    ! - pressure at the lowest model level from atm
    ! - temperature at the lowest model level from atm
    ! - zonal wind at the lowest model level from atm
    ! - meridional wind at the lowest model level from atm
    ! - specific humidity at the lowest model level from atm
    allocate(flds(6))
    flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', 'Sa_shum'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compice)) then
             call addfld_from(compatm , fldname)
             call addfld_to(compice   , fldname)
          endif
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
             call addmap_from(compatm, fldname, compice, maptype, 'one', 'unset')
             call addmrg_to(compice, fldname, mrg_from=compatm, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to ice: states and fluxes from ocn
    ! - sea surface temperature from ocn
    ! - sea surface salinity from ocn
    ! - zonal sea water velocity from ocn
    ! - meridional sea water velocity from ocn
    ! - zonal sea surface slope from ocn
    ! - meridional sea surface slope from ocn
    ! - ocean melt and freeze potential from ocn
    allocate(flds(7))
    flds = (/'So_t   ', 'So_s   ', 'So_u   ', 'So_v   ','So_dhdx', 'So_dhdy', 'Fioo_q '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compice)) then
             call addfld_from(compocn , fldname)
             call addfld_to(compice   , fldname)
          endif
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compocn,compocn), fldname, rc=rc)) then
             call addmap_from(compocn, fldname, compice, mapfcopy , 'unset', 'unset')
             call addmrg_to(compice, fldname, mrg_from=compocn, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compwav)) then
          call addfld_from(compwav , 'Sw_elevation_spectrum')
          call addfld_to(compice   , 'Sw_elevation_spectrum')
       end if
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Sw_elevation_spectrum', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_elevation_spectrum', rc=rc)) then
          call addmap_from(compwav, 'Sw_elevation_spectrum', compice, mapbilnr_nstod, 'one', 'unset')
          call addmrg_to(compice, 'Sw_elevation_spectrum', mrg_from=compwav, &
               mrg_fld='Sw_elevation_spectrum', mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO WAV (compwav)
    !=====================================================================

    ! to wav: states from atm
    ! - 10m meridonal and zonal winds
    ! - bottom temperature from atm
    allocate(flds(3))
    flds = (/'Sa_u10m', 'Sa_v10m', 'Sa_tbot'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compwav)) then
             call addfld_from(compatm , fldname)
             call addfld_to(compwav   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
             call addmap_from(compatm, fldname, compwav, mapbilnr_nstod, 'one', 'unset')
             call addmrg_to(compwav, fldname, mrg_from=compatm, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
     end do
     deallocate(flds)

     ! to wav: states from ice
     ! - sea ice fraction
     ! - sea ice thickness
     ! - sea ice floe diameter
     allocate(flds(3))
     flds = (/'Si_ifrac   ', 'Si_floediam', 'Si_thick   '/)
     do n = 1,size(flds)
        fldname = trim(flds(n))
        if (phase == 'advertise') then
           if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compwav)) then
              call addfld_from(compice , fldname)
              call addfld_to(compwav   , fldname)
        end if
        else
           if ( fldchk(is_local%wrap%FBexp(compwav)        , fldname, rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compice,compice), fldname, rc=rc)) then
              call addmap_from(compice, fldname, compwav, mapbilnr_nstod , 'one', 'unset')
              call addmrg_to(compwav, fldname, mrg_from=compice, mrg_fld=fldname, mrg_type='copy')
           end if
        end if
      end do
      deallocate(flds)

      ! to wav: states from ocn
      ! - zonal sea water velocity from ocn
      ! - meridional sea water velocity from ocn
      ! - surface temperature from ocn
      allocate(flds(3))
      flds = (/'So_u', 'So_v', 'So_t'/)
      do n = 1,size(flds)
         fldname = trim(flds(n))
         if (phase == 'advertise') then
            if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compwav)) then
               call addfld_from(compocn , fldname)
               call addfld_to(compwav   , fldname)
            end if
         else
            if ( fldchk(is_local%wrap%FBexp(compwav)        , fldname, rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compocn,compocn), fldname, rc=rc)) then
               call addmap_from(compocn, fldname, compwav, mapbilnr_nstod , 'one', 'unset')
              call addmrg_to(compwav, fldname, mrg_from=compocn, mrg_fld=fldname, mrg_type='copy')
           end if
        end if
     end do
     deallocate(flds)

    !=====================================================================
    ! FIELDS TO LAND (complnd)
    !=====================================================================

    ! to lnd - states and fluxes from atm
    if ( trim(coupling_mode) == 'ufs.nfrac.aoflux') then
       allocate(flds(21))
       flds = (/'Sa_z      ', 'Sa_topo   ', 'Sa_tbot   ', 'Sa_pbot   ', &
                'Sa_shum   ', 'Sa_u      ', 'Sa_v      ', 'Sa_pslv   ', &
                'Faxa_lwdn ', 'Faxa_swdn ', 'Faxa_snowc', 'Faxa_snowl', &
                'Faxa_rainc', 'Faxa_rainl', 'Faxa_rain ', 'Faxa_swnet'/)
    else
       allocate(flds(18))
       flds = (/'Sa_z      ', 'Sa_ta     ', 'Sa_pslv   ', 'Sa_qa     ', &
                'Sa_u      ', 'Sa_v      ', 'Faxa_swdn ', 'Faxa_lwdn ', &
                'Faxa_swnet', 'Faxa_rain ', 'Sa_prsl   ', 'Sa_vfrac  ', &
                'Faxa_snow ', 'Faxa_rainc', 'Sa_tskn   ', 'Sa_exner  ', &
                'Sa_ustar  ', 'Sa_zorl   ' /)
    end if
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(complnd)) then
             call addfld_from(compatm , fldname)
             call addfld_to(complnd   , fldname)
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)        , fldname, rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), fldname, rc=rc)) then
             call addmap_from(compatm, fldname, complnd, maptype, 'one', 'unset')
             call addmrg_to(complnd, fldname, mrg_from=compatm, mrg_fld=fldname, mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

  end subroutine esmFldsExchange_ufs

end module esmFldsExchange_ufs_mod
