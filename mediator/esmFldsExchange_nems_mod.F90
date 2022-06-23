module esmFldsExchange_nems_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_nems

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange_nems(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use med_internalstate_mod , only : mastertask, logunit
    use med_internalstate_mod , only : compmed, compatm, compocn, compice, complnd, compwav, ncomps
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use med_internalstate_mod , only : mapconsf_aofrac, mapbilnr_nstod
    use med_internalstate_mod , only : coupling_mode, mapnames
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use med_internalstate_mod , only : InternalState, mastertask, logunit

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: i, n, maptype
    character(len=CX)   :: msgString
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS), allocatable :: flds(:), oflds(:), aflds(:), iflds(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_nems)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set maptype according to coupling_mode
    if (trim(coupling_mode) == 'nems_orig' .or. trim(coupling_mode) == 'nems_orig_data') then
      maptype = mapnstod_consf
    else
      maptype = mapconsf
    end if
    write(msgString,'(A,i6,A)') trim(subname)//': maptype is ',maptype,', '//mapnames(maptype)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    !=====================================================================
    ! scalar information
    !=====================================================================

    if (phase == 'advertise') then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld(fldListFr(n)%flds, trim(cvalue))
          call addfld(fldListTo(n)%flds, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    ! masks from components
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice)) call addfld(fldListFr(compice)%flds, 'Si_imask')
       if (is_local%wrap%comp_present(compocn)) call addfld(fldListFr(compocn)%flds, 'So_omask')
       if (is_local%wrap%comp_present(complnd)) call addfld(fldListFr(complnd)%flds, 'Sl_lfrin')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
       end if
    end if

    if ( trim(coupling_mode) == 'nems_orig_data') then
       ! atm fields required for atm/ocn flux calculation
       allocate(flds(10))
       flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', &
                'Sa_shum', 'Sa_u10m', 'Sa_v10m', 'Sa_t2m ', 'Sa_q2m '/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             if (is_local%wrap%comp_present(compatm) )then
                call addfld(fldListFr(compatm)%flds, trim(fldname))
             end if
          else
            if ( fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
               call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
            end if
          end if
       end do
       deallocate(flds)

       ! fields returned by the atm/ocn flux computation which are otherwise unadvertised
       allocate(flds(8))
       flds = (/'So_tref  ', 'So_qref  ', 'So_ustar ', 'So_re    ','So_ssq   ', &
                'So_u10   ', 'So_duu10n', 'Faox_lat '/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             call addfld(fldListMed_aoflux%flds, trim(fldname))
          end if
       end do
       deallocate(flds)
    end if

    if (trim(coupling_mode) == 'nems_frac_aoflux' .or. trim(coupling_mode) == 'nems_frac_aoflux_sbs') then
       allocate(flds(12))
       flds = (/'Sa_u     ', 'Sa_v     ', 'Sa_z     ', 'Sa_tbot  ', 'Sa_pbot  ', &
                'Sa_pslv  ', 'Sa_shum  ', 'Sa_ptem  ', 'Sa_dens  ', 'Sa_u10m  ', &
                'Sa_v10m  ', 'Faxa_lwdn'/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             if (is_local%wrap%comp_present(compatm) )then
                call addfld(fldListFr(compatm)%flds, trim(fldname))
             end if
          else
            if ( fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
               call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
            end if
          end if
       end do
       deallocate(flds)

       ! fields returned by the atm/ocn flux computation which are otherwise unadvertised
       allocate(flds(13))
       flds = (/'So_tref  ', 'So_qref  ','So_u10   ', 'So_ustar ','So_ssq   ', &
                'So_re    ', 'So_duu10n','Faox_lwup', 'Faox_sen ','Faox_lat ', &
                'Faox_evap', 'Faox_taux','Faox_tauy'/)
       do n = 1,size(flds)
          fldname = trim(flds(n))
          if (phase == 'advertise') then
             call addfld(fldListMed_aoflux%flds, trim(fldname))
          end if
       end do
       deallocate(flds)
    end if

    ! TODO: unused, but required to maintain B4B repro for mediator restarts; should be removed
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'mean_sw_pen_to_ocn')
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE (compatm)
    !=====================================================================

    ! to atm: fractions (computed in med_phases_prep_atm)
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compice)%flds, 'Si_ifrac')
          call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
       end if
       ! ofrac used by atm
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_ofrac')
       end if
       ! lfrac used by atm
       if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListTo(compatm)%flds, 'Sl_lfrac')
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
    flds = (/'Faii_taux', 'Faii_tauy', 'Faii_lat ', 'Faii_sen ', 'Faii_lwup', &
             'Faii_evap', 'Si_vice  ', 'Si_vsno  ', 'Si_t     '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
             call addfld(fldListFr(compice)%flds, trim(fldname))
             call addfld(fldListTo(compatm)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compatm, maptype, 'ifrac', 'unset')
             call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    allocate(flds(4))
    flds = (/'Si_avsdr', 'Si_avsdf', 'Si_anidr', 'Si_anidf'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm)) then
             call addfld(fldListFr(compice)%flds, trim(fldname))
             call addfld(fldListTo(compatm)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compatm, maptype, 'ifrac', 'unset')
             call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to atm: unmerged surface temperatures from ocn
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compocn)%flds, 'So_t')
          call addfld(fldListTo(compatm)%flds, 'So_t')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, maptype, 'ofrac', 'unset')
          call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! to atm: unmerged surface temperatures from lnd
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(complnd)%flds, 'Sl_t')
          call addfld(fldListTo(compatm)%flds, 'Sl_t')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sl_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_t', compatm, maptype, 'lfrin', 'unset')
          call addmrg(fldListTo(compatm)%flds, 'Sl_t', mrg_from=complnd, mrg_fld='Sl_t', mrg_type='copy')
       end if
    end if

    ! to atm: unmerged from mediator, merge will be done under FV3/CCPP composite step
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux,
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux
    ! - evaporation water flux from water, not in the list do we need to send it to atm?
    if (trim(coupling_mode) == 'nems_frac_aoflux') then
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          allocate(flds(5))
          flds = (/ 'lat ', 'sen ', 'lwup', 'taux', 'tauy' /)
          if (phase == 'advertise') then
             do n = 1,size(flds)
                call addfld(fldListMed_aoflux%flds , 'Faox_'//trim(flds(n)))
                call addfld(fldListTo(compatm)%flds, 'Faox_'//trim(flds(n)))
             end do
          else
             do n = 1,size(flds)
                if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_'//trim(flds(n)), rc=rc)) then
                   if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                      call addmap(fldListMed_aoflux%flds, 'Faox_'//trim(flds(n)), compatm, maptype, 'ofrac', 'unset')
                   end if
                   call addmrg(fldListTo(compatm)%flds, 'Faox_'//trim(flds(n)), mrg_from=compmed, mrg_fld='Faox_'//trim(flds(n)), mrg_type='copy')
                end if
             end do
          end if
          deallocate(flds)
       end if
    end if

    ! to atm: surface roughness length from wav
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compwav)%flds, 'Sw_z0')
          call addfld(fldListTo(compatm)%flds, 'Sw_z0')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sw_z0', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_z0', rc=rc)) then
          call addmap(fldListFr(compwav)%flds, 'Sw_z0', compatm, mapnstod_consf, 'one', 'unset')
          call addmrg(fldListTo(compatm)%flds, 'Sw_z0', mrg_from=compwav, mrg_fld='Sw_z0', mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! to ocn: sea level pressure from atm
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
          call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Sa_pslv', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_pslv', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, maptype, 'one', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', mrg_from=compatm, mrg_fld='Sa_pslv', mrg_type='copy')
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
             call addfld(fldListFr(compatm)%flds, trim(aflds(n)))
             call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'one', 'unset')
          end if
       end if
    end do

    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compice)%flds, trim(iflds(n)))
             call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(iflds(n)), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(iflds(n)), compocn, mapfcopy, 'unset', 'unset')
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
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(compocn)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end do
    deallocate(flds)

    if (trim(coupling_mode) == 'nems_orig' .or. trim(coupling_mode) == 'nems_frac' .or. &
        trim(coupling_mode) == 'nems_frac_aoflux_sbs') then
       ! to ocn: merge surface stress (custom merge calculation in med_phases_prep_ocn)
       allocate(oflds(2))
       allocate(aflds(2))
       allocate(iflds(2))
       oflds = (/'Foxx_taux', 'Foxx_tauy'/)
       aflds = (/'Faxa_taux', 'Faxa_tauy'/)
       iflds = (/'Fioi_taux', 'Fioi_tauy'/)
       do n = 1,size(oflds)
          if (phase == 'advertise') then
             if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compatm) &
                .and. is_local%wrap%comp_present(compocn)) then
                   call addfld(fldListFr(compice)%flds, trim(iflds(n)))
                   call addfld(fldListFr(compatm)%flds, trim(aflds(n)))
                   call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
             end if
          else
             if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compice,compice), trim(iflds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
                call addmap(fldListFr(compice)%flds, trim(iflds(n)), compocn, mapfcopy, 'unset', 'unset')
                call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, mapconsf_aofrac, 'aofrac', 'unset')
             end if
          end if
       end do
       deallocate(oflds)
       deallocate(aflds)
       deallocate(iflds)

       ! to ocn: net long wave via auto merge
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compatm)%flds, 'Faxa_lwnet')
             call addfld(fldListTo(compocn)%flds, 'Faxa_lwnet')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Faxa_lwnet', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwnet', rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapconsf_aofrac, 'aofrac', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Faxa_lwnet', &
                  mrg_from=compatm, mrg_fld='Faxa_lwnet', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if

       ! to ocn: merged sensible heat flux (custom merge in med_phases_prep_ocn)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compatm)%flds, 'Faxa_sen')
             call addfld(fldListTo(compocn)%flds, 'Faxa_sen')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Faxa_sen', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_sen', rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapconsf_aofrac, 'aofrac', 'unset')
          end if
       end if

       ! to ocn: evaporation water flux (custom merge in med_phases_prep_ocn)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compatm)%flds, 'Faxa_lat')
             call addfld(fldListTo(compocn)%flds, 'Faxa_evap')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Faxa_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lat' , rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, mapconsf_aofrac, 'aofrac', 'unset')
          end if
       end if
    else if (trim(coupling_mode) == 'nems_orig_data' .or. trim(coupling_mode) == 'nems_frac_aoflux') then
       ! nems_orig_data
       ! to ocn: surface stress from mediator and ice stress via auto merge
       allocate(flds(2))
       flds = (/'taux', 'tauy'/)
       do n = 1,size(flds)
          if (phase == 'advertise') then
             if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
                call addfld(fldListMed_aoflux%flds  , 'Faox_'//trim(flds(n)))
                call addfld(fldListFr(compice)%flds , 'Fioi_'//trim(flds(n)))
                call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(flds(n)))
             end if
          else
             if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_'//trim(flds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_'//trim(flds(n)), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_'//trim(flds(n)), rc=rc)) then
                call addmap(fldListFr(compice)%flds,  'Fioi_'//trim(flds(n)), compocn, mapfcopy, 'unset', 'unset')
                call addmrg(fldListTo(compocn)%flds,  'Foxx_'//trim(flds(n)), &
                     mrg_from=compmed, mrg_fld='Faox_'//trim(flds(n)), mrg_type='merge', mrg_fracname='ofrac')
                call addmrg(fldListTo(compocn)%flds,  'Foxx_'//trim(flds(n)), &
                     mrg_from=compice, mrg_fld='Fioi_'//trim(flds(n)), mrg_type='merge', mrg_fracname='ifrac')
             end if
          end if
       end do
       deallocate(flds)

       ! to ocn: long wave net via auto merge
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListMed_aoflux%flds , 'Faox_lwup')
             call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
             call addfld(fldListTo(compocn)%flds, 'Foxx_lwnet')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, maptype, 'one', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
                  mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
             call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
                  mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if

       ! to ocn: sensible heat flux from mediator via auto merge
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compocn)) then
             call addfld(fldListMed_aoflux%flds , 'Faox_sen')
             call addfld(fldListTo(compocn)%flds, 'Faox_sen')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Faox_sen', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_sen' , rc=rc)) then
             call addmrg(fldListTo(compocn)%flds, 'Faox_sen', &
                  mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if

       ! to ocn: evaporation water flux from mediator via auto merge
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compocn)) then
             call addfld(fldListMed_aoflux%flds , 'Faox_evap')
             call addfld(fldListTo(compocn)%flds, 'Faox_evap')
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Faox_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap' , rc=rc)) then
             call addmrg(fldListTo(compocn)%flds, 'Faox_evap', &
                  mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! to ocn: water flux due to melting ice from ice
    ! to ocn: heat flux from melting ice from ice
    ! to ocn: salt flux from ice
    allocate(flds(3))
    flds = (/'Fioi_meltw', 'Fioi_melth', 'Fioi_salt '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compice)%flds, trim(fldname))
             call addfld(fldListTo(compocn)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end do
    deallocate(flds)

    ! to ocn: partitioned stokes drift from wav
    allocate(flds(6))
    flds = (/'Sw_ustokes1', 'Sw_ustokes2', 'Sw_ustokes3', &
             'Sw_vstokes1', 'Sw_vstokes2', 'Sw_vstokes3'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compwav)%flds, trim(fldname))
             call addfld(fldListTo(compocn)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav,compwav), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compwav)%flds, trim(fldname), compocn, mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from=compwav, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! to ice - fluxes from atm
    ! to ice: downward longwave heat flux from atm
    ! to ice: downward direct near-infrared incident solar radiation  from atm
    ! to ice: downward direct visible incident solar radiation        from atm
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! to ice: downward Diffuse visible incident solar radiation       from atm
    ! to ice: rain from atm
    ! to ice: snow from atm

    allocate(flds(7))
    flds = (/'Faxa_lwdn ', 'Faxa_swndr', 'Faxa_swvdr', 'Faxa_swndf', 'Faxa_swvdf', &
             'Faxa_rain ', 'Faxa_snow '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compice)) then
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(compice)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, maptype, 'one', 'unset')
             call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to ice - state from atm
    ! to ice: height at the lowest model level from atm
    ! to ice: pressure at the lowest model level from atm
    ! to ice: temperature at the lowest model level from atm
    ! to ice: zonal wind at the lowest model level from atm
    ! to ice: meridional wind at the lowest model level from atm
    ! to ice: specific humidity at the lowest model level from atm
    allocate(flds(6))
    flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', &
             'Sa_shum'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compice)) then
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(compice)%flds, trim(fldname))
          endif
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, maptype, 'one', 'unset')
             call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! to ice - states and fluxes from ocn
    ! to ice: sea surface temperature from ocn
    ! to ice: sea surface salinity from ocn
    ! to ice: zonal sea water velocity from ocn
    ! to ice: meridional sea water velocity from ocn
    ! to ice: zonal sea surface slope from ocn
    ! to ice: meridional sea surface slope from ocn
    ! to ice: ocean melt and freeze potential from ocn
    allocate(flds(7))
    flds = (/'So_t   ', 'So_s   ', 'So_u   ', 'So_v   ','So_dhdx', &
             'So_dhdy', 'Fioo_q '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compice)) then
             call addfld(fldListFr(compocn)%flds, trim(fldname))
             call addfld(fldListTo(compice)%flds, trim(fldname))
          endif
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy , 'unset', 'unset')
             call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO WAV (compwav)
    !=====================================================================

    ! to wav - 10m winds and bottom temperature from atm
    allocate(flds(3))
    flds = (/'Sa_u10m', 'Sa_v10m', 'Sa_tbot'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compwav)) then
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(compwav)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compwav, mapnstod_consf, 'one', 'unset')
             call addmrg(fldListTo(compwav)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
     end do
     deallocate(flds)

     ! to wav: sea ice fraction
     if (phase == 'advertise') then
        if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compwav)) then
           call addfld(fldListFr(compice)%flds, 'Si_ifrac')
           call addfld(fldListTo(compwav)%flds, 'Si_ifrac')
        end if
     else
        if ( fldchk(is_local%wrap%FBexp(compwav)        , 'Si_ifrac', rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compice,compice), 'Si_ifrac', rc=rc)) then
            call addmap(fldListFr(compice)%flds, 'Si_ifrac', compwav, mapfcopy , 'unset', 'unset')
            call addmrg(fldListTo(compwav)%flds, 'Si_ifrac', mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')
        end if
     end if

     ! to wav: zonal sea water velocity from ocn
     ! to wav: meridional sea water velocity from ocn
     ! to wav: surface temperature from ocn
     allocate(flds(3))
     flds = (/'So_u', 'So_v', 'So_t'/)
     do n = 1,size(flds)
        fldname = trim(flds(n))
        if (phase == 'advertise') then
           if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compwav)) then
              call addfld(fldListFr(compocn)%flds, trim(fldname))
              call addfld(fldListTo(compwav)%flds, trim(fldname))
           end if
        else
           if ( fldchk(is_local%wrap%FBexp(compwav)        , trim(fldname), rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
              call addmap(fldListFr(compocn)%flds, trim(fldname), compwav, mapfcopy , 'unset', 'unset')
              call addmrg(fldListTo(compwav)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
           end if
        end if
     end do
     deallocate(flds)

    !=====================================================================
    ! FIELDS TO LAND (complnd)
    !=====================================================================

    ! to lnd - states and fluxes from atm
    if ( trim(coupling_mode) == 'nems_orig_data') then
       allocate(flds(21))
       flds = (/'Sa_z      ', 'Sa_topo   ', 'Sa_tbot   ', 'Sa_pbot   ', &
                'Sa_shum   ', 'Sa_u      ', 'Sa_v      ', 'Faxa_lwdn ', &
                'Sa_ptem   ', 'Sa_dens   ', 'Faxa_swdn ', 'Faxa_swnet', &
                'Faxa_snowc', 'Faxa_snowl', 'Faxa_rainc', 'Faxa_rainl', & 
                'Sa_pslv   ', &
                'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swvdf'/)
    else
       allocate(flds(9))
       flds = (/'Sa_z      ', 'Sa_tbot   ', 'Sa_pbot   ', 'Sa_shum   ', &
                'Sa_u      ', 'Sa_v      ', 'Faxa_swdn ', 'Faxa_lwdn ', &
                'Faxa_rain ' /)
    end if
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(complnd)) then
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(complnd)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             print*, "i am here !!!"
             call addmap(fldListFr(compatm)%flds, trim(fldname), complnd, maptype, 'one', 'unset')
             call addmrg(fldListTo(complnd)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

  end subroutine esmFldsExchange_nems

end module esmFldsExchange_nems_mod
