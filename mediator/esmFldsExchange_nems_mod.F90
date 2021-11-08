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
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : compmed, compatm, compocn, compice, comprof, ncomps
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use esmflds               , only : mapconsf_aofrac
    use esmflds               , only : coupling_mode, mapnames
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use med_internalstate_mod , only : mastertask, logunit

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    integer             :: i, n, maptype
    character(len=CX)   :: msgString
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS), allocatable :: flds(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_nems)'
    !--------------------------------------

    rc = ESMF_SUCCESS

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

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,ncomps
       call addfld(fldListFr(n)%flds, trim(cvalue))
       call addfld(fldListTo(n)%flds, trim(cvalue))
    end do

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    ! masks from components
    call addfld(fldListFr(compice)%flds, 'Si_imask')
    call addfld(fldListFr(compocn)%flds, 'So_omask')
    call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')

    if ( trim(coupling_mode) == 'nems_orig_data') then
      ! atm and ocn fields required for atm/ocn flux calculation'
      allocate(flds(10))
      flds = (/'Sa_u   ','Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', 'Sa_shum', &
               'Sa_u10m','Sa_v10m', 'Sa_t2m ', 'Sa_q2m '/)
      do n = 1,size(flds)
         fldname = trim(flds(n))
         call addfld(fldListFr(compatm)%flds, trim(fldname))
         call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
      end do
      deallocate(flds)

      ! unused fields needed by the atm/ocn flux computation
      allocate(flds(13))
      flds = (/'So_tref  ', 'So_qref  ','So_u10   ', 'So_ustar ','So_ssq   ', &
               'So_re    ', 'So_duu10n','Faox_lwup', 'Faox_sen ','Faox_lat ', &
               'Faox_evap', 'Faox_taux','Faox_tauy'/)
      do n = 1,size(flds)
         fldname = trim(flds(n))
         call addfld(fldListMed_aoflux%flds, trim(fldname))
      end do
      deallocate(flds)
    end if

    ! unused fields from ice - but that are needed to be realized by the cice cap
    call addfld(fldListFr(compice)%flds, 'Faii_evap')
    call addfld(fldListFr(compice)%flds, 'mean_sw_pen_to_ocn')

    !=====================================================================
    ! FIELDS TO ATMOSPHERE (compatm)
    !=====================================================================

    ! to atm: fractions (computed in med_phases_prep_atm)
    call addfld(fldListFr(compice)%flds, 'Si_ifrac')
    call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
    ! ofrac used by atm
    call addfld(fldListFr(compatm)%flds, 'Sa_ofrac')

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
    flds = (/'Faii_taux', 'Faii_tauy', 'Faii_lat ', &
             'Faii_sen ', 'Faii_lwup', 'Faii_evap', &
             'Si_vice  ', 'Si_vsno  ', 'Si_t     '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compice)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
       call addmap(fldListFr(compice)%flds, trim(fldname), compatm, maptype, 'ifrac', 'unset')
       call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
    end do
    deallocate(flds)

    allocate(flds(4))
    flds = (/'avsdr    ', 'avsdf    ', &
             'anidr    ', 'anidf    '/)
    do n = 1,size(flds)
       fldname = 'Si_'//trim(flds(n))
       call addfld(fldListFr(compice)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
       call addmap(fldListFr(compice)%flds, trim(fldname), compatm, maptype, 'ifrac', 'unset')
       call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
    end do
    deallocate(flds)

    ! to atm: unmerged surface temperatures from ocn
    call addfld(fldListFr(compocn)%flds, 'So_t')
    call addfld(fldListTo(compatm)%flds, 'So_t')
    call addmap(fldListFr(compocn)%flds, 'So_t', compatm, maptype, 'ofrac', 'unset')
    call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! to ocn: sea level pressure from atm
    call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
    call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
    call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, maptype, 'one', 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', mrg_from=compatm, mrg_fld='Sa_pslv', mrg_type='copy')

    ! to ocn: from atm (custom merge in med_phases_prep_ocn)
    ! - downward direct  near-infrared incident solar radiation
    ! - downward diffuse near-infrared incident solar radiation
    ! - downward dirrect visible incident solar radiation
    ! - downward diffuse visible incident solar radiation
    allocate(flds(4))
    flds = (/'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swvdf'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
    end do
    deallocate(flds)

    ! to ocn: from ice net shortwave radiation (custom merge in med_phases_prep_ocn)
    ! - downward direct  near-infrared incident solar radiation
    ! - downward diffuse near-infrared incident solar radiation
    ! - downward dirrect visible incident solar radiation
    ! - downward diffuse visible incident solar radiation
    allocate(flds(4))
    flds = (/'vdr', 'vdf', 'idr', 'idf'/)
    do n = 1,size(flds)
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_'//trim(flds(n)))
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_'//trim(flds(n)))
       call addmap(fldListFr(compice)%flds, 'Fioi_swpen_'//trim(flds(n)), compocn, mapfcopy, 'unset', 'unset')
    end do
    deallocate(flds)

    ! to ocn: rain and snow via auto merge
    allocate(flds(2))
    flds = (/'Faxa_rain', 'Faxa_snow'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, maptype, 'one', 'unset')
       call addmrg(fldListTo(compocn)%flds, trim(fldname), &
            mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ofrac')
    end do
    deallocate(flds)

    if (trim(coupling_mode) == 'nems_orig' .or. trim(coupling_mode) == 'nems_frac') then
       ! to ocn: merge surface stress (custom merge calculation in med_phases_prep_ocn)
       allocate(flds(2))
       flds = (/'taux', 'tauy'/)
       do n = 1,size(flds)
          call addfld(fldListTo(compocn)%flds, 'Foxx_'//trim(flds(n)))
          call addfld(fldListFr(compice)%flds, 'Fioi_'//trim(flds(n)))
          call addfld(fldListFr(compatm)%flds, 'Faxa_'//trim(flds(n)))
          call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(flds(n)), compocn, mapconsf_aofrac, 'aofrac', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_'//trim(flds(n)), compocn, mapfcopy, 'unset', 'unset')
       end do
       deallocate(flds)

       ! to ocn: net long wave via auto merge
       call addfld(fldListTo(compocn)%flds, 'Faxa_lwnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwnet')
       call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapconsf_aofrac, 'aofrac', 'unset')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_lwnet', &
            mrg_from=compatm, mrg_fld='Faxa_lwnet', mrg_type='copy_with_weights', mrg_fracname='ofrac')

       ! to ocn: merged sensible heat flux (custom merge in med_phases_prep_ocn)
       call addfld(fldListTo(compocn)%flds, 'Faxa_sen')
       call addfld(fldListFr(compatm)%flds, 'Faxa_sen')
       call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapconsf_aofrac, 'aofrac', 'unset')

       ! to ocn: evaporation water flux (custom merge in med_phases_prep_ocn)
       call addfld(fldListTo(compocn)%flds, 'Faxa_evap')
       call addfld(fldListFr(compatm)%flds, 'Faxa_lat')
       call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, mapconsf_aofrac, 'aofrac', 'unset')
    else
       ! nems_orig_data
       ! to ocn: surface stress from mediator and ice stress via auto merge
       allocate(flds(2))
       flds = (/'taux', 'tauy'/)
       do n = 1,size(flds)
          call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(flds(n)))
          call addfld(fldListFr(compice)%flds , 'Fioi_'//trim(flds(n)))
          call addmap(fldListFr(compice)%flds,  'Fioi_'//trim(flds(n)), compocn, mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds,  'Foxx_'//trim(flds(n)), &
               mrg_from=compmed, mrg_fld='Faox_'//trim(flds(n)), mrg_type='merge', mrg_fracname='ofrac')
          call addmrg(fldListTo(compocn)%flds,  'Foxx_'//trim(flds(n)), &
               mrg_from=compice, mrg_fld='Fioi_'//trim(flds(n)), mrg_type='merge', mrg_fracname='ifrac')
       end do
       deallocate(flds)

       ! to ocn: long wave net via auto merge
       call addfld(fldListTo(compocn)%flds, 'Foxx_lwnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
       call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, maptype, 'one', 'unset')
       call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
             mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
       call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
             mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='merge', mrg_fracname='ofrac')

       ! to ocn: sensible heat flux from mediator via auto merge
       call addfld(fldListTo(compocn)%flds, 'Faox_sen')
       call addmrg(fldListTo(compocn)%flds, 'Faox_sen', &
          mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='copy_with_weights', mrg_fracname='ofrac')

       ! to ocn: evaporation water flux from mediator via auto merge
       call addfld(fldListTo(compocn)%flds, 'Faox_evap')
       call addmrg(fldListTo(compocn)%flds, 'Faox_evap', &
          mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='copy_with_weights', mrg_fracname='ofrac')
    end if

    ! to ocn: water flux due to melting ice from ice
    ! to ocn: heat flux from melting ice from ice
    ! to ocn: salt flux from ice
    allocate(flds(3))
    flds = (/'Fioi_meltw', 'Fioi_melth', 'Fioi_salt '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compice)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
       call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
       call addmrg(fldListTo(compocn)%flds, trim(fldname), &
            mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ifrac')
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
    flds = (/'Faxa_lwdn  '    , 'Faxa_swndr '   , 'Faxa_swvdr '   , 'Faxa_swndf ' , 'Faxa_swvdf ', &
             'Faxa_rain  '    , 'Faxa_snow  '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addmap(fldListFr(compatm)%flds, trim(fldname), compice, maptype, 'one', 'unset')
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
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
    flds = (/'Sa_z        ', 'Sa_pbot     ', 'Sa_tbot     ','Sa_u        ','Sa_v        ','Sa_shum     '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addmap(fldListFr(compatm)%flds, trim(fldname), compice, maptype, 'one', 'unset')
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
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
    flds = (/'So_t   ', 'So_s   ', 'So_u   ', 'So_v   ','So_dhdx', 'So_dhdy', 'Fioo_q '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addfld(fldListFr(compocn)%flds, trim(fldname))
       call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy , 'unset', 'unset')
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
    end do
    deallocate(flds)

  end subroutine esmFldsExchange_nems

end module esmFldsExchange_nems_mod
