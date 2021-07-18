module esmFldsExchange_cesm_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !
  ! Merging arguments:
  ! mrg_fromN = source component index that for the field to be merged
  ! mrg_fldN  = souce field name to be merged
  ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum', 'sum_with_weights', 'merge')
  ! NOTE:
  ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn fluxes or for ocn albedos
  !
  ! NOTE:
  ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in the
  ! atm/ocn flux calculations. Input fields required from either the atm or the ocn for
  ! these computation will use the logical 'use_med_aoflux' below. This is used to determine
  ! mappings between the atm and ocn needed for these computations.
  !--------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : logunit, mastertask

  implicit none
  public

  public :: esmFldsExchange_cesm

  ! currently require mapping files 
  character(len=CX)   :: glc2ice_rmap     ='unset'
  character(len=CX)   :: glc2ocn_liq_rmap ='unset'
  character(len=CX)   :: glc2ocn_ice_rmap ='unset'
  character(len=CX)   :: rof2ocn_fmap     ='unset'
  character(len=CX)   :: rof2ocn_ice_rmap ='unset'
  character(len=CX)   :: rof2ocn_liq_rmap ='unset'
  character(len=CX)   :: wav2ocn_smap     ='unset'
  character(len=CX)   :: ice2wav_smap     ='unset'
  character(len=CX)   :: ocn2wav_smap     ='unset'

  ! no mapping files (value is 'idmap' or 'unset')
  character(len=CX)   :: atm2ice_map='unset'
  character(len=CX)   :: atm2ocn_map='unset'
  character(len=CX)   :: atm2lnd_map='unset'
  character(len=CX)   :: ice2atm_map='unset'
  character(len=CX)   :: ocn2atm_map='unset'
  character(len=CX)   :: lnd2atm_map='unset'
  character(len=CX)   :: lnd2rof_map='unset'
  character(len=CX)   :: rof2lnd_map='unset'
  character(len=CX)   :: atm2wav_map='unset'

  logical             :: mapuv_with_cart3d
  logical             :: flds_i2o_per_cat
  logical             :: flds_co2a
  logical             :: flds_co2b
  logical             :: flds_co2c

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange_cesm(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState, logunit, mastertask
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : compmed, compatm, complnd, compocn
    use esmflds               , only : compice, comprof, compwav, ncomps
    use esmflds               , only : compglc, num_icesheets, ocn2glc_coupling ! compglc is an array of integers
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch, mappatch_uv3d
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use esmflds               , only : map_glc2ocn_ice, map_glc2ocn_liq, map_rof2ocn_ice, map_rof2ocn_liq
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use esmFlds               , only : coupling_mode

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: n, ns
    logical             :: is_lnd, is_glc
    character(len=5)    :: iso(2)
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname=' (esmFldsExchange_cesm) '
    !--------------------------------------

    rc = ESMF_SUCCESS

    iso(1) = '     '
    iso(2) = '_wiso'

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    if (phase /= 'advertise') then
       nullify(is_local%wrap)
       call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (phase == 'advertise') then

       ! mapping to atm
       call NUOPC_CompAttributeGet(gcomp, name='ice2atm_map', value=ice2atm_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'ice2atm_map = '// trim(ice2atm_map)
       call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_map', value=lnd2atm_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'lnd2atm_map = '// trim(lnd2atm_map)
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_map', value=ocn2atm_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'ocn2atm_map = '// trim(ocn2atm_map)

       ! mapping to lnd
       call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_map', value=atm2lnd_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'atm2lnd_map = '// trim(atm2lnd_map)
       call NUOPC_CompAttributeGet(gcomp, name='rof2lnd_map', value=rof2lnd_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'rof2lnd_map = '// trim(rof2lnd_map)

       ! mapping to ice
       call NUOPC_CompAttributeGet(gcomp, name='atm2ice_map', value=atm2ice_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'atm2ice_map = '// trim(atm2ice_map)
       call NUOPC_CompAttributeGet(gcomp, name='glc2ice_rmapname', value=glc2ice_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'glc2ice_rmapname = '// trim(glc2ice_rmap)

       ! mapping to ocn
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_map', value=atm2ocn_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'atm2ocn_map = '// trim(atm2ocn_map)
       call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_liq_rmapname', value=glc2ocn_liq_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'glc2ocn_liq_rmapname = '// trim(glc2ocn_liq_rmap)
       call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_ice_rmapname', value=glc2ocn_ice_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'glc2ocn_ice_rmapname = '// trim(glc2ocn_ice_rmap)
       call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smapname', value=wav2ocn_smap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'wav2ocn_smapname = '// trim(wav2ocn_smap)

       call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_fmapname', value=rof2ocn_fmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'rof2ocn_fmapname = '// trim(rof2ocn_fmap)

       call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_liq_rmapname', value=rof2ocn_liq_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'rof2ocn_liq_rmapname = '// trim(rof2ocn_liq_rmap)
       call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_ice_rmapname', value=rof2ocn_ice_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'rof2ocn_ice_rmapname = '// trim(rof2ocn_ice_rmap)

       ! mapping to rof
       call NUOPC_CompAttributeGet(gcomp, name='lnd2rof_map', value=lnd2rof_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit, '(a)') trim(subname)//'lnd2rof_map = '// trim(lnd2rof_map)

       ! mapping to wav
       call NUOPC_CompAttributeGet(gcomp, name='atm2wav_map', value=atm2wav_map, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit,'(a)') trim(subname)//'atm2wav_map = '// trim(atm2wav_map)

       call NUOPC_CompAttributeGet(gcomp, name='ice2wav_smapname', value=ice2wav_smap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit,'(a)') trim(subname)//'ice2wav_smapname = '// trim(ice2wav_smap)
       call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', value=ocn2wav_smap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit,'(a)') trim(subname)//'ocn2wav_smapname = '// trim(ocn2wav_smap)

       ! uv cart3d mapping
       call NUOPC_CompAttributeGet(gcomp, name='mapuv_with_cart3d', value=cvalue,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) write(logunit,'(a)') trim(subname)//'mapuv_with_cart3d = '// trim(cvalue)
       read(cvalue,*) mapuv_with_cart3d

       ! co2 transfer between componetns
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2a
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2b
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2c

       ! are multiple ice categories being sent from the ice and ocn?
       call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_i2o_per_cat

       ! are multiple ocean depths for temperature and salinity sent from the ocn to glc?
       call NUOPC_CompAttributeGet(gcomp, name='ocn2glc_coupling', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) ocn2glc_coupling

       ! write diagnostic output
       if (mastertask) then
          write(logunit,'(a)') trim(subname)//' flds_co2a = '// trim(cvalue)
          write(logunit,'(a)') trim(subname)//' flds_co2b = '// trim(cvalue)
          write(logunit,'(a)') trim(subname)//' flds_co2c = '// trim(cvalue)
          write(logunit,'(a)') trim(subname)//' flds_i2o_per_cat = '// trim(cvalue)
          write(logunit,'(a)') trim(subname)//' ocn2glc_coupling = '// trim(cvalue)
       end if

    end if

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
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_lfrin')
       call addfld(fldListFr(compocn)%flds, 'So_omask')
       call addfld(fldListFr(compice)%flds, 'Si_imask')
       do ns = 1,num_icesheets
          call addfld(fldlistFr(compglc(ns))%flds, 'Sg_area')
       end do
    else
       call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
    end if

    ! ---------------------------------------------------------------------
    ! to med: atm and ocn fields required for atm/ocn flux calculation'
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_u')
       call addfld(fldListFr(compatm)%flds, 'Sa_v')
       call addfld(fldListFr(compatm)%flds, 'Sa_z')
       call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
       call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
       call addfld(fldListFr(compatm)%flds, 'Sa_shum')
       call addfld(fldListFr(compatm)%flds, 'Sa_ptem')
       call addfld(fldListFr(compatm)%flds, 'Sa_dens')
       call addfld(fldListFr(compatm)%flds, 'Sa_shum_wiso')
    else
       if (is_local%wrap%aoflux_grid == 'ogrid') then
          if (mapuv_with_cart3d) then
             call addmap(fldListFr(compatm)%flds, 'Sa_u' , compocn, mappatch_uv3d, 'one', 'unset')
             call addmap(fldListFr(compatm)%flds, 'Sa_v' , compocn, mappatch_uv3d, 'one', 'unset')
          else
             call addmap(fldListFr(compatm)%flds, 'Sa_u' , compocn, mappatch, 'one', 'unset')
             call addmap(fldListFr(compatm)%flds, 'Sa_v' , compocn, mappatch, 'one', 'unset')
          end if
          call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, mapbilnr, 'one', 'unset')
          call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapbilnr, 'one', 'unset')
          call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, mapbilnr, 'one', 'unset')
          call addmap(fldListFr(compatm)%flds, 'Sa_shum', compocn, mapbilnr, 'one', 'unset')
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, mapbilnr, 'one', 'unset')
          call addmap(fldListFr(compatm)%flds, 'Sa_dens', compocn, mapbilnr, 'one', 'unset')
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum_wiso', rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_shum_wiso', compocn, mapbilnr, 'one', 'unset')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-11): budget implemention needs to be done in CMEPS
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Fall_swnet')
       call addfld(fldListFr(compice)%flds, 'Faii_swnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compice, mapconsf, 'one'  , atm2ice_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one'  , atm2ocn_map)
       end if
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_swnet', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Faii_swnet', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if

    !=====================================================================
    ! FIELDS TO LAND
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! from atm:
    ! to lnd: height at the lowest model level from atm
    ! to lnd: surface height from atm
    ! to lnd: zonal wind at the lowest model level from atm
    ! to lnd: meridional wind at the lowest model level from atm
    ! to lnd: Temperature at the lowest model level from atm
    ! to lnd: potential temperature at the lowest model level from atm
    ! to lnd: Pressure at the lowest model level from atm
    ! to lnd: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------

    allocate(flds(9))
    flds = (/'Sa_z        ',&
             'Sa_topo     ',&
             'Sa_u        ',&
             'Sa_v        ',&
             'Sa_tbot     ',&
             'Sa_ptem     ',&
             'Sa_pbot     ',&
             'Sa_shum     ',&
             'Sa_shum_wiso'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), &
                  complnd, mapbilnr, 'one', atm2lnd_map)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: convective and large scale precipitation rate water equivalent from atm
    ! to lnd: convective and large-scale (stable) snow rate from atm
    ! to lnd: downward longwave heat flux from atm
    ! to lnd: downward direct near-infrared incident solar radiation  from atm
    ! to lnd: downward direct visible incident solar radiation        from atm
    ! to lnd: downward diffuse near-infrared incident solar radiation from atm
    ! to lnd: downward Diffuse visible incident solar radiation       from atm
    ! to lnd: black carbon deposition fluxes from atm
    !   - hydrophylic black carbon dry deposition flux
    !   - hydrophobic black carbon dry deposition flux
    !   - hydrophylic black carbon wet deposition flux
    ! to lnd: organic carbon deposition fluxes from atm
    !   - hydrophylic organic carbon dry deposition flux
    !   - hydrophobic organic carbon dry deposition flux
    !   - hydrophylic organic carbon wet deposition flux
    ! to lnd: dust wet deposition flux (sizes 1-4) from atm
    ! to lnd: dust dry deposition flux (sizes 1-4) from atm
    ! to lnd: nitrogen deposition fields from atm
    ! ---------------------------------------------------------------------

    ! TODO (mvertens, 2018-12-13): the nitrogen deposition fluxes here
    ! are not treated the same was as in cesm2.0 release
    ! TODO (mvertens, 2019-03-10): add water isotopes from atm

    allocate(flds(14))
    flds = (/'Faxa_rainc ',&
             'Faxa_rainl ',&
             'Faxa_snowc ',&
             'Faxa_snowl ',&
             'Faxa_lwdn  ',&
             'Faxa_swndr ',&
             'Faxa_swvdr ',&
             'Faxa_swndf ',&
             'Faxa_swvdf ',&
             'Faxa_bcph  ',&
             'Faxa_ocph  ',&
             'Faxa_dstwet',&
             'Faxa_dstdry',&
             'Faxa_ndep  ' /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if (fldchk(is_local%wrap%FBexp(complnd)         , trim(fldname), rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), &
                  complnd, mapconsf, 'one', atm2lnd_map)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: river channel total water volume from rof
    ! to lnd: river channel main channel water volume from rof
    ! to lnd: river water flux back to land due to flooding
    ! ---------------------------------------------------------------------
    allocate(flds(6))
    flds = (/'Flrr_volr        ',&
             'Flrr_volr_wiso   ',&
             'Flrr_volrmch     ',&
             'Flrr_volrmch_wiso',&
             'Flrr_flood       ',&
             'Flrr_flood_wiso  '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), trim(fldname), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, trim(fldname), &
                  complnd, mapconsf, 'one', rof2lnd_map)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from=comprof, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet grid coverage on global grid from glc
    ! to lnd: ice sheet mask where we are potentially sending non-zero fluxes from glc
    ! to lnd: fields with multiple elevation classes from glc
    ! ---------------------------------------------------------------------

    ! The suffix _elev on land fields designates fields with elevation classes
    ! fields from glc->med do not have elevation classes whereas
    ! fields from med->lnd are in multiple elevation classes

    if (phase == 'advertise') then
       do ns = 1, num_icesheets
          call addfld(fldListFr(compglc(ns))%flds, 'Sg_icemask')     ! ice sheet grid coverage
          call addfld(fldListFr(compglc(ns))%flds, 'Sg_icemask_coupled_fluxes')
          call addfld(fldListFr(compglc(ns))%flds, 'Sg_ice_covered') ! fraction of glacier area
          call addfld(fldListFr(compglc(ns))%flds, 'Sg_topo')        ! surface height of glacer
          call addfld(fldListFr(compglc(ns))%flds, 'Flgg_hflx')      ! downward heat flux from glacier interior
       end do
       call addfld(fldListTo(complnd)%flds, 'Sg_icemask')
       call addfld(fldListTo(complnd)%flds, 'Sg_icemask_coupled_fluxes')
       call addfld(fldListTo(complnd)%flds, 'Sg_ice_covered_elev')
       call addfld(fldListTo(complnd)%flds, 'Sg_topo_elev')
       call addfld(fldListTo(complnd)%flds, 'Flgg_hflx_elev')
    else
       ! custom merge in med_phases_prep_lnd for Sg_icemask and Sg_icemask_coupled_fluxes
       ! custom map merge in med_phases_prep_lnd for Sg_ice_covered_elev, Sg_topo_elev and Flgg_hflx_elev
       if ( fldchk(is_local%wrap%FBExp(complnd), 'Sg_icemask', rc=rc)) then
          do ns = 1, num_icesheets
             if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Sg_icemask', rc=rc)) then
                call addmap(fldListFr(compglc(ns))%flds, 'Sg_icemask', &
                     complnd,  mapconsd, 'one', 'unset')
             end if
          end do
       end if
       if ( fldchk(is_local%wrap%FBExp(complnd), 'Sg_icemask_coupled_fluxes', rc=rc)) then
          do ns = 1, num_icesheets
             if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Sg_icemask_coupled_fluxes', rc=rc)) then
                call addmap(fldListFr(compglc(ns))%flds, 'Sg_icemask_coupled_fluxes', &
                     complnd, mapconsd, 'one', 'unset')
             end if
          end do
       end if
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    !----------------------------------------------------------
    ! to atm: Fractions
    !----------------------------------------------------------
    if (phase == 'advertise') then
       ! the following are computed in med_phases_prep_atm
       call addfld(fldListTo(compatm)%flds, 'Sl_lfrac')
       call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
       call addfld(fldListTo(compatm)%flds, 'So_ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged direct  albedo (visible radiation)
    ! to atm: merged diffuse albedo (visible radiation)
    ! to atm: merged direct  albedo (near-infrared radiation)
    ! to atm: merged diffuse albedo (near-infrared radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_avsdr')
       call addfld(fldListFr(compice)%flds, 'Si_avsdr')
       call addfld(fldListMed_ocnalb%flds , 'So_avsdr')
       call addfld(fldListTo(compatm)%flds, 'Sx_avsdr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdr', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_avsdr', rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_avsdr', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdr', &
                  mrg_from=complnd, mrg_fld='Sl_avsdr', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_avsdr', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_avsdr', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdr', &
                  mrg_from=compice, mrg_fld='Si_avsdr', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_avsdr', rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_avsdr', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdr', &
                  mrg_from=compmed, mrg_fld='So_avsdr', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_avsdf')
       call addfld(fldListFr(compice)%flds, 'Si_avsdf')
       call addfld(fldListMed_ocnalb%flds , 'So_avsdf')
       call addfld(fldListTo(compatm)%flds, 'Sx_avsdf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdf', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_avsdf', rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_avsdf', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdf', &
                  mrg_from=complnd, mrg_fld='Sl_avsdf', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_avsdf', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_avsdf', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdf', &
                  mrg_from=compice, mrg_fld='Si_avsdf', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_avsdf', rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_avsdf', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_avsdf', &
                  mrg_from=compmed, mrg_fld='So_avsdf', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_anidr')
       call addfld(fldListFr(compice)%flds, 'Si_anidr')
       call addfld(fldListMed_ocnalb%flds , 'So_anidr')
       call addfld(fldListTo(compatm)%flds, 'Sx_anidr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidr', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_anidr', rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_anidr', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidr', &
                  mrg_from=complnd, mrg_fld='Sl_anidr', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_anidr', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_anidr', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidr', &
                  mrg_from=compice, mrg_fld='Si_anidr', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_anidr', rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_anidr', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidr', &
                  mrg_from=compmed, mrg_fld='So_anidr', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_anidf')
       call addfld(fldListFr(compice)%flds, 'Si_anidf')
       call addfld(fldListMed_ocnalb%flds , 'So_anidf')
       call addfld(fldListTo(compatm)%flds, 'Sx_anidf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidf', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_anidf', rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_anidf', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidf', &
                  mrg_from=complnd, mrg_fld='Sl_anidf', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_anidf', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_anidf', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidf', &
                  mrg_from=compice, mrg_fld='Si_anidf', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_anidf', rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_anidf', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_anidf', &
                  mrg_from=compmed, mrg_fld='So_anidf', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged reference temperature at 2 meters
    ! to atm: merged 10m wind speed
    ! to atm: merged reference specific humidity at 2 meters
    ! to atm: merged reference specific water isoptope humidity at 2 meters
    ! ---------------------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_tref')
       call addfld(fldListFr(compice)%flds , 'Si_tref')
       call addfld(fldListMed_aoflux%flds  , 'So_tref')
       call addfld(fldListTo(compatm)%flds , 'Sx_tref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_tref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_tref', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_tref', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=complnd, mrg_fld='Sl_tref', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_tref', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_tref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=compice, mrg_fld='Si_tref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_tref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'So_tref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=compmed, mrg_fld='So_tref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_u10')
       call addfld(fldListFr(compice)%flds , 'Si_u10')
       call addfld(fldListMed_aoflux%flds  , 'So_u10')
       call addfld(fldListTo(compatm)%flds , 'Sx_u10')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_u10', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_u10', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_u10', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=complnd, mrg_fld='Sl_u10', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_u10', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_u10', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=compice, mrg_fld='Si_u10', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'So_u10', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=compmed, mrg_fld='So_u10', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_qref')
       call addfld(fldListFr(compice)%flds , 'Si_qref')
       call addfld(fldListMed_aoflux%flds  , 'So_qref')
       call addfld(fldListTo(compatm)%flds , 'Sx_qref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=complnd, mrg_fld='Sl_qref', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_qref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=compice, mrg_fld='Si_qref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'So_qref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=compmed, mrg_fld='So_qref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_qref_wiso')
       call addfld(fldListFr(compice)%flds , 'Si_qref_wiso')
       call addfld(fldListMed_aoflux%flds  , 'So_qref_wiso')
       call addfld(fldListTo(compatm)%flds , 'Sx_qref_wiso')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref_wiso', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref_wiso', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref_wiso', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=complnd, mrg_fld='Sl_qref_wiso', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref_wiso', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_qref_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=compice, mrg_fld='Si_qref_wiso', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref_wiso', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'So_qref_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map) ! map ocn->atm
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=compmed, mrg_fld='So_qref_wiso', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged reference temperature at 2 meters
    ! to atm: merged 10m wind speed
    ! to atm: merged reference specific humidity at 2 meters
    ! to atm: merged reference specific water isoptope humidity at 2 meters
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_tref')
       call addfld(fldListFr(compice)%flds , 'Si_tref')
       call addfld(fldListMed_aoflux%flds  , 'So_tref')
       call addfld(fldListTo(compatm)%flds , 'Sx_tref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_tref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_tref', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_tref', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=complnd, mrg_fld='Sl_tref', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_tref', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_tref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=compice, mrg_fld='Si_tref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_tref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'So_tref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
                  mrg_from=compmed, mrg_fld='So_tref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_u10')
       call addfld(fldListFr(compice)%flds , 'Si_u10')
       call addfld(fldListMed_aoflux%flds  , 'So_u10')
       call addfld(fldListTo(compatm)%flds , 'Sx_u10')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_u10', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_u10', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_u10', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=complnd, mrg_fld='Sl_u10', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_u10', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_u10', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=compice, mrg_fld='Si_u10', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'So_u10', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
                  mrg_from=compmed, mrg_fld='So_u10', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_qref')
       call addfld(fldListFr(compice)%flds , 'Si_qref')
       call addfld(fldListMed_aoflux%flds  , 'So_qref')
       call addfld(fldListTo(compatm)%flds , 'Sx_qref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=complnd, mrg_fld='Sl_qref', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_qref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=compice, mrg_fld='Si_qref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'So_qref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref', &
                  mrg_from=compmed, mrg_fld='So_qref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_qref_wiso')
       call addfld(fldListFr(compice)%flds , 'Si_qref_wiso')
       call addfld(fldListMed_aoflux%flds  , 'So_qref_wiso')
       call addfld(fldListTo(compatm)%flds , 'Sx_qref_wiso')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref_wiso', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref_wiso', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref_wiso', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=complnd, mrg_fld='Sl_qref_wiso', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref_wiso', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Si_qref_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=compice, mrg_fld='Si_qref_wiso', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref_wiso', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'So_qref_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref_wiso', &
                  mrg_from=compmed, mrg_fld='So_qref_wiso', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged zonal surface stress
    ! to atm: merged meridional surface stress
    ! to atm: merged surface latent heat flux
    ! to atm: merged surface sensible heat flux
    ! to atm: merged surface upward longwave heat flux
    ! to atm: evaporation water flux from water
    ! to atm: evaporation water flux from water isotopes
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_taux')
       call addfld(fldListFr(complnd)%flds, 'Fall_taux')
       call addfld(fldListFr(compice)%flds, 'Faii_taux')
       call addfld(fldListMed_aoflux%flds , 'Faox_taux')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_taux', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_taux', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_taux', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_taux', &
                  mrg_from=complnd, mrg_fld='Fall_taux', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_taux', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_taux', &
                  mrg_from=compice, mrg_fld='Faii_taux', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'Faox_taux', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_taux', &
                  mrg_from=compmed, mrg_fld='Faox_taux', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_tauy')
       call addfld(fldListFr(complnd)%flds, 'Fall_tauy')
       call addfld(fldListFr(compice)%flds, 'Faii_tauy')
       call addfld(fldListMed_aoflux%flds , 'Faox_tauy')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_tauy', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_tauy', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_tauy', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_tauy', &
                  mrg_from=complnd, mrg_fld='Fall_tauy', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_tauy', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_tauy', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_tauy', &
                  mrg_from=compice, mrg_fld='Faii_tauy', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'Faox_tauy', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_tauy', &
                  mrg_from=compmed, mrg_fld='Faox_tauy', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_lat')
       call addfld(fldListFr(complnd)%flds, 'Fall_lat')
       call addfld(fldListFr(compice)%flds, 'Faii_lat')
       call addfld(fldListMed_aoflux%flds , 'Faox_lat')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_lat', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_lat', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_lat', &
                  mrg_from=complnd, mrg_fld='Fall_lat', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lat', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_lat', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_lat', &
                  mrg_from=compice, mrg_fld='Faii_lat', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'Faox_lat', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_lat', &
                  mrg_from=compmed, mrg_fld='Faox_lat', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_sen')
       call addfld(fldListFr(complnd)%flds, 'Fall_sen')
       call addfld(fldListFr(compice)%flds, 'Faii_sen')
       call addfld(fldListMed_aoflux%flds , 'Faox_sen')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_sen', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_sen', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_sen', &
                  mrg_from=complnd, mrg_fld='Fall_sen', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_sen', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_sen', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_sen', &
                  mrg_from=compice, mrg_fld='Faii_sen', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'Faox_sen', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_sen', &
                  mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_evap')
       call addfld(fldListFr(complnd)%flds, 'Fall_evap')
       call addfld(fldListFr(compice)%flds, 'Faii_evap')
       call addfld(fldListMed_aoflux%flds , 'Faox_evap')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_evap', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_evap', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap', &
                  mrg_from=complnd, mrg_fld='Fall_evap', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_evap', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap', &
                  mrg_from=compice, mrg_fld='Faii_evap', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds  , 'Faox_evap', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap', &
                  mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_lwup')
       call addfld(fldListFr(complnd)%flds, 'Fall_lwup')
       call addfld(fldListFr(compice)%flds, 'Faii_lwup')
       call addfld(fldListMed_aoflux%flds , 'Faox_lwup')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_lwup', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_lwup', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_lwup', &
                  mrg_from=complnd, mrg_fld='Fall_lwup', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_lwup', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_lwup', &
                  mrg_from=compice, mrg_fld='Faii_lwup', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'Faox_lwup', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds, 'Faxx_lwup', &
                  mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Faxx_evap_wiso')
       call addfld(fldListFr(complnd)%flds, 'Fall_evap_wiso')
       call addfld(fldListFr(compice)%flds, 'Faii_evap_wiso')
       call addfld(fldListMed_aoflux%flds , 'Faox_evap_wiso')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_evap_wiso', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap_wiso', rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Fall_evap_wiso', compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap_wiso', &
                  mrg_from=complnd, mrg_fld='Fall_evap_wiso', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap_wiso', rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_evap_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap_wiso', &
                  mrg_from=compice, mrg_fld='Faii_evap_wiso', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap_wiso', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap(fldListMed_aoflux%flds, 'Faox_evap_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap_wiso', &
                  mrg_from=compmed, mrg_fld='Faox_evap_wiso', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged surface temperature and unmerged temperatures from ice and ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_t')
       call addfld(fldListFr(compice)%flds, 'Si_t')
       call addfld(fldListFr(compocn)%flds, 'So_t')
       call addfld(fldListTo(compatm)%flds, 'So_t')
       call addfld(fldListTo(compatm)%flds, 'Sx_t')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Sx_t', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_t', compatm, mapconsf , 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
                  mrg_from=complnd, mrg_fld='Sl_t', mrg_type='merge', mrg_fracname='lfrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
                  mrg_from=compice, mrg_fld='Si_t', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
             call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
                  mrg_from=compocn, mrg_fld='So_t', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
       if (fldchk(is_local%wrap%FBexp(compatm), 'So_t', rc=rc)) then
          call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface snow depth             from ice (needed for cam)
    ! to atm: mean ice volume per unit area  from ice
    ! to atm: mean snow volume per unit area from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_snowh')
       call addfld(fldListTo(compatm)%flds, 'Si_snowh')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_snowh', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_snowh', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_snowh', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg(fldListTo(compatm)%flds, 'Si_snowh', mrg_from=compice, mrg_fld='Si_snowh', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_vice')
       call addfld(fldListTo(compatm)%flds, 'Si_vice')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_vice', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_vice', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_vice', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg(fldListTo(compatm)%flds, 'Si_vice', mrg_from=compice, mrg_fld='Si_vice', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_vsno')
       call addfld(fldListTo(compatm)%flds, 'Si_vsno')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_vsno', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_vsno', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_vsno', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg(fldListTo(compatm)%flds, 'Si_vsno', mrg_from=compice, mrg_fld='Si_vsno', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface saturation specific humidity in ocean from med aoflux
    ! to atm: square of exch. coeff (tracers)               from med aoflux
    ! to atm: surface fraction velocity                     from med aoflux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'So_ssq')
       call addfld(fldListTo(compatm)%flds , 'So_ssq')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_ssq', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_ssq', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap(fldListMed_aoflux%flds , 'So_ssq', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg(fldListTo(compatm)%flds   , 'So_ssq', mrg_from=compmed, mrg_fld='So_ssq', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'So_re')
       call addfld(fldListTo(compatm)%flds , 'So_re')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_re', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_re', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap(fldListMed_aoflux%flds , 'So_re', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg(fldListTo(compatm)%flds   , 'So_re', mrg_from=compmed, mrg_fld='So_re', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'So_ustar')
       call addfld(fldListTo(compatm)%flds , 'So_ustar')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_ustar', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_ustar', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap(fldListMed_aoflux%flds , 'So_ustar', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg(fldListTo(compatm)%flds   , 'So_ustar', mrg_from=compmed, mrg_fld='So_ustar', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface fraction velocity     from land
    ! to atm: aerodynamic resistance        from land
    ! to atm: surface snow water equivalent from land
    ! ---------------------------------------------------------------------
    allocate(flds(3))
    flds = (/'Sl_fv   ',&
             'Sl_ram1 ',&
             'Sl_snowh'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), &
                  compatm, mapconsf, 'lfrin', lnd2atm_map)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to atm: dust fluxes from land (4 sizes)
    ! ---------------------------------------------------------------------
    fldname = 'Fall_flxdst'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
          call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'lfrin', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='lfrac')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: MEGAN emissions fluxes from land
    !-----------------------------------------------------------------------------
    fldname = 'Fall_voc'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
          call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', atm2lnd_map)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='merge', mrg_fracname='lfrac')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: fire emissions fluxes from land
    !-----------------------------------------------------------------------------
    ! 'wild fire emission fluxes'
    fldname  = 'Fall_fire'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
          call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='merge', mrg_fracname='lfrac')
       end if
    end if

    ! 'wild fire plume height'
    fldname  = 'Sl_fztop'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
          call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: dry deposition velocities from land
    !-----------------------------------------------------------------------------
    fldname  = 'Sl_ddvel'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
          call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
    else
       call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')
       call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_lwdn ',&
             'Faxa_swndr',&
             'Faxa_swndf',&
             'Faxa_swvdr',&
             'Faxa_swvdf'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: surface upward longwave heat flux from mediator
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup')
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwup')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lwup', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwup', &
               mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merged longwave net heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_lwdn')
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup' )
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwnet')
    else
       ! (mom6) (send longwave net to ocn via auto merge)
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_map)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward shortwave heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swdn')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swdn')
    else
       if (fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_swdn', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(compocn)         , 'Faxa_swdn', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swdn', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swdn', &
               mrg_from=compatm, mrg_fld='Faxa_swdn', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')

       call addfld(fldListFr(compice)%flds, 'Fioi_swpen')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdf')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idf')

       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdf')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idf')
    else
       ! Net shortwave ocean (custom calculation in prep_phases_ocn_mod.F90)

       ! import swpen from ice without bands
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen',  compocn, mapfcopy, 'unset', 'unset')
       end if

       ! import swpen from ice by bands
       if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idf', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdr', compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdf', compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idr', compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idf', compocn, mapfcopy, 'unset', 'unset')
       end if

       ! import sw from atm by bands
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc) .and. &
            (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet'    , rc=rc)) .or. &
            (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc) .and. &
             fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', rc=rc) .and. &
             fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', rc=rc) .and. &
             fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', rc=rc))) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_map)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: per ice thickness fraction and sw penetrating into ocean from ice
    ! ---------------------------------------------------------------------
    if (flds_i2o_per_cat) then
       if (phase == 'advertise') then
          ! 'fractional ice coverage wrt ocean for each thickness category '
          call addfld(fldListFr(compice)%flds, 'Si_ifrac_n')
          call addfld(fldListTo(compocn)%flds, 'Si_ifrac_n')

          ! net shortwave radiation penetrating into ocean for each thickness category
          call addfld(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n')
          call addfld(fldListTo(compocn)%flds, 'Fioi_swpen_ifrac_n')

          ! 'fractional atmosphere coverage wrt ocean' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Sf_afrac')
          ! 'fractional atmosphere coverage used in radiation computations wrt ocean' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Sf_afracr')
          ! 'net shortwave radiation times atmosphere fraction' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_afracr')
       else
          call addmap(fldListFr(compice)%flds, 'Si_ifrac_n', &
               compocn, mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Si_ifrac_n', &
               mrg_from=compice, mrg_fld='Si_ifrac_n', mrg_type='copy')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n', &
               compocn, mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_swpen_ifrac_n', &
               mrg_from=compice, mrg_fld='Fioi_swpen_ifrac_n', mrg_type='copy')
          ! Note that 'Sf_afrac, 'Sf_afracr' and 'Foxx_swnet_afracr' will have explicit merging in med_phases_prep_ocn
       end if
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate water equivalent from atm
    !  to ocn: snow rate water equivalent from atm
    ! ---------------------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainl')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
       call addfld(fldListTo(compocn)%flds, 'Faxa_rain' )

       call addfld(fldListFr(compatm)%flds, 'Faxa_rainc_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainl_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rain_wiso' )
       call addfld(fldListTo(compocn)%flds, 'Faxa_rain_wiso' )

       call addfld(fldListFr(compatm)%flds, 'Faxa_snowc')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snowl')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snow' )
       call addfld(fldListTo(compocn)%flds, 'Faxa_snow' )

       call addfld(fldListFr(compatm)%flds, 'Faxa_snowc_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snowl_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snow_wiso' )
       call addfld(fldListTo(compocn)%flds, 'Faxa_snow_wiso' )
    else
       do n = 1,2
          ! Note that the mediator atm/ocn flux calculation needs Faxa_rainc for the gustiness parameterization
          ! which by default is not actually used
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain' //iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             if (iso(n) == ' ') then
                call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n) , &
                     mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
                     mrg_type='sum_with_weights', mrg_fracname='ofrac')
             else
                call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n) , &
                     mrg_from=compatm, mrg_fld=trim('Faxa_rainc'//iso(n))//':'//trim('Faxa_rainl'//iso(n)), &
                     mrg_type='sum_with_weights', mrg_fracname='ofrac')
             end if
          else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n), &
                  mrg_from=compatm, mrg_fld='Faxa_rain'//iso(n), mrg_type='copy')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             if (iso(n) == ' ') then
                call addmrg(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n) , &
                     mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', &
                     mrg_type='sum_with_weights', mrg_fracname='ofrac')
             else
                call addmrg(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n) , &
                     mrg_from=compatm, mrg_fld=trim('Faxa_snowc'//iso(n))//':'//trim('Faxa_snowl'//iso(n)), &
                     mrg_type='sum_with_weights', mrg_fracname='ofrac')
             end if
          else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n), &
                  mrg_from=compatm, mrg_fld='Faxa_snow'//iso(n), mrg_type='copy')
          end if
       end do
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merged sensible heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_sen')
       call addfld(fldListMed_aoflux%flds  , 'Faox_sen')
       call addfld(fldListFr(compice)%flds , 'Fioi_melth')
       call addfld(fldListTo(compocn)%flds , 'Foxx_sen')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
               mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then

       call addfld(fldListFr(compatm)%flds, 'Faxa_lat' )
       call addfld(fldListMed_aoflux%flds , 'Faox_lat' )
       call addfld(fldListMed_aoflux%flds , 'Faox_evap')
       call addfld(fldListTo(compocn)%flds, 'Foxx_lat' )
       call addfld(fldListTo(compocn)%flds, 'Foxx_evap')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lat', &
               mrg_from=compmed, mrg_fld='Faox_lat', mrg_type='merge', mrg_fracname='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_evap', &
               mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds , 'Faox_lat_wiso' )
       call addfld(fldListTo(compocn)%flds, 'Foxx_lat_wiso' )
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat_wiso', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lat_wiso', &
               mrg_from=compmed, mrg_fld='Faox_lat_wiso', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: wind speed squared at 10 meters from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds , 'So_duu10n')
       call addfld(fldListTo(compocn)%flds, 'So_duu10n')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'So_duu10n', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'So_duu10n', rc=rc)) then

          call addmap(fldListMed_aoflux%flds , 'So_duu10n', compatm, mapconsf, 'ofrac', ocn2atm_map) ! map ocn->atm
          call addmrg(fldListTo(compocn)%flds, 'So_duu10n', &
               mrg_from=compmed, mrg_fld='So_duu10n', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
       call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Sa_pslv', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Sa_pslv', rc=rc)) then

          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compice, mapbilnr, 'one', atm2ocn_map)

          call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', &
               mrg_from=compatm, mrg_fld='Sa_pslv', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: black  carbon deposition fluxes from atm
    !   - hydrophylic black carbon dry deposition flux
    !   - hydrophobic black carbon dry deposition flux
    !   - hydrophylic black carbon wet deposition flux
    ! to ocn: organic carbon deposition fluxes from atm
    !   - hydrophylic organic carbon dry deposition flux
    !   - hydrophobic organic carbon dry deposition flux
    !   - hydrophylic organic carbon wet deposition flux
    ! to ocn: dust wet deposition flux (sizes 1-4) from atm
    ! to ocn: dust dry deposition flux (sizes 1-4) from atm
    ! to ocn: nitrogen deposition fields (2) from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_bcph  ', 'Faxa_ocph  ', 'Faxa_dstwet' , 'Faxa_dstdry', 'Faxa_ndep  ' /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: merge zonal and meridional surface stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_taux')
       call addfld(fldListFr(compice)%flds , 'Fioi_taux')
       call addfld(fldListFr(compatm)%flds , 'Faxa_taux')
       call addfld(fldListTo(compocn)%flds , 'Foxx_taux')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', &
               mrg_from=compice, mrg_fld='Fioi_taux', mrg_type='merge', mrg_fracname='ifrac')
          call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', &
               mrg_from=compmed, mrg_fld='Faox_taux', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_tauy')
       call addfld(fldListFr(compice)%flds , 'Fioi_tauy')
       call addfld(fldListFr(compatm)%flds , 'Faxa_tauy')
       call addfld(fldListTo(compocn)%flds , 'Foxx_tauy')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_tauy', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', &
               mrg_from=compice, mrg_fld='Fioi_tauy', mrg_type='merge', mrg_fracname='ifrac')
          call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', &
               mrg_from=compmed, mrg_fld='Faox_tauy', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: water flux due to melting ice from ice
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds , 'Fioi_meltw'//iso(n))
          call addfld(fldListTo(compocn)%flds , 'Fioi_meltw'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_meltw'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_meltw'//iso(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_meltw'//iso(n),    compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Fioi_meltw'//iso(n), &
                  mrg_from=compice, mrg_fld='Fioi_meltw'//iso(n), mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: heat flux from melting ice from ice
    ! to ocn: salt flux from ice
    ! to ocn: hydrophylic black carbon deposition flux from ice
    ! to ocn: hydrophobic black carbon deposition flux from ice
    ! to ocn: dust flux from ice
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-07): is fioi_melth being handled here?
    ! Is fd.yaml correctly aliasing Fioi_melth?

    allocate(flds(5))
    flds = (/'Fioi_melth ',&
             'Fioi_salt  ',&
             'Fioi_bcphi ',&
             'Fioi_bcpho ',&
             'Fioi_flxdst'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end do
    deallocate(flds)

    !-----------------------------
    ! to ocn: liquid runoff from rof and glc components
    ! to ocn: frozen runoff flux from rof and glc components
    ! to ocn: waterflux back to ocn due to flooding from rof
    !-----------------------------

    if (phase == 'advertise') then
       do n = 1,size(iso)
          ! Note that Flrr_flood below needs to be added to
          ! fldlistFr(comprof) in order to be mapped correctly but the ocean
          ! does not receive it so it is advertised but it will! not be connected
          do ns = 1, num_icesheets
             call addfld(fldListFr(compglc(ns))%flds, 'Fogg_rofl'//iso(n))
          end do
          call addfld(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Flrr_flood'//iso(n))
          do ns = 1, num_icesheets
             call addfld(fldListFr(compglc(ns))%flds, 'Fogg_rofi'//iso(n))
          end do
          call addfld(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n))
       end do
    else
       do n = 1,size(iso)
          if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl'//iso(n) , rc=rc)) then
             ! liquid from river and possibly flood from river to ocean
             if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n) , rc=rc)) then
                if (trim(rof2ocn_liq_rmap) == 'unset') then
                   call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), &
                        compocn, mapconsd, 'none', 'unset')
                else
                   call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), &
                        compocn, map_rof2ocn_liq, 'none', rof2ocn_liq_rmap)
                end if
                if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood'//iso(n), rc=rc)) then
                   call addmap(fldListFr(comprof)%flds, 'Flrr_flood'//iso(n), &
                        compocn, mapconsd, 'one', rof2ocn_fmap)
                   call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                        mrg_from=comprof, mrg_fld='Forr_rofl:Flrr_flood', mrg_type='sum')
                else
                   call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                        mrg_from=comprof, mrg_fld='Forr_rofl', mrg_type='sum')
                end if
             end if
             ! liquid from glc to ocean
             do ns = 1, num_icesheets
                if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fogg_rofl'//iso(n) , rc=rc)) then
                   ! TODO: this custom map needs to be different for every ice sheet - how will this be handled?
                   call addmap(fldListFr(compglc(ns))%flds, 'Fogg_rofl'//iso(n), &
                        compocn, map_glc2ocn_liq, 'one' , glc2ocn_liq_rmap)
                   call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                        mrg_from=compglc(ns), mrg_fld='Fogg_rofl'//iso(n), mrg_type='sum')
                end if
             end do
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi'//iso(n) , rc=rc)) then
             ! ice from river to ocean
             if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n) , rc=rc)) then
                if (trim(rof2ocn_ice_rmap) == 'unset') then
                   call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), &
                        compocn, mapconsd, 'none', 'unset')
                else
                   call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), &
                        compocn, map_rof2ocn_ice, 'none', rof2ocn_ice_rmap)
                end if
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                     mrg_from=comprof, mrg_fld='Forr_rofi', mrg_type='sum')
             end if
             ! ice from glc to ocean
             do ns = 1, num_icesheets
                if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fogg_rofi'//iso(n) , rc=rc)) then
                   ! TODO: this custom map needs to be different for every ice sheet - how will this be handled?
                   call addmap(fldListFr(compglc(ns))%flds, 'Fogg_rofi'//iso(n), &
                        compocn, map_glc2ocn_ice, 'one', glc2ocn_ice_rmap)
                   call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                        mrg_from=compglc(ns), mrg_fld='Fogg_rofi'//iso(n), mrg_type='sum')
                end if
             end do
          end if
       end do
    end if

    !-----------------------------
    ! to ocn: Langmuir multiplier from wave
    ! to ocn: Stokes drift u component from wave
    ! to ocn: Stokes drift v component from wave
    ! to ocn: Stokes drift depth from wave
    !-----------------------------
    allocate(flds(4))
    flds = (/'Sw_lamult ',&
             'Sw_ustokes',&
             'Sw_vstokes',&
             'Sw_hstokes'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compwav)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav, compwav), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compwav)%flds, trim(fldname), &
                  compocn,  mapbilnr, 'one', wav2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compwav, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: downward longwave heat flux from atm
    ! to ice: downward direct near-infrared incident solar radiation  from atm
    ! to ice: downward direct visible incident solar radiation        from atm
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! to ice: downward Diffuse visible incident solar radiation       from atm
    ! to ice: hydrophylic black carbon dry deposition flux   from atm
    ! to ice: hydrophobic black carbon dry deposition flux   from atm
    ! to ice: hydrophylic black carbon wet deposition flux   from atm
    ! to ice: hydrophylic organic carbon dry deposition flux from atm
    ! to ice: hydrophobic organic carbon dry deposition flux from atm
    ! to ice: hydrophylic organic carbon wet deposition flux from atm
    ! to ice: dust wet deposition flux (size 1) from atm
    ! to ice: dust wet deposition flux (size 2) from atm
    ! to ice: dust wet deposition flux (size 3) from atm
    ! to ice: dust wet deposition flux (size 4) from atm
    ! to ice: dust dry deposition flux (size 1) from atm
    ! to ice: dust dry deposition flux (size 2) from atm
    ! to ice: dust dry deposition flux (size 3) from atm
    ! to ice: dust dry deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    allocate(flds(9))
    flds = (/'Faxa_lwdn  '    , 'Faxa_swndr '   , 'Faxa_swvdr '   , 'Faxa_swndf ' , 'Faxa_swvdf ', &
             'Faxa_bcph  '    , 'Faxa_ocph  '    , 'Faxa_dstwet'  , 'Faxa_dstdry' /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapconsf, 'one', atm2ice_map)
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: convective and large scale precipitation rate water equivalent from atm
    ! to ice: rain and snow rate from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainl')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
       call addfld(fldListTo(compice)%flds, 'Faxa_rain' )

       call addfld(fldListFr(compatm)%flds, 'Faxa_rainc_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainl_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rain_wiso' )
       call addfld(fldListTo(compice)%flds, 'Faxa_rain_wiso' )
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainc', compice, mapconsf, 'one', atm2ice_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainl', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain' , &
               mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain', &
               mrg_from=compatm, mrg_fld='Faxa_rain', mrg_type='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain_wiso' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainc_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainl_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain_wiso' , &
               mrg_from=compatm, mrg_fld='Faxa_rainc_wiso:Faxa_rainl_wiso', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain_wiso', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain_wiso', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_rain_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain_wiso', &
               mrg_from=compatm, mrg_fld='Faxa_rain_wiso', mrg_type='copy')
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_snowc')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snowl')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snow' )
       call addfld(fldListTo(compice)%flds, 'Faxa_snow' )

       call addfld(fldListFr(compatm)%flds, 'Faxa_snowc_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snowl_wiso')
       call addfld(fldListFr(compatm)%flds, 'Faxa_snow_wiso' )
       call addfld(fldListTo(compice)%flds, 'Faxa_snow_wiso' )
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_snowc', compice, mapconsf, 'one', atm2ice_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_snowl', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_snow' , &
               mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_snow', &
               mrg_from=compatm, mrg_fld='Faxa_snow', mrg_type='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc_wiso', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_snowc_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmap(fldListFr(compatm)%flds, 'Faxa_snowl_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_snow_wiso' , &
               mrg_from=compatm, mrg_fld='Faxa_snowc_wiso:Faxa_snowl_wiso', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow_wiso', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow_wiso', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_snow_wiso', compice, mapconsf, 'one', atm2ice_map)
          call addmrg(fldListTo(compice)%flds, 'Faxa_snow_wiso', &
               mrg_from=compatm, mrg_fld='Faxa_snow_wiso', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: height at the lowest model level from atm
    ! to ice: pressure at the lowest model level fromatm
    ! to ice: temperature at the lowest model level from atm
    ! to ice: potential temperature at the lowest model level from atm
    ! to ice: density at the lowest model level from atm
    ! to ice: zonal wind at the lowest model level from atm
    ! to ice: meridional wind at the lowest model level from atm
    ! to ice: specific humidity at the lowest model level from atm
    ! to ice: specific humidity for water isotopes at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(9))
    flds = (/'Sa_z        ', 'Sa_pbot     ', 'Sa_tbot     ', 'Sa_ptem     ', &
             'Sa_dens     ', 'Sa_u        ', 'Sa_v        ', 'Sa_shum     ', 'Sa_shum_wiso'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             if (trim(fldname) == 'Sa_u' .or. trim(fldname) == 'Sa_v') then
                if (mapuv_with_cart3d) then
                   call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mappatch_uv3d, 'one', atm2ice_map)
                else
                   call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mappatch, 'one', atm2ice_map)
                end if
             else
                call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapbilnr, 'one', atm2ice_map)
             end if
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: sea surface temperature from ocn
    ! to ice: sea surface salinity from ocn
    ! to ice: zonal sea water velocity from ocn
    ! to ice: meridional sea water velocity from ocn
    ! to ice: zonal sea surface slope from ocean
    ! to ice: meridional sea surface slope from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(6))
    flds = (/'So_t   ', 'So_s   ', 'So_u   ', 'So_v   ', 'So_dhdx', 'So_dhdy'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy , 'unset', 'unset')
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: ocean melt and freeze potential from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'Fioo_q')
       call addfld(fldListTo(compice)%flds, 'Fioo_q')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'Fioo_q', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'Fioo_q', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'Fioo_q', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'Fioo_q', &
               mrg_from=compocn, mrg_fld='Fioo_q', mrg_type='copy')
       end if
    end if

    !-----------------------------
    ! to ice: Ratio of ocean surface level abund. H2_16O/H2O/Rstd from ocean
    !-----------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_roce_wiso')
       call addfld(fldListTo(compice)%flds, 'So_roce_wiso')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_roce_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_roce_wiso', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_roce_wiso', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_roce_wiso', &
               mrg_from=compocn, mrg_fld='So_roce_wiso', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: frozen runoff from rof and glc
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Firr_rofi'//iso(n)) ! water flux into sea ice due to runoff (frozen)
          do ns = 1, num_icesheets
             call addfld(fldListFr(compglc(ns))%flds, 'Figg_rofi'//iso(n)) ! glc frozen runoff_iceberg flux to ice
          end do
          call addfld(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n)) ! total frozen water flux into sea ice
       else
          if ( fldchk(is_local%wrap%FBExp(compice), 'Fixx_rofi'//iso(n), rc=rc)) then
             if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then
                call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), &
                     compice, mapconsf, 'none', rof2ocn_ice_rmap)
                call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                     mrg_from=comprof, mrg_fld='Firr_rofi'//iso(n), mrg_type='sum')
             end if
             do ns = 1, num_icesheets
                if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Figg_rofi'//iso(n), rc=rc)) then
                   call addmap(fldListFr(compglc(ns))%flds, 'Figg_rofi'//iso(n), &
                        compice, mapconsf, 'one' , glc2ice_rmap)
                   call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                        mrg_from=compglc(ns), mrg_fld='Figg_rofi'//iso(n), mrg_type='sum')
                end if
             end do
          end if
       end if
    end do

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    !----------------------------------------------------------
    ! to wav: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compwav)%flds, 'Si_ifrac')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Si_ifrac', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_ifrac', rc=rc)) then
             ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
          call addmap(fldListFr(compice)%flds, 'Si_ifrac', compwav, mapbilnr, 'one', ice2wav_smap)
          call addmrg(fldListTo(compwav)%flds, 'Si_ifrac', &
               mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: ocean boundary layer depth from ocn
    ! to wav: ocean currents from ocn
    ! to wav: ocean surface temperature from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(4))
    flds = (/'So_t      ', 'So_u      ', 'So_v      ', 'So_bldepth'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, trim(fldname))
          call addfld(fldListTo(compwav)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compwav)         , trim(fldname), rc=rc)) then
             ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
             call addmap(fldListFr(compocn)%flds, trim(fldname), compwav, mapbilnr, 'one', ocn2wav_smap)
             call addmrg(fldListTo(compwav)%flds, trim(fldname), &
                  mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to wav: zonal wind at the lowest model level from atm
    ! to wav: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(3))
    flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_tbot'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compwav)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compwav, mapbilnr, 'one', atm2wav_map)
             call addmrg(fldListTo(compwav)%flds, trim(fldname), &
                  mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO RIVER (comprof)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid surface)
    ! to rof: water flux from land (liquid glacier, wetland, and lake)
    ! to rof: water flux from land (liquid subsurface)
    ! to rof: water flux from land direct to ocean
    ! to rof: irrigation flux from land (withdrawal from rivers)
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-13): the following isotopes have not yet been defined in the NUOPC field dict
    ! allocate(flds(12))
    ! flds = (/'Flrl_rofsur', 'Flrl_rofsur_wiso', 'Flrl_rofgwl', 'Flrl_rofgwl_wiso', &
    !          'Flrl_rofsub', 'Flrl_rofsub_wiso', 'Flrl_rofdto', 'Flrl_rofdto_wiso', &
    !          'Flrl_rofi'  , 'Flrl_rofi_wiso'  , 'Flrl_irrig' , 'Flrl_irrig_wiso'   /)

    allocate(flds(6))
    flds = (/'Flrl_rofsur', 'Flrl_rofgwl', 'Flrl_rofsub', 'Flrl_rofdto', 'Flrl_rofi  ', 'Flrl_irrig '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(comprof)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), comprof, mapconsf, 'lfrac', lnd2rof_map)
             call addmrg(fldListTo(comprof)%flds, trim(fldname), &
                  mrg_from=complnd, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='lfrac')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO LAND-ICE (compglc)
    !=====================================================================

    !-----------------------------
    ! to glc: from land
    !-----------------------------
    ! - fields sent from lnd->med ARE in multiple elevation classes
    ! - fields sent from med->glc do NOT have elevation classes

    ! Sets a coupling field for all glc elevation classes (1:glc_nec) plus bare land (index 0).
    ! Note that, if glc_nec = 0, then we don't create any coupling fields (not even the bare land (0) fldindex)
    ! Note : Sl_topo is sent from lnd -> med, but is NOT sent to glc (only used for the remapping in the mediator)

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_tsrf_elev')   ! surface temperature of glacier (1->glc_nec+1)
       call addfld(fldListFr(complnd)%flds, 'Sl_topo_elev')   ! surface heights of glacier     (1->glc_nec+1)
       call addfld(fldListFr(complnd)%flds, 'Flgl_qice_elev') ! glacier ice flux               (1->glc_nec+1)
       do ns = 1,num_icesheets
          call addfld(fldListTo(compglc(ns))%flds, 'Sl_tsrf')
          call addfld(fldListTo(compglc(ns))%flds, 'Flgl_qice')
       end do
    else
       ! custom mapping, accumulation and merging will be done in prep_glc_mod.F90
       do ns = 1,num_icesheets
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Flgl_qice_elev', rc=rc)) then
             call addmap(FldListFr(complnd)%flds, 'Flgl_qice_elev', compglc(ns), mapbilnr, 'lfrac', 'unset')
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Sl_tsrf_elev'  , rc=rc)) then
             call addmap(FldListFr(complnd)%flds, 'Sl_tsrf_elev', compglc(ns), mapbilnr, 'lfrac', 'unset')
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Sl_topo_elev'  , rc=rc)) then
             ! This is needed just for mappingn to glc - but is not sent as a field
             call addmap(FldListFr(complnd)%flds, 'Sl_topo_elev', compglc(ns), mapbilnr, 'lfrac', 'unset')
          end if
       end do
    end if

    !-----------------------------
    ! to glc: from ocn
    !-----------------------------
    if (ocn2glc_coupling) then
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, 'So_t_depth')
          call addfld(fldListFr(compocn)%flds, 'So_s_depth')
          do ns = 1,num_icesheets
             call addfld(fldListTo(compglc(ns))%flds, 'So_t_depth')
             call addfld(fldListTo(compglc(ns))%flds, 'So_s_depth')
          end do
       else
          ! custom mapping, accumulation and merging will be done in prep_glc_mod.F90
          ! the following is used to create the route handle
          do ns = 1,num_icesheets
             if ( fldchk(is_local%wrap%FBImp(compocn,compocn) , 'So_t_depth', rc=rc)) then
                call addmap(FldListFr(compocn)%flds, 'So_t_depth', compglc(ns), mapbilnr, 'none', 'unset')
             end if
             if ( fldchk(is_local%wrap%FBImp(compocn,compocn) , 'So_s_depth', rc=rc)) then
                call addmap(FldListFr(compocn)%flds, 'So_s_depth', compglc(ns), mapbilnr, 'none', 'unset')
             end if
          end do
       end if
    end if

    !=====================================================================
    ! CO2 EXCHANGE
    !=====================================================================

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO)

    if (flds_co2a) then
       ! ---------------------------------------------------------------------
       ! to lnd and ocn: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', compocn, mapbilnr, 'one', atm2ocn_map)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2prog', &
               mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', compocn, mapbilnr, 'one', atm2ocn_map)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2diag', &
               mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
       end if

    else if (flds_co2b) then

       ! ---------------------------------------------------------------------
       ! to lnd: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from land
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Fall_fco2_lnd')
          call addfld(fldListTo(compatm)%flds, 'Fall_fco2_lnd')
       else
          call addmap(fldListFr(complnd)%flds, 'Fall_fco2_lnd', compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
               mrg_from=complnd, mrg_fld='Fall_fco2_lnd', mrg_type='copy_with_weights', mrg_fracname='lfrac')
       end if

    else if (flds_co2c) then

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', compocn, mapbilnr, 'one', atm2ocn_map)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2prog', &
               mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', compocn, mapbilnr, 'one', atm2ocn_map)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2diag', &
               mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from land
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Fall_fco2_lnd')
          call addfld(fldListTo(compatm)%flds, 'Fall_fco2_lnd')
       else
          call addmap(fldListFr(complnd)%flds, 'Fall_fco2_lnd', compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
               mrg_from=complnd, mrg_fld='Fall_fco2_lnd', mrg_type='copy_with_weights', mrg_fracname='lfrac')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from ocn
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, 'Faoo_fco2_ocn')
          call addfld(fldListTo(compatm)%flds, 'Faoo_fco2_ocn')
       else
          call addmap(fldListFr(compocn)%flds, 'Faoo_fco2_ocn', compatm, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    endif

    !-----------------------------------------------------------------------------
    ! CARMA fields (volumetric soil water)
    !-----------------------------------------------------------------------------
    ! TODO: add this

  end subroutine esmFldsExchange_cesm

end module esmFldsExchange_cesm_mod
