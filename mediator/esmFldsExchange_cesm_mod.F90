module esmFldsExchange_cesm_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !
  ! -----------------------------------------------------------------------------------------
  ! subroutine med_fldList_addmrg_to(index, fldname, mrg_from, mrg_fld, mrg_type, mrg_fracname, rc)
  !   integer         , intent(in)            :: index
  !   character(len=*), intent(in)            :: fldname
  !   integer         , intent(in)            :: mrg_from
  !   character(len=*), intent(in)            :: mrg_fld
  !   character(len=*), intent(in)            :: mrg_type
  !   character(len=*), intent(in) , optional :: mrg_fracname
  !   integer         , intent(out), optional :: rc
  !
  ! index        : destination component index that merging will occur to
  ! fldname      : field name in mediator export field bundle for destination component
  ! mrg_from     : source component index that will contribute to the merge
  ! mrg_fld      : field name fom source component field bundle that will be used in merge
  ! mrg_type     : one of ['copy', 'copy_with_weights', 'sum', 'sum_with_weights', 'merge']
  ! mrg_fracname : if mrg_type is copy_with_weights or merge -
  !                fraction name in fraction field bundle to use in merge
  !
  ! -----------------------------------------------------------------------------------------
  ! subroutine med_fldList_addmap_from(index, fldname, destcomp, maptype, mapnorm, mapfile)
  !   integer          , intent(in)           :: index
  !   character(len=*) , intent(in)           :: fldname
  !   integer          , intent(in)           :: destcomp
  !   integer          , intent(in)           :: maptype
  !   character(len=*) , intent(in)           :: mapnorm
  !   character(len=*) , intent(in), optional :: mapfile
  !
  ! index        : source component index that mapping will occur from
  ! fldname      : field name in mediator import field for source component
  ! destcomp     : destination component index
  ! maptype      : mapping type (see med_internal_state_mod.F90 for the supported mapping types)
  !                if maptype is mapfcopy - create a redistribution route handle
  ! mapnorm      : normalization type, one of ['unset', 'one', 'none', fracname]
  !                fracname -  is the field name of the field in the fraction field bundle corresponding to the
  !                source field that will be used for normalization
  !                'one' - implies that the mapped field is divided by mapping 'one' from the source to the
  !                destination mesh
  !                'none'  - do not use any normalization - use if maytype is not mapfcopy
  !                'unset' - do not use any normalization - only used if maptype is mapfcopy
  ! mapfile      : if mapfile is idmap - create a redistribution route nhandle
  !                if mapfile is unset then create the mapping route handle at run time
  !
  ! -----------------------------------------------------------------------------------------
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
  use med_internalstate_mod , only : logunit, maintask, samegrid_atmlnd
  use med_internalstate_mod , only : mrg_fracname_lnd2atm_state, mrg_fracname_lnd2atm_flux, map_fracname_lnd2atm
  use med_internalstate_mod , only : mrg_fracname_lnd2rof, map_fracname_lnd2rof
  use med_internalstate_mod , only : mrg_fracname_lnd2glc, map_fracname_lnd2glc

  implicit none
  public

  public :: esmFldsExchange_cesm

  ! currently required mapping files
  character(len=CX)   :: rof2ocn_ice_rmap ='unset'
  character(len=CX)   :: rof2ocn_liq_rmap ='unset'
  character(len=CX)   :: rof2lnd_map = 'unset'
  character(len=CX)   :: lnd2rof_map = 'unset'

  ! no mapping files (value is 'idmap' or 'unset')
  character(len=CX)   :: atm2ice_map = 'unset'
  character(len=CX)   :: atm2ocn_map = 'unset'
  character(len=CX)   :: atm2lnd_map = 'unset'
  character(len=CX)   :: atm2wav_map = 'unset'
  character(len=CX)   :: ice2atm_map = 'unset'
  character(len=CX)   :: ice2wav_map = 'unset'
  character(len=CX)   :: lnd2atm_map = 'unset'
  character(len=CX)   :: ocn2atm_map = 'unset'
  character(len=CX)   :: ocn2wav_map = 'unset'
  character(len=CX)   :: rof2ocn_map = 'unset'
  character(len=CX)   :: wav2ocn_map = 'unset'

  logical             :: mapuv_with_cart3d              ! Map U/V vector wind fields from ATM to OCN/ICE by rotating in Cartesian 3D space and then back
  logical             :: flds_i2o_per_cat               ! Ice thickness category fields passed to OCN
  logical             :: flds_co2a                      ! Pass CO2 from ATM to surface components
  logical             :: flds_co2b                      ! Pass CO2 from ATM to LND and back from LND to ATM
  logical             :: flds_co2c                      ! Pass CO2 from ATM to surface (OCN/LND) and back from them to ATM
  logical             :: flds_wiso                      ! Pass water isotop fields
  logical             :: flds_r2l_stream_channel_depths ! Pass channel depths from ROF to LND

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
    use med_internalstate_mod , only : InternalState, logunit, maintask
    use med_internalstate_mod , only : compmed, compatm, complnd, compocn
    use med_internalstate_mod , only : compice, comprof, compwav, compglc, ncomps
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mappatch_uv3d, mapbilnr_nstod
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use med_internalstate_mod , only : map_rof2ocn_ice, map_rof2ocn_liq
    use esmFlds               , only : addfld_ocnalb => med_fldList_addfld_ocnalb
    use esmFlds               , only : addfld_aoflux => med_fldList_addfld_aoflux
    use esmFlds               , only : addmap_aoflux => med_fldList_addmap_aoflux
    use esmFlds               , only : addmap_ocnalb => med_fldList_addmap_ocnalb
    use esmFlds               , only : addfld_to => med_fldList_addfld_to
    use esmFlds               , only : addfld_from => med_fldList_addfld_from
    use esmFlds               , only : addmap_from => med_fldList_addmap_from
    use esmFlds               , only : addmrg_to => med_fldList_addmrg_to

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: n, ns
    character(len=CL)   :: atm_mesh_name
    character(len=CL)   :: lnd_mesh_name
    character(len=CL)   :: ice_mesh_name
    character(len=CL)   :: ocn_mesh_name
    character(len=CL)   :: cvalue
    character(len=CS)   :: mrgfld_source
    logical             :: wav_coupling_to_cice
    logical             :: ocn2glc_coupling
    character(len=*) , parameter   :: subname=' (esmFldsExchange_cesm) '
    !--------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='wav_coupling_to_cice', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wav_coupling_to_cice

    call NUOPC_CompAttributeGet(gcomp, name='ocn2glc_coupling', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn2glc_coupling

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (phase == 'advertise') then

      ! determine if atm and lnd have the same mesh
      call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=atm_mesh_name, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call NUOPC_CompAttributeGet(gcomp, name='mesh_lnd', value=lnd_mesh_name, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=ice_mesh_name, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=ocn_mesh_name, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (trim(atm_mesh_name) == trim(lnd_mesh_name)) then
        atm2lnd_map = 'idmap'
        lnd2atm_map = 'idmap'
      end if
      if (trim(atm_mesh_name) == trim(ocn_mesh_name)) then
        atm2ocn_map = 'idmap'
        ocn2atm_map = 'idmap'
      end if
      if (trim(atm_mesh_name) == trim(ice_mesh_name)) then
        atm2ice_map = 'idmap'
        ice2atm_map = 'idmap'
      end if

       ! mapping rof=>lnd and lnd=>rof - the following two maps are needed for MIZUROUTE
       call NUOPC_CompAttributeGet(gcomp, name='rof2lnd_map', value=rof2lnd_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (maintask) write(logunit, '(a)') trim(subname)//'rof2lnd_map = '// trim(rof2lnd_map)
       call NUOPC_CompAttributeGet(gcomp, name='lnd2rof_map', value=lnd2rof_map,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (maintask) write(logunit, '(a)') trim(subname)//'lnd2rof_map = '// trim(lnd2rof_map)

       ! mapping to rof => ocn with custom mapping
       call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_liq_rmapname', value=rof2ocn_liq_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (maintask) write(logunit, '(a)') trim(subname)//'rof2ocn_liq_rmapname = '// trim(rof2ocn_liq_rmap)
       call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_ice_rmapname', value=rof2ocn_ice_rmap,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (maintask) write(logunit, '(a)') trim(subname)//'rof2ocn_ice_rmapname = '// trim(rof2ocn_ice_rmap)

       ! uv cart3d mapping
       call NUOPC_CompAttributeGet(gcomp, name='mapuv_with_cart3d', value=cvalue,  rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) mapuv_with_cart3d

       ! is co2 transfer between components enabled?
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

       ! are water isotope exchanges enabled?
       call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_wiso

       call NUOPC_CompAttributeGet(gcomp, name='flds_r2l_stream_channel_depths', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_r2l_stream_channel_depths

       ! write diagnostic output
       if (maintask) then
          write(logunit,'(a)'   ) ' flds_co2a: prognostic and diagnostic CO2 at lowest atm level is sent to lnd and ocn'
          write(logunit,'(a)'   ) ' flds_co2b: prognostic and diagnostic CO2 at lowest atm level is sent to lnd and ocn'
          write(logunit,'(a)'   ) '            and surface flux of CO2 from lnd is sent back to atm'
          write(logunit,'(a)'   ) ' flds_co2c: prognostic and diagnostic CO2 at lowest atm level is sent to lnd and ocn'
          write(logunit,'(a)'   ) '            and surface flux of CO2 from lnd is sent back to atm'
          write(logunit,'(a)'   ) '            and surface flux of CO2 from ocn is sent back to atm'
          write(logunit,'(a,l7)') trim(subname)//' flds_co2a                       = ',flds_co2a
          write(logunit,'(a,l7)') trim(subname)//' flds_co2b                       = ',flds_co2b
          write(logunit,'(a,l7)') trim(subname)//' flds_co2c                       = ',flds_co2c
          write(logunit,'(a,l7)') trim(subname)//' flds_wiso                       = ',flds_wiso
          write(logunit,'(a,l7)') trim(subname)//' flds_i2o_per_cat                = ',flds_i2o_per_cat
          write(logunit,'(a,l7)') trim(subname)//' flds_r2l_stream_channel_depths  = ',flds_r2l_stream_channel_depths
          write(logunit,'(a,l7)') trim(subname)//' mapuv_with_cart3d               = ',mapuv_with_cart3d
       end if

    end if

    !=====================================================================
    ! scalar information
    !=====================================================================

    if (phase == 'advertise') then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld_from(n, trim(cvalue))
          call addfld_to(n, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_lfrin')
       call addfld_from(compocn, 'So_omask')
       call addfld_from(compice, 'Si_imask')
       do ns = 1,is_local%wrap%num_icesheets
          call addfld_from(compglc(ns), 'Sg_area')
       end do
    else
       call addmap_from(compocn, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
    end if

    ! ---------------------------------------------------------------------
    ! to med: atm and ocn fields required for atm/ocn flux calculation'
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_u')
       call addfld_from(compatm, 'Sa_v')
       call addfld_from(compatm, 'Sa_z')
       call addfld_from(compatm, 'Sa_tbot')
       call addfld_from(compatm, 'Sa_pbot')
       call addfld_from(compatm, 'Sa_shum')
       call addfld_from(compatm, 'Sa_ptem')
       call addfld_from(compatm, 'Sa_dens')
       call addfld_from(compatm, 'Faxa_rainc')
       if (flds_wiso) then
          call addfld_from(compatm, 'Sa_shum_wiso')
       end if
    else
       if (is_local%wrap%aoflux_grid == 'ogrid') then
          if (mapuv_with_cart3d) then
             call addmap_from(compatm, 'Sa_u' , compocn, mappatch_uv3d, 'one', atm2ocn_map)
             call addmap_from(compatm, 'Sa_v' , compocn, mappatch_uv3d, 'one', atm2ocn_map)
          else
             call addmap_from(compatm, 'Sa_u' , compocn, mappatch, 'one', atm2ocn_map)
             call addmap_from(compatm, 'Sa_v' , compocn, mappatch, 'one', atm2ocn_map)
          end if
          call addmap_from(compatm, 'Faxa_rainc', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_z'   , compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_tbot', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_pbot', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_shum', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_ptem', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_dens', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_pslv', compocn, mapbilnr, 'one', atm2ocn_map)
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum_wiso', rc=rc)) then
             call addmap_from(compatm, 'Sa_shum_wiso', compocn, mapbilnr, 'one', atm2ocn_map)
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Fall_swnet')
       call addfld_from(compice, 'Faii_swnet')
       call addfld_from(compatm, 'Faxa_swnet')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swnet', compice, mapconsf, 'one'  , atm2ice_map)
          call addmap_from(compatm, 'Faxa_swnet', compocn, mapconsf, 'one'  , atm2ocn_map)
       end if
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_swnet', rc=rc)) then
          call addmap_from(compice, 'Faii_swnet', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if

    !=====================================================================
    ! FIELDS TO LAND
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to lnd: height at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_z')
       call addfld_to(complnd, 'Sa_z')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_z', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_z', rc=rc)) then
          call addmap_from(compatm, 'Sa_z', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_z', mrg_from=compatm, mrg_fld='Sa_z', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: surface height from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_topo')
       call addfld_to(complnd, 'Sa_topo')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_topo', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_topo', rc=rc)) then
          call addmap_from(compatm, 'Sa_topo', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_topo', mrg_from=compatm, mrg_fld='Sa_topo', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: zonal wind at the lowest model level from atm
    ! to lnd: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_u')
       call addfld_to(complnd, 'Sa_u')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u', rc=rc)) then
          call addmap_from(compatm, 'Sa_u', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_u', mrg_from=compatm, mrg_fld='Sa_u', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_v')
       call addfld_to(complnd, 'Sa_v')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v', rc=rc)) then
          call addmap_from(compatm, 'Sa_v', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_v', mrg_from=compatm, mrg_fld='Sa_v', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: pressure at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_pbot')
       call addfld_to(complnd, 'Sa_pbot')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_pbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_pbot', rc=rc)) then
          call addmap_from(compatm, 'Sa_pbot', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_pbot', mrg_from=compatm, mrg_fld='Sa_pbot', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: o3 at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_o3')
       call addfld_to(complnd, 'Sa_o3')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_o3', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_o3', rc=rc)) then
          call addmap_from(compatm, 'Sa_o3', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_o3', mrg_from=compatm, mrg_fld='Sa_o3', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: cld to grnd lightning flash freq
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_lightning')
       call addfld_to(complnd, 'Sa_lightning')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_lightning', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_lightning', rc=rc)) then
          call addmap_from(compatm, 'Sa_lightning', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_lightning', mrg_from=compatm, mrg_fld='Sa_lightning', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: temperature at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_tbot')
       call addfld_to(complnd, 'Sa_tbot')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_tbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_tbot', rc=rc)) then
          call addmap_from(compatm, 'Sa_tbot', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_tbot', mrg_from=compatm, mrg_fld='Sa_tbot', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: potential temperature at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_ptem')
       call addfld_to(complnd, 'Sa_ptem')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_ptem', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_ptem', rc=rc)) then
          call addmap_from(compatm, 'Sa_ptem', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_ptem', mrg_from=compatm, mrg_fld='Sa_ptem', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_shum')
       call addfld_to(complnd, 'Sa_shum')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_shum', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_shum', rc=rc)) then
          call addmap_from(compatm, 'Sa_shum', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_shum', mrg_from=compatm, mrg_fld='Sa_shum', mrg_type='copy')
       end if
    end if
    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compatm, 'Sa_shum_wiso')
          call addfld_to(complnd, 'Sa_shum_wiso')
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_shum_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_shum_wiso', rc=rc)) then
             call addmap_from(compatm, 'Sa_shum_wiso', complnd, mapbilnr, 'one', atm2lnd_map)
             call addmrg_to(complnd, 'Sa_shum_wiso', mrg_from=compatm, mrg_fld='Sa_shum_wiso', mrg_type='copy')
          end if
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: prognostic CO2 at the lowest atm model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_co2prog')
       call addfld_to(complnd, 'Sa_co2prog')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_co2prog', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_co2prog', rc=rc)) then
          call addmap_from(compatm, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_co2prog', mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: diagnostic CO2 at the lowest atm model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_co2diag')
       call addfld_to(complnd, 'Sa_co2diag')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_co2diag', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_co2diag', rc=rc)) then
          call addmap_from(compatm, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Sa_co2diag', mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: convective and large scale precipitation rate water equivalent from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_rainc')
       call addfld_to(complnd, 'Faxa_rainc')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_rainc', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_rainc', rc=rc)) then
          call addmap_from(compatm, 'Faxa_rainc', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_rainc', mrg_from=compatm, mrg_fld='Faxa_rainc', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_rainl')
       call addfld_to(complnd, 'Faxa_rainl')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_rainl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_rainl', rc=rc)) then
          call addmap_from(compatm, 'Faxa_rainl', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_rainl', mrg_from=compatm, mrg_fld='Faxa_rainl', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: convective and large-scale (stable) snow rate from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_snowc')
       call addfld_to(complnd, 'Faxa_snowc')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_snowc', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_snowc', rc=rc)) then
          call addmap_from(compatm, 'Faxa_snowc', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_snowc', mrg_from=compatm, mrg_fld='Faxa_snowc', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_snowl')
       call addfld_to(complnd, 'Faxa_snowl')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_snowl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_snowl', rc=rc)) then
          call addmap_from(compatm, 'Faxa_snowl', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_snowl', mrg_from=compatm, mrg_fld='Faxa_snowl', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: downward longwave heat flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_lwdn')
       call addfld_to(complnd, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_lwdn', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_lwdn', rc=rc)) then
          call addmap_from(compatm, 'Faxa_lwdn', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_lwdn', mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: downward direct near-infrared incident solar radiation  from atm
    ! to lnd: downward direct visible incident solar radiation        from atm
    ! to lnd: downward diffuse near-infrared incident solar radiation from atm
    ! to lnd: downward Diffuse visible incident solar radiation       from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndr')
       call addfld_to(complnd, 'Faxa_swndr')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_swndr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndr', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_swndr', mrg_from=compatm, mrg_fld='Faxa_swndr', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdr')
       call addfld_to(complnd, 'Faxa_swvdr')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_swvdr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdr', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_swvdr', mrg_from=compatm, mrg_fld='Faxa_swvdr', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndf')
       call addfld_to(complnd, 'Faxa_swndf')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_swndf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_swndf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndf', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_swndf', mrg_from=compatm, mrg_fld='Faxa_swndf', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdf')
       call addfld_to(complnd, 'Faxa_swvdf')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_swvdf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdf', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_swvdf', mrg_from=compatm, mrg_fld='Faxa_swvdf', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_bcph')
       call addfld_to(complnd, 'Faxa_bcph')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_bcph', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_bcph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_bcph', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_bcph', mrg_from=compatm, mrg_fld='Faxa_bcph', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: black carbon deposition fluxes from atm
    !   - hydrophylic black carbon dry deposition flux
    !   - hydrophobic black carbon dry deposition flux
    !   - hydrophylic black carbon wet deposition flux
    ! to lnd: organic carbon deposition fluxes from atm
    !   - hydrophylic organic carbon dry deposition flux
    !   - hydrophobic organic carbon dry deposition flux
    !   - hydrophylic organic carbon wet deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_ocph')
       call addfld_to(complnd, 'Faxa_ocph')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_ocph', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_ocph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_ocph', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_ocph', mrg_from=compatm, mrg_fld='Faxa_ocph', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: dust wet deposition flux (sizes 1-4) from atm
    ! to lnd: dust dry deposition flux (sizes 1-4) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_dstwet')
       call addfld_to(complnd, 'Faxa_dstwet')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_dstwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_dstwet', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstwet', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_dstwet', mrg_from=compatm, mrg_fld='Faxa_dstwet', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_dstdry')
       call addfld_to(complnd, 'Faxa_dstdry')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_dstdry', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_dstdry', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstdry', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_dstdry', mrg_from=compatm, mrg_fld='Faxa_dstdry', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to lnd: nitrogen deposition fields from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_ndep')
       call addfld_to(complnd, 'Faxa_ndep')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Faxa_ndep', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Faxa_ndep', rc=rc)) then
          call addmap_from(compatm, 'Faxa_ndep', complnd, mapconsf, 'one', atm2lnd_map)
          call addmrg_to(complnd, 'Faxa_ndep', mrg_from=compatm, mrg_fld='Faxa_ndep', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd: river channel total water volume from rof
    ! to lnd: river channel main channel water volume from rof
    ! to lnd: river water flux back to land due to flooding
    ! to lnd: tributary water depth
    ! to lnd: tributary channel depth
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(comprof, 'Flrr_volr')
       call addfld_to(complnd, 'Flrr_volr')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_volr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volr', rc=rc)) then
          call addmap_from(comprof, 'Flrr_volr', complnd, mapconsf, 'one', rof2lnd_map)
          call addmrg_to(complnd, 'Flrr_volr', mrg_from=comprof, mrg_fld='Flrr_volr', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(comprof, 'Flrr_volrmch')
       call addfld_to(complnd, 'Flrr_volrmch')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_volrmch', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volrmch', rc=rc)) then
          call addmap_from(comprof, 'Flrr_volrmch', complnd, mapconsf, 'one', rof2lnd_map)
          call addmrg_to(complnd, 'Flrr_volrmch', mrg_from=comprof, mrg_fld='Flrr_volrmch', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(comprof, 'Flrr_flood')
       call addfld_to(complnd, 'Flrr_flood')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_flood', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood', rc=rc)) then
          call addmap_from(comprof, 'Flrr_flood', complnd, mapconsf, 'one', rof2lnd_map)
          call addmrg_to(complnd, 'Flrr_flood', mrg_from=comprof, mrg_fld='Flrr_flood', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(comprof, 'Sr_tdepth')
       call addfld_to(complnd, 'Sr_tdepth')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Sr_tdepth', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(comprof, comprof), 'Sr_tdepth', rc=rc)) then
          call addmap_from(comprof, 'Sr_tdepth', complnd, mapconsf, 'one', rof2lnd_map)
          call addmrg_to(complnd, 'Sr_tdepth', mrg_from=comprof, mrg_fld='Sr_tdepth', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(comprof, 'Sr_tdepth_max')
       call addfld_to(complnd, 'Sr_tdepth_max')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Sr_tdepth_max', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(comprof, comprof), 'Sr_tdepth_max', rc=rc)) then
          call addmap_from(comprof, 'Sr_tdepth_max', complnd, mapconsf, 'one', rof2lnd_map)
          call addmrg_to(complnd, 'Sr_tdepth_max', mrg_from=comprof, mrg_fld='Sr_tdepth_max', mrg_type='copy')
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(comprof, 'Flrr_volr_wiso')
          call addfld_to(complnd, 'Flrr_volr_wiso')
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_volr_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volr_wiso', rc=rc)) then
             call addmap_from(comprof, 'Flrr_volr_wiso', complnd, mapconsf, 'one', rof2lnd_map)
             call addmrg_to(complnd, 'Flrr_volr_wiso', &
                  mrg_from=comprof, mrg_fld='Flrr_volr_wiso', mrg_type='copy')
          end if
       end if
       if (phase == 'advertise') then
          call addfld_from(comprof, 'Flrr_volrmch_wiso')
          call addfld_to(complnd, 'Flrr_volrmch_wiso')
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_volrmch_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volrmch_wiso', rc=rc)) then
             call addmap_from(comprof, 'Flrr_volrmch_wiso', complnd, mapconsf, 'one', rof2lnd_map)
             call addmrg_to(complnd, 'Flrr_volrmch_wiso', &
                  mrg_from=comprof, mrg_fld='Flrr_volrmch_wiso', mrg_type='copy')
          end if
       end if
       if (phase == 'advertise') then
          call addfld_from(comprof, 'Flrr_flood_wiso')
          call addfld_to(complnd, 'Flrr_flood_wiso')
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flrr_flood_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood_wiso', rc=rc)) then
             call addmap_from(comprof, 'Flrr_flood_wiso', complnd, mapconsf, 'one', rof2lnd_map)
             call addmrg_to(complnd, 'Flrr_flood_wiso', &
                  mrg_from=comprof, mrg_fld='Flrr_flood_wiso', mrg_type='copy')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet grid coverage on global grid from glc
    ! to lnd: ice sheet mask where we are potentially sending non-zero fluxes from glc
    ! to lnd: fields with multiple elevation classes from glc
    ! ---------------------------------------------------------------------

    ! The suffix _elev on land fields designates fields with elevation classes
    ! fields from glc->med do not have elevation classes whereas
    ! fields from med->lnd are in multiple elevation classes

    if (phase == 'advertise') then
       do ns = 1, is_local%wrap%num_icesheets
          call addfld_from(compglc(ns), 'Sg_icemask')     ! ice sheet grid coverage
          call addfld_from(compglc(ns), 'Sg_icemask_coupled_fluxes')
          call addfld_from(compglc(ns), 'Sg_ice_covered') ! fraction of glacier area
          call addfld_from(compglc(ns), 'Sg_topo')        ! surface height of glacer
          call addfld_from(compglc(ns), 'Flgg_hflx')      ! downward heat flux from glacier interior
       end do
       call addfld_to(complnd, 'Sg_icemask')
       call addfld_to(complnd, 'Sg_icemask_coupled_fluxes')
       call addfld_to(complnd, 'Sg_ice_covered_elev')
       call addfld_to(complnd, 'Sg_topo_elev')
       call addfld_to(complnd, 'Flgg_hflx_elev')
    else
       ! custom merge in med_phases_prep_lnd for Sg_icemask and Sg_icemask_coupled_fluxes
       ! custom map merge in med_phases_prep_lnd for Sg_ice_covered_elev, Sg_topo_elev and Flgg_hflx_elev
       if ( fldchk(is_local%wrap%FBExp(complnd), 'Sg_icemask', rc=rc)) then
          do ns = 1, is_local%wrap%num_icesheets
             if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Sg_icemask', rc=rc)) then
                call addmap_from(compglc(ns), 'Sg_icemask', &
                     complnd,  mapconsd, 'one', 'unset')
             end if
          end do
       end if
       if ( fldchk(is_local%wrap%FBExp(complnd), 'Sg_icemask_coupled_fluxes', rc=rc)) then
          do ns = 1, is_local%wrap%num_icesheets
             if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Sg_icemask_coupled_fluxes', rc=rc)) then
                call addmap_from(compglc(ns), 'Sg_icemask_coupled_fluxes', &
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
       call addfld_to(compatm, 'Sl_lfrac')
       call addfld_to(compatm, 'Si_ifrac')
       call addfld_to(compatm, 'So_ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged direct  albedo (visible radiation)
    ! to atm: merged diffuse albedo (visible radiation)
    ! to atm: merged direct  albedo (near-infrared radiation)
    ! to atm: merged diffuse albedo (near-infrared radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_avsdr')
       call addfld_from(compice, 'Si_avsdr')
       call addfld_ocnalb('So_avsdr')
       call addfld_to(compatm, 'Sx_avsdr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdr', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_avsdr', rc=rc)) then
             call addmap_from(complnd, 'Sl_avsdr', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm, 'Sx_avsdr', &
                  mrg_from=complnd, mrg_fld='Sl_avsdr', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_avsdr', rc=rc)) then
             call addmap_from(compice, 'Si_avsdr', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm, 'Sx_avsdr', &
                  mrg_from=compice, mrg_fld='Si_avsdr', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_avsdr', rc=rc)) then
             call addmap_ocnalb( 'So_avsdr', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg_to(compatm, 'Sx_avsdr', &
                  mrg_from=compmed, mrg_fld='So_avsdr', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_avsdf')
       call addfld_from(compice, 'Si_avsdf')
       call addfld_ocnalb( 'So_avsdf')
       call addfld_to(compatm, 'Sx_avsdf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdf', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_avsdf', rc=rc)) then
             call addmap_from(complnd, 'Sl_avsdf', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm, 'Sx_avsdf', &
                  mrg_from=complnd, mrg_fld='Sl_avsdf', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_avsdf', rc=rc)) then
             call addmap_from(compice, 'Si_avsdf', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm, 'Sx_avsdf', &
                  mrg_from=compice, mrg_fld='Si_avsdf', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_avsdf', rc=rc)) then
             call addmap_ocnalb( 'So_avsdf', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg_to(compatm, 'Sx_avsdf', &
                  mrg_from=compmed, mrg_fld='So_avsdf', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_anidr')
       call addfld_from(compice, 'Si_anidr')
       call addfld_ocnalb( 'So_anidr')
       call addfld_to(compatm, 'Sx_anidr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidr', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_anidr', rc=rc)) then
             call addmap_from(complnd, 'Sl_anidr', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm, 'Sx_anidr', &
                  mrg_from=complnd, mrg_fld='Sl_anidr', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_anidr', rc=rc)) then
             call addmap_from(compice, 'Si_anidr', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm, 'Sx_anidr', &
                  mrg_from=compice, mrg_fld='Si_anidr', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_anidr', rc=rc)) then
             call addmap_ocnalb( 'So_anidr', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg_to(compatm, 'Sx_anidr', &
                  mrg_from=compmed, mrg_fld='So_anidr', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_anidf')
       call addfld_from(compice, 'Si_anidf')
       call addfld_ocnalb( 'So_anidf')
       call addfld_to(compatm, 'Sx_anidf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidf', rc=rc)) then
          ! Note that for aqua-plant there will be no import from complnd or compice - and the
          ! current logic below takes care of this.
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_anidf', rc=rc)) then
             call addmap_from(complnd, 'Sl_anidf', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm, 'Sx_anidf', &
                  mrg_from=complnd, mrg_fld='Sl_anidf', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_anidf', rc=rc)) then
             call addmap_from(compice, 'Si_anidf', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm, 'Sx_anidf', &
                  mrg_from=compice, mrg_fld='Si_anidf', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_anidf', rc=rc)) then
             call addmap_ocnalb( 'So_anidf', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg_to(compatm, 'Sx_anidf', &
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
       call addfld_from(complnd , 'Sl_tref')
       call addfld_from(compice , 'Si_tref')
       call addfld_aoflux('So_tref')
       call addfld_to(compatm , 'Sx_tref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_tref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_tref', rc=rc)) then
             call addmap_from(complnd , 'Sl_tref', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=complnd, mrg_fld='Sl_tref', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_tref', rc=rc)) then
             call addmap_from(compice , 'Si_tref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=compice, mrg_fld='Si_tref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_tref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_tref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=compmed, mrg_fld='So_tref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd , 'Sl_u10')
       call addfld_from(compice , 'Si_u10')
       call addfld_aoflux('So_u10')
       call addfld_to(compatm , 'Sx_u10')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_u10', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_u10', rc=rc)) then
             call addmap_from(complnd , 'Sl_u10', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=complnd, mrg_fld='Sl_u10', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_u10', rc=rc)) then
             call addmap_from(compice , 'Si_u10', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=compice, mrg_fld='Si_u10', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_u10', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=compmed, mrg_fld='So_u10', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd , 'Sl_qref')
       call addfld_from(compice , 'Si_qref')
       call addfld_aoflux('So_qref')
       call addfld_to(compatm , 'Sx_qref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref', rc=rc)) then
             call addmap_from(complnd , 'Sl_qref', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=complnd, mrg_fld='Sl_qref', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref', rc=rc)) then
             call addmap_from(compice , 'Si_qref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=compice, mrg_fld='Si_qref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_qref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=compmed, mrg_fld='So_qref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(complnd , 'Sl_qref_wiso')
          call addfld_from(compice , 'Si_qref_wiso')
          call addfld_aoflux('So_qref_wiso')
          call addfld_to(compatm , 'Sx_qref_wiso')
       else
          if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref_wiso', rc=rc)) then
             if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref_wiso', rc=rc)) then
                call addmap_from(complnd , 'Sl_qref_wiso', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=complnd, mrg_fld='Sl_qref_wiso', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
             end if
             if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref_wiso', rc=rc)) then
                call addmap_from(compice , 'Si_qref_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=compice, mrg_fld='Si_qref_wiso', mrg_type='merge', mrg_fracname='ifrac')
             end if
             if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref_wiso', rc=rc)) then
                if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                   call addmap_aoflux( 'So_qref_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map) ! map ocn->atm
                end if
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=compmed, mrg_fld='So_qref_wiso', mrg_type='merge', mrg_fracname='ofrac')
             end if
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
       call addfld_from(complnd , 'Sl_tref')
       call addfld_from(compice , 'Si_tref')
       call addfld_aoflux('So_tref')
       call addfld_to(compatm , 'Sx_tref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_tref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_tref', rc=rc)) then
             call addmap_from(complnd , 'Sl_tref', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=complnd, mrg_fld='Sl_tref', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_tref', rc=rc)) then
             call addmap_from(compice , 'Si_tref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=compice, mrg_fld='Si_tref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_tref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_tref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_tref', &
                  mrg_from=compmed, mrg_fld='So_tref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd , 'Sl_u10')
       call addfld_from(compice , 'Si_u10')
       call addfld_aoflux('So_u10')
       call addfld_to(compatm , 'Sx_u10')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_u10', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_u10', rc=rc)) then
             call addmap_from(complnd , 'Sl_u10', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=complnd, mrg_fld='Sl_u10', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_u10', rc=rc)) then
             call addmap_from(compice , 'Si_u10', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=compice, mrg_fld='Si_u10', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_u10', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_u10', &
                  mrg_from=compmed, mrg_fld='So_u10', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_from(complnd , 'Sl_qref')
       call addfld_from(compice , 'Si_qref')
       call addfld_aoflux('So_qref')
       call addfld_to(compatm , 'Sx_qref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref', rc=rc)) then
             call addmap_from(complnd , 'Sl_qref', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=complnd, mrg_fld='Sl_qref', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref', rc=rc)) then
             call addmap_from(compice , 'Si_qref', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=compice, mrg_fld='Si_qref', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_qref', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Sx_qref', &
                  mrg_from=compmed, mrg_fld='So_qref', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(complnd , 'Sl_qref_wiso')
          call addfld_from(compice , 'Si_qref_wiso')
          call addfld_aoflux('So_qref_wiso')
          call addfld_to(compatm , 'Sx_qref_wiso')
       else
          if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_qref_wiso', rc=rc)) then
             if (fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref_wiso', rc=rc)) then
                call addmap_from(complnd , 'Sl_qref_wiso', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=complnd, mrg_fld='Sl_qref_wiso', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
             end if
             if (fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref_wiso', rc=rc)) then
                call addmap_from(compice , 'Si_qref_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=compice, mrg_fld='Si_qref_wiso', mrg_type='merge', mrg_fracname='ifrac')
             end if
             if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref_wiso', rc=rc)) then
                if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                   call addmap_aoflux('So_qref_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map)
                end if
                call addmrg_to(compatm , 'Sx_qref_wiso', &
                     mrg_from=compmed, mrg_fld='So_qref_wiso', mrg_type='merge', mrg_fracname='ofrac')
             end if
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
       call addfld_to(compatm, 'Faxx_taux')
       call addfld_from(complnd, 'Fall_taux')
       call addfld_from(compice, 'Faii_taux')
       call addfld_aoflux( 'Faox_taux')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_taux', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_taux', rc=rc)) then
             call addmap_from(complnd , 'Fall_taux', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_taux', &
                  mrg_from=complnd, mrg_fld='Fall_taux', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc)) then
             call addmap_from(compice , 'Faii_taux', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_taux', &
                  mrg_from=compice, mrg_fld='Faii_taux', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_taux', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Faxx_taux', &
                  mrg_from=compmed, mrg_fld='Faox_taux', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_to(compatm, 'Faxx_tauy')
       call addfld_from(complnd, 'Fall_tauy')
       call addfld_from(compice, 'Faii_tauy')
       call addfld_aoflux( 'Faox_tauy')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_tauy', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_tauy', rc=rc)) then
             call addmap_from(complnd , 'Fall_tauy', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_tauy', &
                  mrg_from=complnd, mrg_fld='Fall_tauy', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_tauy', rc=rc)) then
             call addmap_from(compice , 'Faii_tauy', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_tauy', &
                  mrg_from=compice, mrg_fld='Faii_tauy', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_tauy', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Faxx_tauy', &
                  mrg_from=compmed, mrg_fld='Faox_tauy', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_to(compatm, 'Faxx_lat')
       call addfld_from(complnd, 'Fall_lat')
       call addfld_from(compice, 'Faii_lat')
       call addfld_aoflux( 'Faox_lat')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_lat', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', rc=rc)) then
             call addmap_from(complnd , 'Fall_lat', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_lat', &
                  mrg_from=complnd, mrg_fld='Fall_lat', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lat', rc=rc)) then
             call addmap_from(compice , 'Faii_lat', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_lat', &
                  mrg_from=compice, mrg_fld='Faii_lat', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_lat', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Faxx_lat', &
                  mrg_from=compmed, mrg_fld='Faox_lat', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_to(compatm, 'Faxx_sen')
       call addfld_from(complnd, 'Fall_sen')
       call addfld_from(compice, 'Faii_sen')
       call addfld_aoflux( 'Faox_sen')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_sen', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', rc=rc)) then
             call addmap_from(complnd , 'Fall_sen', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_sen', &
                  mrg_from=complnd, mrg_fld='Fall_sen', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_sen', rc=rc)) then
             call addmap_from(compice , 'Faii_sen', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_sen', &
                  mrg_from=compice, mrg_fld='Faii_sen', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_sen', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Faxx_sen', &
                  mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_to(compatm, 'Faxx_evap')
       call addfld_from(complnd, 'Fall_evap')
       call addfld_from(compice, 'Faii_evap')
       call addfld_aoflux( 'Faox_evap')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_evap', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap', rc=rc)) then
             call addmap_from(complnd , 'Fall_evap', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_evap', &
                  mrg_from=complnd, mrg_fld='Fall_evap', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap', rc=rc)) then
             call addmap_from(compice , 'Faii_evap', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_evap', &
                  mrg_from=compice, mrg_fld='Faii_evap', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_evap', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'Faxx_evap', &
                  mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_to(compatm, 'Faxx_lwup')
       call addfld_from(complnd, 'Fall_lwup')
       call addfld_from(compice, 'Faii_lwup')
       call addfld_aoflux( 'Faox_lwup')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_lwup', rc=rc)) then
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', rc=rc)) then
             call addmap_from(complnd , 'Fall_lwup', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm , 'Faxx_lwup', &
                  mrg_from=complnd, mrg_fld='Fall_lwup', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', rc=rc)) then
             call addmap_from(compice , 'Faii_lwup', compatm, mapconsf, 'ifrac', ice2atm_map)
             call addmrg_to(compatm , 'Faxx_lwup', &
                  mrg_from=compice, mrg_fld='Faii_lwup', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('Faox_lwup', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm, 'Faxx_lwup', &
                  mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_to(compatm, 'Faxx_evap_wiso')
          call addfld_from(complnd, 'Fall_evap_wiso')
          call addfld_from(compice, 'Faii_evap_wiso')
          call addfld_aoflux( 'Faox_evap_wiso')
       else
          if (fldchk(is_local%wrap%FBexp(compatm), 'Faxx_evap_wiso', rc=rc)) then
             if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap_wiso', rc=rc)) then
                call addmap_from(complnd , 'Fall_evap_wiso', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
                call addmrg_to(compatm , 'Faxx_evap_wiso', &
                     mrg_from=complnd, mrg_fld='Fall_evap_wiso', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_flux)
             end if
             if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap_wiso', rc=rc)) then
                call addmap_from(compice , 'Faii_evap_wiso', compatm, mapconsf, 'ifrac', ice2atm_map)
                call addmrg_to(compatm , 'Faxx_evap_wiso', &
                     mrg_from=compice, mrg_fld='Faii_evap_wiso', mrg_type='merge', mrg_fracname='ifrac')
             end if
             if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap_wiso', rc=rc)) then
                if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                   call addmap_aoflux('Faox_evap_wiso', compatm, mapconsf, 'ofrac', ocn2atm_map)
                end if
                call addmrg_to(compatm , 'Faxx_evap_wiso', &
                     mrg_from=compmed, mrg_fld='Faox_evap_wiso', mrg_type='merge', mrg_fracname='ofrac')
             end if
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged surface temperature and unmerged temperatures from ice and ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_t')
       call addfld_from(compice, 'Si_t')
       call addfld_from(compocn, 'So_t')
       call addfld_to(compatm, 'So_t')
       call addfld_to(compatm, 'Sx_t')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Sx_t', rc=rc)) then
          if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc)) then
             call addmap_from(complnd, 'Sl_t', compatm, mapconsf , map_fracname_lnd2atm, lnd2atm_map)
             call addmrg_to(compatm, 'Sx_t', &
                  mrg_from=complnd, mrg_fld='Sl_t', mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
          end if
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
             call addmap_from(compice, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_map)
             call addmrg_to(compatm, 'Sx_t', &
                  mrg_from=compice, mrg_fld='Si_t', mrg_type='merge', mrg_fracname='ifrac')
          end if
          if (fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
             call addmap_from(compocn, 'So_t', compatm, mapconsf, 'ofrac', ocn2atm_map)
             call addmrg_to(compatm, 'Sx_t', &
                  mrg_from=compocn, mrg_fld='So_t', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
       if (fldchk(is_local%wrap%FBexp(compatm), 'So_t', rc=rc)) then
          call addmrg_to(compatm, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: unmerged ugust_out from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_aoflux('So_ugustOut')
       call addfld_to(compatm, 'So_ugustOut')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'So_ugustOut', rc=rc)) then
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_ugustOut', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_ugustOut', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'So_ugustOut', &
                  mrg_from=compmed, mrg_fld='So_ugustOut', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: 10 m winds including/excluding gust component
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_aoflux('So_u10withGust')
       call addfld_to(compatm, 'So_u10withGust')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'So_u10withGust', rc=rc)) then
          if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10withGust', rc=rc)) then
             if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
                call addmap_aoflux('So_u10withGust', compatm, mapconsf, 'ofrac', ocn2atm_map)
             end if
             call addmrg_to(compatm , 'So_u10withGust', &
                  mrg_from=compmed, mrg_fld='So_u10withGust', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    if (phase == 'advertise') then
       call addfld_aoflux('So_u10res')
       call addfld_to(compatm, 'So_u10res')
    else
      if ( fldchk(is_local%wrap%FBexp(compatm), 'So_u10res', rc=rc)) then
         if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_u10res', rc=rc)) then
            if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
               call addmap_aoflux('So_u10res', compatm, mapconsf, 'ofrac', ocn2atm_map)
            end if
            call addmrg_to(compatm , 'So_u10res', &
                 mrg_from=compmed, mrg_fld='So_u10res', mrg_type='merge', mrg_fracname='ofrac')
         end if
      end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface snow depth             from ice (needed for cam)
    ! to atm: mean ice volume per unit area  from ice
    ! to atm: mean snow volume per unit area from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Si_snowh')
       call addfld_to(compatm, 'Si_snowh')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_snowh', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_snowh', rc=rc)) then
          call addmap_from(compice, 'Si_snowh', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg_to(compatm, 'Si_snowh', mrg_from=compice, mrg_fld='Si_snowh', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compice, 'Si_vice')
       call addfld_to(compatm, 'Si_vice')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_vice', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_vice', rc=rc)) then
          call addmap_from(compice, 'Si_vice', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg_to(compatm, 'Si_vice', mrg_from=compice, mrg_fld='Si_vice', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compice, 'Si_vsno')
       call addfld_to(compatm, 'Si_vsno')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_vsno', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_vsno', rc=rc)) then
          call addmap_from(compice, 'Si_vsno', compatm, mapconsf, 'ifrac', ice2atm_map)
          call addmrg_to(compatm, 'Si_vsno', mrg_from=compice, mrg_fld='Si_vsno', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface saturation specific humidity in ocean from med aoflux
    ! to atm: square of exch. coeff (tracers)               from med aoflux
    ! to atm: surface fraction velocity                     from med aoflux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_aoflux('So_ssq')
       call addfld_to(compatm , 'So_ssq')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_ssq', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_ssq', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap_aoflux( 'So_ssq', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg_to(compatm   , 'So_ssq', mrg_from=compmed, mrg_fld='So_ssq', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_aoflux('So_re')
       call addfld_to(compatm , 'So_re')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_re', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_re', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap_aoflux( 'So_re', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg_to(compatm   , 'So_re', mrg_from=compmed, mrg_fld='So_re', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_aoflux('So_ustar')
       call addfld_to(compatm , 'So_ustar')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm) , 'So_ustar', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o , 'So_ustar', rc=rc)) then
          if (trim(is_local%wrap%aoflux_grid) == 'ogrid') then
             call addmap_aoflux( 'So_ustar', compatm, mapconsf, 'ofrac', ocn2atm_map)
          end if
          call addmrg_to(compatm   , 'So_ustar', mrg_from=compmed, mrg_fld='So_ustar', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface fraction velocity     from land
    ! to atm: aerodynamic resistance        from land
    ! to atm: surface snow water equivalent from land
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_fv')
       call addfld_to(compatm, 'Sl_fv')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_fv', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_fv', rc=rc)) then
          call addmap_from(complnd, 'Sl_fv', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_fv', mrg_from=complnd, mrg_fld='Sl_fv', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_ram1')
       call addfld_to(compatm, 'Sl_ram1')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_ram1', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_ram1', rc=rc)) then
          call addmap_from(complnd, 'Sl_ram1', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_ram1', mrg_from=complnd, mrg_fld='Sl_ram1', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_snowh')
       call addfld_to(compatm, 'Sl_snowh')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_snowh', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_snowh', rc=rc)) then
          call addmap_from(complnd, 'Sl_snowh', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_snowh', mrg_from=complnd, mrg_fld='Sl_snowh', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: CARMA fields (volumetric soil water) from land
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_soilw')
       call addfld_to(compatm, 'Sl_soilw')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_soilw', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_soilw', rc=rc)) then
          call addmap_from(complnd, 'Sl_soilw', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_soilw', mrg_from=complnd, mrg_fld='Sl_soilw', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: dust fluxes from land (4 sizes)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Fall_flxdst')
       call addfld_to(compatm, 'Fall_flxdst')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Fall_flxdst', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Fall_flxdst', rc=rc)) then
          call addmap_from(complnd, 'Fall_flxdst', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Fall_flxdst', &
               mrg_from=complnd, mrg_fld='Fall_flxdst', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2atm_flux)
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: MEGAN emissions fluxes from land
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Fall_voc')
       call addfld_to(compatm, 'Fall_voc')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Fall_voc', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Fall_voc', rc=rc)) then
          call addmap_from(complnd, 'Fall_voc', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Fall_voc', &
               mrg_from=complnd, mrg_fld='Fall_voc', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2atm_flux)
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: fire emissions fluxes from land
    !-----------------------------------------------------------------------------
    ! 'wild fire emission fluxes'
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Fall_fire')
       call addfld_to(compatm, 'Fall_fire')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Fall_fire', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Fall_fire', rc=rc)) then
          call addmap_from(complnd, 'Fall_fire', compatm, mapconsf, 'lfrin', lnd2atm_map)
          call addmrg_to(compatm, 'Fall_fire', &
               mrg_from=complnd, mrg_fld='Fall_fire', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2atm_flux)
       end if
    end if
    ! 'wild fire plume height'
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_fztop')
       call addfld_to(compatm, 'Sl_fztop')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Sl_fztop', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Sl_fztop', rc=rc)) then
          call addmap_from(complnd, 'Sl_fztop', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_fztop', &
               mrg_from=complnd, mrg_fld='Sl_fztop', mrg_type='copy')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: dry deposition velocities from land
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Sl_ddvel')
       call addfld_to(compatm, 'Sl_ddvel')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Sl_ddvel', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Sl_ddvel', rc=rc)) then
          call addmap_from(complnd, 'Sl_ddvel', compatm, mapconsf, map_fracname_lnd2atm, lnd2atm_map)
          call addmrg_to(compatm, 'Sl_ddvel', mrg_from=complnd, mrg_fld='Sl_ddvel', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface flux of CO2 from land
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Fall_fco2_lnd')
       call addfld_to(compatm, 'Fall_fco2_lnd')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_co2_lnd', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_co2_lnd', rc=rc)) then
          call addmap_from(complnd, 'Fall_fco2_lnd', compatm, mapconsf, 'one', lnd2atm_map)
          call addmrg_to(compatm, 'Fall_fco2_lnd', &
               mrg_from=complnd, mrg_fld='Fall_fco2_lnd', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2atm_flux)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface flux of CO2 from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Faoo_fco2_ocn')
       call addfld_to(compatm, 'Faoo_fco2_ocn')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fco2_ocn', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_fco2_ocn', rc=rc)) then
          call addmap_from(compocn, 'Faoo_fco2_ocn', compatm, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: surface flux of dms from ocean
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Faoo_fdms_ocn')
       call addfld_to(compatm, 'Faoo_fdms_ocn')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fdms_ocn', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_fdms_ocn', rc=rc)) then
          call addmap_from(compocn, 'Faoo_fdms_ocn', compocn, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: surface flux of bromoform from ocean
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Faoo_fbrf_ocn')
       call addfld_to(compatm, 'Faoo_fbrf_ocn')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fbrf_ocn', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_fbrf_ocn', rc=rc)) then
          call addmap_from(compocn, 'Faoo_fbrf_ocn', compocn, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: surface flux of n2o from ocean
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Faoo_fn2o_ocn')
       call addfld_to(compatm, 'Faoo_fn2o_ocn')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fn2o_ocn', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_fn2o_ocn', rc=rc)) then
          call addmap_from(compocn, 'Faoo_fn2o_ocn', compocn, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: surface flux of nh3 from ocean
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Faoo_fnh3_ocn')
       call addfld_to(compatm, 'Faoo_fnh3_ocn')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn,compocn), 'Faoo_fnh3_ocn', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faoo_fnh3_ocn', rc=rc)) then
          call addmap_from(compocn, 'Faoo_fnh3_ocn', compocn, mapconsd, 'one', ocn2atm_map)
          ! custom merge in med_phases_prep_atm
       end if
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Si_ifrac')
       call addfld_to(compocn, 'Si_ifrac')
    else
       call addmap_from(compice, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')
       call addmrg_to(compocn, 'Si_ifrac', mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_lwdn')
       call addfld_to(compocn, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_lwdn', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', rc=rc)) then
          call addmap_from(compatm, 'Faxa_lwdn', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_lwdn', &
               mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndr')
       call addfld_to(compocn, 'Faxa_swndr')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_swndr', &
               mrg_from=compatm, mrg_fld='Faxa_swndr', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndf')
       call addfld_to(compocn, 'Faxa_swndf')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swndf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_swndf', &
               mrg_from=compatm, mrg_fld='Faxa_swndf', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdr')
       call addfld_to(compocn, 'Faxa_swvdr')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_swvdr', &
               mrg_from=compatm, mrg_fld='Faxa_swvdr', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdf')
       call addfld_to(compocn, 'Faxa_swvdf')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_swvdf', &
               mrg_from=compatm, mrg_fld='Faxa_swvdf', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface upward longwave heat flux from mediator
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_aoflux('Faox_lwup')
       call addfld_to(compocn , 'Foxx_lwup')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lwup', rc=rc)) then
          call addmrg_to(compocn, 'Foxx_lwup', &
               mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: merged longwave net heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm , 'Faxa_lwdn')
       call addfld_aoflux('Faox_lwup' )
       call addfld_to(compocn , 'Foxx_lwnet')
    else
       ! (mom6) (send longwave net to ocn via auto merge)
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc)) then
          call addmap_from(compatm, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_map)
          call addmrg_to(compocn, 'Foxx_lwnet', &
               mrg_from=compmed, mrg_fld='Faox_lwup', mrg_type='merge', mrg_fracname='ofrac')
          call addmrg_to(compocn, 'Foxx_lwnet', &
               mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: downward shortwave heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swdn')
       call addfld_to(compocn, 'Faxa_swdn')
    else
       if (fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_swdn', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(compocn)         , 'Faxa_swdn', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swdn', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_swdn', &
               mrg_from=compatm, mrg_fld='Faxa_swdn', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdr')
       call addfld_from(compatm, 'Faxa_swndr')
       call addfld_from(compatm, 'Faxa_swvdf')
       call addfld_from(compatm, 'Faxa_swndf')

       call addfld_from(compice, 'Fioi_swpen')
       call addfld_from(compice, 'Fioi_swpen_vdr')
       call addfld_from(compice, 'Fioi_swpen_vdf')
       call addfld_from(compice, 'Fioi_swpen_idr')
       call addfld_from(compice, 'Fioi_swpen_idf')

       call addfld_to(compocn, 'Foxx_swnet')
       call addfld_to(compocn, 'Foxx_swnet_vdr')
       call addfld_to(compocn, 'Foxx_swnet_vdf')
       call addfld_to(compocn, 'Foxx_swnet_idr')
       call addfld_to(compocn, 'Foxx_swnet_idf')
    else
       ! Net shortwave ocean (custom calculation in prep_phases_ocn_mod.F90)

       ! import swpen from ice without bands
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen', rc=rc)) then
          call addmap_from(compice, 'Fioi_swpen',  compocn, mapfcopy, 'unset', 'unset')
       end if

       ! import swpen from ice by bands
       if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idf', rc=rc)) then
          call addmap_from(compice, 'Fioi_swpen_vdr', compocn, mapfcopy, 'unset', 'unset')
          call addmap_from(compice, 'Fioi_swpen_vdf', compocn, mapfcopy, 'unset', 'unset')
          call addmap_from(compice, 'Fioi_swpen_idr', compocn, mapfcopy, 'unset', 'unset')
          call addmap_from(compice, 'Fioi_swpen_idf', compocn, mapfcopy, 'unset', 'unset')
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
          call addmap_from(compatm, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_map)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: per ice thickness fraction and sw penetrating into ocean from ice
    ! ---------------------------------------------------------------------
    if (flds_i2o_per_cat) then
       if (phase == 'advertise') then
          ! 'fractional ice coverage wrt ocean for each thickness category '
          call addfld_from(compice, 'Si_ifrac_n')
          call addfld_to(compocn, 'Si_ifrac_n')

          ! net shortwave radiation penetrating into ocean for each thickness category
          call addfld_from(compice, 'Fioi_swpen_ifrac_n')
          call addfld_to(compocn, 'Fioi_swpen_ifrac_n')

          ! 'fractional atmosphere coverage wrt ocean' (computed in med_phases_prep_ocn)
          call addfld_to(compocn, 'Sf_afrac')
          ! 'fractional atmosphere coverage used in radiation computations wrt ocean' (computed in med_phases_prep_ocn)
          call addfld_to(compocn, 'Sf_afracr')
          ! 'net shortwave radiation times atmosphere fraction' (computed in med_phases_prep_ocn)
          call addfld_to(compocn, 'Foxx_swnet_afracr')
       else
          call addmap_from(compice, 'Si_ifrac_n', &
               compocn, mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Si_ifrac_n', &
               mrg_from=compice, mrg_fld='Si_ifrac_n', mrg_type='copy')
          call addmap_from(compice, 'Fioi_swpen_ifrac_n', &
               compocn, mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_swpen_ifrac_n', &
               mrg_from=compice, mrg_fld='Fioi_swpen_ifrac_n', mrg_type='copy')
          ! Note that 'Sf_afrac, 'Sf_afracr' and 'Foxx_swnet_afracr' will have explicit merging in med_phases_prep_ocn
       end if
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate water equivalent from atm
    !  to ocn: snow rate water equivalent from atm
    ! ---------------------------------------------------------------------

    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_rainc')
       call addfld_from(compatm, 'Faxa_rainl')
       call addfld_to(compocn, 'Faxa_rain' )
       call addfld_from(compatm, 'Faxa_snowc')
       call addfld_from(compatm, 'Faxa_snowl')
       call addfld_to(compocn, 'Faxa_snow' )
    else
       ! TODO: why are we not merging Faxa_rain and Faxa_snow if they are sent from atm with ofrac
       ! Note that the mediator atm/ocn flux calculation needs Faxa_rainc for the gustiness parameterization
       ! which by default is not actually used
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain' , rc=rc)) then
          call addmap_from(compatm, 'Faxa_rainl', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Faxa_rainc', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
               mrg_type='sum_with_weights', mrg_fracname='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', rc=rc)) then
          call addmap_from(compatm, 'Faxa_snowl', compocn, mapconsf, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Faxa_snowc', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_snow'  , &
               mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', mrg_type='sum_with_weights', mrg_fracname='ofrac')
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compatm, 'Faxa_rainc_wiso')
          call addfld_from(compatm, 'Faxa_rainl_wiso')
          call addfld_to(compocn, 'Faxa_rain_wiso' )
          call addfld_from(compatm, 'Faxa_snowc_wiso')
          call addfld_from(compatm, 'Faxa_snowl_wiso')
          call addfld_from(compatm, 'Faxa_snow_wiso' )
       else
          ! Note that the mediator atm/ocn flux calculation needs Faxa_rainc for the gustiness parameterization
          ! which by default is not actually used
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain_wiso' , rc=rc)) then
             call addmap_from(compatm, 'Faxa_rainl_wiso', compocn, mapconsf, 'one', atm2ocn_map)
             call addmap_from(compatm, 'Faxa_rainc_wiso', compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg_to(compocn, 'Faxa_rain_wiso' , &
                  mrg_from=compatm, mrg_fld=trim('Faxa_rainc_wiso')//':'//trim('Faxa_rainl_wiso'), &
                  mrg_type='sum_with_weights', mrg_fracname='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc_wiso', rc=rc)) then
             call addmap_from(compatm, 'Faxa_snowl_wiso', compocn, mapconsf, 'one', atm2ocn_map)
             call addmap_from(compatm, 'Faxa_snowc_wiso', compocn, mapconsf, 'one', atm2ocn_map)
             call addmrg_to(compocn, 'Faxa_snow_wiso', &
                  mrg_from=compatm, mrg_fld=trim('Faxa_snowc_wiso')//':'//trim('Faxa_snowl_wiso'), &
                  mrg_type='sum_with_weights', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merged sensible heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm , 'Faxa_sen')
       call addfld_aoflux('Faox_sen')
       call addfld_from(compice , 'Fioi_melth')
       call addfld_to(compocn , 'Foxx_sen')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
          call addmrg_to(compocn, 'Foxx_sen', &
               mrg_from=compmed, mrg_fld='Faox_sen', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_lat' )
       call addfld_aoflux( 'Faox_lat' )
       call addfld_aoflux( 'Faox_evap')
       call addfld_to(compocn, 'Foxx_lat' )
       call addfld_to(compocn, 'Foxx_evap')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat', rc=rc)) then
          call addmrg_to(compocn, 'Foxx_lat', &
               mrg_from=compmed, mrg_fld='Faox_lat', mrg_type='merge', mrg_fracname='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc)) then
          call addmrg_to(compocn, 'Foxx_evap', &
               mrg_from=compmed, mrg_fld='Faox_evap', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_aoflux( 'Faox_lat_wiso' )
          call addfld_to(compocn, 'Foxx_lat_wiso' )
       else
          if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat_wiso', rc=rc)) then
             call addmrg_to(compocn, 'Foxx_lat_wiso', &
                  mrg_from=compmed, mrg_fld='Faox_lat_wiso', mrg_type='merge', mrg_fracname='ofrac')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: wind speed squared at 10 meters from med
    ! ---------------------------------------------------------------------
    ! Note that this is a field output by the atm/ocn flux computation
    ! If the aoflux grid is ogrid - then nothing needs to be done to send to the ocean
    ! All other mappings are set in med_phases_aoflux_mod.F90
    if (phase == 'advertise') then
       call addfld_aoflux( 'So_duu10n')
       call addfld_to(compocn, 'So_duu10n')
    else
       if (fldchk(is_local%wrap%FBExp(compocn), 'So_duu10n', rc=rc)) then
          call addmrg_to(compocn, 'So_duu10n', mrg_from=compmed, mrg_fld='So_duu10n', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_pslv')
       call addfld_to(compocn, 'Sa_pslv')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Sa_pslv', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Sa_pslv', rc=rc)) then
          call addmap_from(compatm, 'Sa_pslv', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmap_from(compatm, 'Sa_pslv', compice, mapbilnr, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Sa_pslv', &
               mrg_from=compatm, mrg_fld='Sa_pslv', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: prognostic CO2 at the lowest atm model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_co2prog')
       call addfld_to(compocn, 'Sa_co2prog')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Sa_co2prog', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Sa_co2prog', rc=rc)) then
          call addmap_from(compatm, 'Sa_co2prog', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Sa_co2prog', mrg_from=compatm, mrg_fld='Sa_co2prog', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: diagnostic CO2 at the lowest atm model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_co2diag')
       call addfld_to(compocn, 'Sa_co2diag')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Sa_co2diag', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Sa_co2diag', rc=rc)) then
          call addmap_from(compatm, 'Sa_co2diag', compocn, mapbilnr, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Sa_co2diag', mrg_from=compatm, mrg_fld='Sa_co2diag', mrg_type='copy')
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
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Faxa_bcph')
       call addfld_from(compatm, 'Faxa_bcph')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcph', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_bcph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_bcph', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_bcph', &
               mrg_from=compatm, mrg_fld='Faxa_bcph', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Faxa_ocph')
       call addfld_from(compatm, 'Faxa_ocph')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocph', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_ocph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_ocph', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_ocph', &
               mrg_from=compatm, mrg_fld='Faxa_ocph', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Faxa_dstwet')
       call addfld_from(compatm, 'Faxa_dstwet')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_dstwet', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstwet', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_dstwet', &
               mrg_from=compatm, mrg_fld='Faxa_dstwet', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Faxa_dstdry')
       call addfld_from(compatm, 'Faxa_dstdry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_dstdry', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstdry', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_dstdry', &
               mrg_from=compatm, mrg_fld='Faxa_dstdry', mrg_type='copy_with_weights', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Faxa_ndep')
       call addfld_from(compatm, 'Faxa_ndep')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ndep', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_ndep', rc=rc)) then
          call addmap_from(compatm, 'Faxa_ndep', compocn, mapconsf, 'one', atm2ocn_map)
          call addmrg_to(compocn, 'Faxa_ndep', &
               mrg_from=compatm, mrg_fld='Faxa_ndep', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: enthalpy from atm rain, snow, evaporation
    ! to ocn: enthalpy from liquid and ice river runoff
    ! to ocn: enthalpy from ice melt
    ! ---------------------------------------------------------------------
    ! Note - do not need to add addmap or addmrg for the following since they
    ! will be computed directly in med_phases_prep_ocn
    if (phase == 'advertise') then
       call addfld_to(compocn, 'Foxx_hrain')
       call addfld_to(compocn, 'Foxx_hsnow')
       call addfld_to(compocn, 'Foxx_hevap')
       call addfld_to(compocn, 'Foxx_hcond')
       call addfld_to(compocn, 'Foxx_hrofl')
       call addfld_to(compocn, 'Foxx_hrofi')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merge zonal and meridional surface stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_to(compocn , 'Foxx_taux')
       call addfld_from(compice , 'Fioi_taux')
       call addfld_aoflux('Faox_taux')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_taux', rc=rc)) then
          if (fldchk(is_local%wrap%FBimp(compice,compice), 'Fioi_taux', rc=rc)) then
             call addmap_from(compice, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
             call addmrg_to(compocn, 'Foxx_taux', &
                  mrg_from=compice, mrg_fld='Fioi_taux', mrg_type='merge', mrg_fracname='ifrac')
          end if
          call addmrg_to(compocn, 'Foxx_taux', &
               mrg_from=compmed, mrg_fld='Faox_taux', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_to(compocn , 'Foxx_tauy')
       call addfld_from(compice , 'Fioi_tauy')
       call addfld_aoflux('Faox_tauy')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_tauy', rc=rc)) then
          if (fldchk(is_local%wrap%FBimp(compice,compice), 'Fioi_tauy', rc=rc)) then
             call addmap_from(compice, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
             call addmrg_to(compocn, 'Foxx_tauy', &
                  mrg_from=compice, mrg_fld='Fioi_tauy', mrg_type='merge', mrg_fracname='ifrac')
          end if
          call addmrg_to(compocn, 'Foxx_tauy', &
               mrg_from=compmed, mrg_fld='Faox_tauy', mrg_type='merge', mrg_fracname='ofrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: water flux due to melting ice from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice , 'Fioi_meltw')
       call addfld_to(compocn , 'Fioi_meltw')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_meltw', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_meltw', rc=rc)) then
          call addmap_from(compice, 'Fioi_meltw',    compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_meltw', &
               mrg_from=compice, mrg_fld='Fioi_meltw', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if
    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compice , 'Fioi_meltw_wiso')
          call addfld_to(compocn , 'Fioi_meltw_wiso')
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_meltw_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_meltw_wiso', rc=rc)) then
             call addmap_from(compice, 'Fioi_meltw_wiso',    compocn,  mapfcopy, 'unset', 'unset')
             call addmrg_to(compocn, 'Fioi_meltw_wiso', &
                  mrg_from=compice, mrg_fld='Fioi_meltw_wiso', mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: heat flux from melting ice from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Fioi_melth')
       call addfld_to(compocn, 'Fioi_melth')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_melth', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth', rc=rc)) then
          call addmap_from(compice, 'Fioi_melth', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_melth', &
               mrg_from=compice, mrg_fld='Fioi_melth', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: salt flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Fioi_salt')
       call addfld_to(compocn, 'Fioi_salt')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_salt', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_salt', rc=rc)) then
          call addmap_from(compice, 'Fioi_salt', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_salt', &
               mrg_from=compice, mrg_fld='Fioi_salt', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: hydrophylic black carbon deposition flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Fioi_bcphi')
       call addfld_to(compocn, 'Fioi_bcphi')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_bcphi', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_bcphi', rc=rc)) then
          call addmap_from(compice, 'Fioi_bcphi', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_bcphi', &
               mrg_from=compice, mrg_fld='Fioi_bcphi', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: hydrophobic black carbon deposition flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Fioi_bcpho')
       call addfld_to(compocn, 'Fioi_bcpho')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_bcpho', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_bcpho', rc=rc)) then
          call addmap_from(compice, 'Fioi_bcpho', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_bcpho', &
               mrg_from=compice, mrg_fld='Fioi_bcpho', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ocn: dust flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Fioi_flxdst')
       call addfld_to(compocn, 'Fioi_flxdst')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_flxdst', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_flxdst', rc=rc)) then
          call addmap_from(compice, 'Fioi_flxdst', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compocn, 'Fioi_flxdst', &
               mrg_from=compice, mrg_fld='Fioi_flxdst', mrg_type='copy_with_weights', mrg_fracname='ifrac')
       end if
    end if

    !-----------------------------
    ! to ocn: liquid runoff from rof originating from lnd
    ! to ocn: liquid runoff from rof originating from glc
    ! to ocn: ice runoff from rof originating from lnd
    ! to ocn: ice runoff from rof originating from glc
    ! to ocn: waterflux back to ocn due to flooding from rof
    !-----------------------------

    if (phase == 'advertise') then
       ! Note that Flrr_flood below needs to be added to
       ! fldlistFr(comprof) in order to be mapped correctly to the ocean but the ocean
       ! does not receive it so it is advertised but it will not be connected
       call addfld_from(comprof, 'Forr_rofl')
       call addfld_from(comprof, 'Forr_rofi')
       call addfld_from(comprof, 'Forr_rofl_glc')
       call addfld_from(comprof, 'Forr_rofi_glc')
       call addfld_to(compocn, 'Foxx_rofl')
       call addfld_to(compocn, 'Foxx_rofi')
       call addfld_to(compocn, 'Forr_rofl_glc')
       call addfld_to(compocn, 'Forr_rofi_glc')
       call addfld_to(compocn, 'Flrr_flood')
    else
      ! Liquid runoff from land and glc - mapping
      if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl' , rc=rc)) then
        if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rc=rc)) then
          if (trim(rof2ocn_liq_rmap) == 'unset') then
            call addmap_from(comprof, 'Forr_rofl', compocn, mapconsd, 'one', 'unset')
          else
            call addmap_from(comprof, 'Forr_rofl', compocn, map_rof2ocn_liq, 'none', rof2ocn_liq_rmap)
          end if
        end if
      end if
      if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood', rc=rc)) then
        if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rc=rc)) then
          call addmap_from(comprof, 'Flrr_flood', compocn, mapconsd, 'one', rof2ocn_map)
        end if
      end if
      if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl_glc', rc=rc)) then
        if (fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc', rc=rc) .or. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl', rc=rc)) then
          if (trim(rof2ocn_liq_rmap) == 'unset') then
            call addmap_from(comprof, 'Forr_rofl_glc', compocn, mapconsd, 'one', 'unset')
          else
            call addmap_from(comprof, 'Forr_rofl_glc', compocn, map_rof2ocn_liq, 'none', rof2ocn_liq_rmap)
          end if
        end if
      end if

      ! Liquid runoff from land and glc - merging
      if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rc=rc)) then
        mrgfld_source = 'Forr_rofl'
        if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood', rc=rc)) then
          mrgfld_source = trim(mrgfld_source) //':Flrr_flood'
        end if
        if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl_glc', rc=rc)) then
          mrgfld_source = trim(mrgfld_source) //':Forr_rofl_glc'
        end if
        call addmrg_to(compocn, 'Foxx_rofl', mrg_from=comprof, mrg_fld=trim(mrgfld_source), mrg_type='sum')
      end if

      ! Frozen runoff from land and glc - mapping
      if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi' , rc=rc)) then
        if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , rc=rc)) then
          if (trim(rof2ocn_ice_rmap) == 'unset') then
            call addmap_from(comprof, 'Forr_rofi', compocn, mapconsd, 'one', 'unset')
          else
            call addmap_from(comprof, 'Forr_rofi', compocn, map_rof2ocn_ice, 'none', rof2ocn_ice_rmap)
          end if
        end if
      end if
      if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi_glc', rc=rc)) then
        if (fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc', rc=rc) .or. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi', rc=rc)) then
          if (trim(rof2ocn_ice_rmap) == 'unset') then
            call addmap_from(comprof, 'Forr_rofi_glc', compocn, mapconsd, 'one', 'unset')
          else
            call addmap_from(comprof, 'Forr_rofi_glc', compocn, map_rof2ocn_ice, 'none', rof2ocn_ice_rmap)
          end if
        end if
      end if

      ! Frozen runoff from land and glc - merging
      if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , rc=rc)) then
        mrgfld_source = 'Forr_rofi'
        if (fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi_glc', rc=rc)) then
          mrgfld_source = trim(mrgfld_source) //':Forr_rofi_glc'
        end if
        call addmrg_to(compocn, 'Foxx_rofi', mrg_from=comprof, mrg_fld=trim(mrgfld_source), mrg_type='sum')
      end if
    end if

    !-----------------------------
    ! from wav: for daily averaged fields for
    ! output to auxiliary file only
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_ustokes_avg')
       call addfld_from(compwav, 'Sw_vstokes_avg')
       call addfld_from(compwav, 'Sw_hs_avg')
       call addfld_from(compwav, 'Sw_phs0_avg')
       call addfld_from(compwav, 'Sw_phs1_avg')
       call addfld_from(compwav, 'Sw_pdir0_avg')
       call addfld_from(compwav, 'Sw_pdir1_avg')
       call addfld_from(compwav, 'Sw_pTm10_avg')
       call addfld_from(compwav, 'Sw_pTm11_avg')
       call addfld_from(compwav, 'Sw_Tm1_avg')
       call addfld_from(compwav, 'Sw_thm_avg')
       call addfld_from(compwav, 'Sw_thp0_avg')
       call addfld_from(compwav, 'Sw_fp0_avg')
       call addfld_from(compwav, 'Sw_u_avg')
       call addfld_from(compwav, 'Sw_v_avg')
       call addfld_from(compwav, 'Sw_tusx_avg')
       call addfld_from(compwav, 'Sw_tusy_avg')
    end if

    !-----------------------------
    ! to ocn: Langmuir multiplier from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_lamult')
       call addfld_to(compocn, 'Sw_lamult')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_lamult', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_lamult', rc=rc)) then
          call addmap_from(compwav, 'Sw_lamult', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_lamult', mrg_from=compwav, mrg_fld='Sw_lamult', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ocn: Stokes drift u component from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_ustokes')
       call addfld_to(compocn, 'Sw_ustokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_ustokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_ustokes', rc=rc)) then
          call addmap_from(compwav, 'Sw_ustokes', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_ustokes', mrg_from=compwav, mrg_fld='Sw_ustokes', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ocn: Stokes drift v component from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_vstokes')
       call addfld_to(compocn, 'Sw_vstokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_vstokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_vstokes', rc=rc)) then
          call addmap_from(compwav, 'Sw_vstokes', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_vstokes', mrg_from=compwav, mrg_fld='Sw_vstokes', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ocn: Stokes drift depth from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_hstokes')
       call addfld_to(compocn, 'Sw_hstokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_hstokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_hstokes', rc=rc)) then
          call addmap_from(compwav, 'Sw_hstokes', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_hstokes', mrg_from=compwav, mrg_fld='Sw_hstokes', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ocn: Partitioned stokes drift components in x-direction
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_pstokes_x')
       call addfld_to(compocn, 'Sw_pstokes_x')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_pstokes_x', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_pstokes_x', rc=rc)) then
          call addmap_from(compwav, 'Sw_pstokes_x', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_pstokes_x', mrg_from=compwav, mrg_fld='Sw_pstokes_x', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ocn: Partitioned stokes drift components in y-direction
    !-----------------------------
    if (phase == 'advertise') then
       call addfld_from(compwav, 'Sw_pstokes_y')
       call addfld_to(compocn, 'Sw_pstokes_y')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_pstokes_y', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_pstokes_y', rc=rc)) then
          call addmap_from(compwav, 'Sw_pstokes_y', compocn,  mapbilnr_nstod, 'one', wav2ocn_map)
          call addmrg_to(compocn, 'Sw_pstokes_y', mrg_from=compwav, mrg_fld='Sw_pstokes_y', mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: downward longwave heat flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_lwdn')
       call addfld_to(compice, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_lwdn', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', rc=rc)) then
          call addmap_from(compatm, 'Faxa_lwdn', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_lwdn', mrg_from=compatm, mrg_fld='Faxa_lwdn', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: downward direct near-infrared incident solar radiation  from atm
    ! to ice: downward direct visible incident solar radiation        from atm
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! to ice: downward Diffuse visible incident solar radiation       from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndr')
       call addfld_to(compice, 'Faxa_swndr')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndr', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_swndr', mrg_from=compatm, mrg_fld='Faxa_swndr', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdr')
       call addfld_to(compice, 'Faxa_swvdr')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdr', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_swvdr', mrg_from=compatm, mrg_fld='Faxa_swvdr', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swndf')
       call addfld_to(compice, 'Faxa_swndf')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swndf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swndf', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_swndf', mrg_from=compatm, mrg_fld='Faxa_swndf', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_swvdf')
       call addfld_to(compice, 'Faxa_swvdf')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc)) then
          call addmap_from(compatm, 'Faxa_swvdf', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_swvdf', mrg_from=compatm, mrg_fld='Faxa_swvdf', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: hydrophylic black carbon dry deposition flux   from atm
    ! to ice: hydrophobic black carbon dry deposition flux   from atm
    ! to ice: hydrophylic black carbon wet deposition flux   from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_bcph')
       call addfld_to(compice, 'Faxa_bcph')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcph', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_bcph', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_bcph', mrg_from=compatm, mrg_fld='Faxa_bcph', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: hydrophylic organic carbon dry deposition flux from atm
    ! to ice: hydrophobic organic carbon dry deposition flux from atm
    ! to ice: hydrophylic organic carbon wet deposition flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_ocph')
       call addfld_to(compice, 'Faxa_ocph')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocph', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocph', rc=rc)) then
          call addmap_from(compatm, 'Faxa_ocph', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_ocph', mrg_from=compatm, mrg_fld='Faxa_ocph', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: dust wet deposition flux (size 1) from atm
    ! to ice: dust wet deposition flux (size 2) from atm
    ! to ice: dust wet deposition flux (size 3) from atm
    ! to ice: dust wet deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_dstwet')
       call addfld_to(compice, 'Faxa_dstwet')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstwet', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_dstwet', mrg_from=compatm, mrg_fld='Faxa_dstwet', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: dust dry deposition flux (size 1) from atm
    ! to ice: dust dry deposition flux (size 2) from atm
    ! to ice: dust dry deposition flux (size 3) from atm
    ! to ice: dust dry deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_dstdry')
       call addfld_to(compice, 'Faxa_dstdry')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry', rc=rc)) then
          call addmap_from(compatm, 'Faxa_dstdry', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_dstdry', mrg_from=compatm, mrg_fld='Faxa_dstdry', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: convective and large scale precipitation rate water equivalent from atm
    ! to ice: rain and snow rate from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_rainc')
       call addfld_from(compatm, 'Faxa_rainl')
       call addfld_from(compatm, 'Faxa_rain' )
       call addfld_to(compice, 'Faxa_rain' )
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', rc=rc)) then
          call addmap_from(compatm, 'Faxa_rainc', compice, mapconsf, 'one', atm2ice_map)
          call addmap_from(compatm, 'Faxa_rainl', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain', rc=rc)) then
          call addmap_from(compatm, 'Faxa_rain', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_rain', mrg_from=compatm, mrg_fld='Faxa_rain', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Faxa_snowc')
       call addfld_from(compatm, 'Faxa_snowl')
       call addfld_from(compatm, 'Faxa_snow' )
       call addfld_to(compice, 'Faxa_snow' )
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', rc=rc)) then
          call addmap_from(compatm, 'Faxa_snowc', compice, mapconsf, 'one', atm2ice_map)
          call addmap_from(compatm, 'Faxa_snowl', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_snow' , &
               mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', mrg_type='sum')
       else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow', rc=rc)) then
          call addmap_from(compatm, 'Faxa_snow', compice, mapconsf, 'one', atm2ice_map)
          call addmrg_to(compice, 'Faxa_snow', &
               mrg_from=compatm, mrg_fld='Faxa_snow', mrg_type='copy')
       end if
    end if

    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compatm, 'Faxa_rainc_wiso')
          call addfld_from(compatm, 'Faxa_rainl_wiso')
          call addfld_from(compatm, 'Faxa_rain_wiso' )
          call addfld_to(compice, 'Faxa_rain_wiso' )
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain_wiso' , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', rc=rc)) then
             call addmap_from(compatm, 'Faxa_rainc_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmap_from(compatm, 'Faxa_rainl_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmrg_to(compice, 'Faxa_rain_wiso' , &
                  mrg_from=compatm, mrg_fld='Faxa_rainc_wiso:Faxa_rainl_wiso', mrg_type='sum')
          else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain_wiso', rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain_wiso', rc=rc)) then
             call addmap_from(compatm, 'Faxa_rain_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmrg_to(compice, 'Faxa_rain_wiso', &
                  mrg_from=compatm, mrg_fld='Faxa_rain_wiso', mrg_type='copy')
          end if
       end if

       if (phase == 'advertise') then
          call addfld_from(compatm, 'Faxa_snowc_wiso')
          call addfld_from(compatm, 'Faxa_snowl_wiso')
          call addfld_from(compatm, 'Faxa_snow_wiso' )
          call addfld_to(compice, 'Faxa_snow_wiso' )
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc_wiso', rc=rc)) then
             call addmap_from(compatm, 'Faxa_snowc_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmap_from(compatm, 'Faxa_snowl_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmrg_to(compice, 'Faxa_snow_wiso' , &
                  mrg_from=compatm, mrg_fld='Faxa_snowc_wiso:Faxa_snowl_wiso', mrg_type='sum')
          else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow_wiso', rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow_wiso', rc=rc)) then
             call addmap_from(compatm, 'Faxa_snow_wiso', compice, mapconsf, 'one', atm2ice_map)
             call addmrg_to(compice, 'Faxa_snow_wiso', mrg_from=compatm, mrg_fld='Faxa_snow_wiso', mrg_type='copy')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: height at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_z')
       call addfld_to(compice, 'Sa_z')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_z', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_z', rc=rc)) then
          call addmap_from(compatm, 'Sa_z', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_z', mrg_from=compatm, mrg_fld='Sa_z', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: pressure at the lowest model level fromatm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_pbot')
       call addfld_to(compice, 'Sa_pbot')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_pbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_pbot', rc=rc)) then
          call addmap_from(compatm, 'Sa_pbot', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_pbot', mrg_from=compatm, mrg_fld='Sa_pbot', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: temperature at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_tbot')
       call addfld_to(compice, 'Sa_tbot')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_tbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_tbot', rc=rc)) then
          call addmap_from(compatm, 'Sa_tbot', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_tbot', mrg_from=compatm, mrg_fld='Sa_tbot', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: potential temperature at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_ptem')
       call addfld_to(compice, 'Sa_ptem')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_ptem', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_ptem', rc=rc)) then
          call addmap_from(compatm, 'Sa_ptem', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_ptem', mrg_from=compatm, mrg_fld='Sa_ptem', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: density at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_dens')
       call addfld_to(compice, 'Sa_dens')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_dens', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_dens', rc=rc)) then
          call addmap_from(compatm, 'Sa_dens', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_dens', mrg_from=compatm, mrg_fld='Sa_dens', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: zonal wind at the lowest model level from atm
    ! to ice: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_u')
       call addfld_to(compice, 'Sa_u')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u', rc=rc)) then
          if (mapuv_with_cart3d) then
             call addmap_from(compatm, 'Sa_u', compice, mappatch_uv3d, 'one', atm2ice_map)
          else
             call addmap_from(compatm, 'Sa_u', compice, mappatch, 'one', atm2ice_map)
          end if
          call addmrg_to(compice, 'Sa_u', mrg_from=compatm, mrg_fld='Sa_u', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_v')
       call addfld_to(compice, 'Sa_v')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v', rc=rc)) then
          if (mapuv_with_cart3d) then
             call addmap_from(compatm, 'Sa_v', compice, mappatch_uv3d, 'one', atm2ice_map)
          else
             call addmap_from(compatm, 'Sa_v', compice, mappatch, 'one', atm2ice_map)
          end if
          call addmrg_to(compice, 'Sa_v', mrg_from=compatm, mrg_fld='Sa_v', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: specific humidity at the lowest model level from atm
    ! to ice: specific humidity for water isotopes at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_shum')
       call addfld_to(compice, 'Sa_shum')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_shum', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_shum', rc=rc)) then
          call addmap_from(compatm, 'Sa_shum', compice, mapbilnr, 'one', atm2ice_map)
          call addmrg_to(compice, 'Sa_shum', mrg_from=compatm, mrg_fld='Sa_shum', mrg_type='copy')
       end if
    end if
    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compatm, 'Sa_shum_wiso')
          call addfld_to(compice, 'Sa_shum_wiso')
       else
          if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_shum_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_shum_wiso', rc=rc)) then
             call addmap_from(compatm, 'Sa_shum_wiso', compice, mapbilnr, 'one', atm2ice_map)
             call addmrg_to(compice, 'Sa_shum_wiso', mrg_from=compatm, mrg_fld='Sa_shum_wiso', mrg_type='copy')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: sea surface temperature from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_t')
       call addfld_to(compice, 'So_t')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap_from(compocn, 'So_t', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: sea surface salinity from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_s')
       call addfld_to(compice, 'So_s')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_s', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_s', rc=rc)) then
          call addmap_from(compocn, 'So_s', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_s', mrg_from=compocn, mrg_fld='So_s', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: zonal sea water velocity from ocn
    ! to ice: meridional sea water velocity from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_u')
       call addfld_to(compice, 'So_u')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_u', rc=rc)) then
          call addmap_from(compocn, 'So_u', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_u', mrg_from=compocn, mrg_fld='So_u', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_v')
       call addfld_to(compice, 'So_v')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_v', rc=rc)) then
          call addmap_from(compocn, 'So_v', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_v', mrg_from=compocn, mrg_fld='So_v', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: zonal sea surface slope from ocn
    ! to ice: meridional sea surface slope from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_dhdx')
       call addfld_to(compice, 'So_dhdx')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_dhdx', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_dhdx', rc=rc)) then
          call addmap_from(compocn, 'So_dhdx', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_dhdx', mrg_from=compocn, mrg_fld='So_dhdx', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_dhdy')
       call addfld_to(compice, 'So_dhdy')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_dhdy', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_dhdy', rc=rc)) then
          call addmap_from(compocn, 'So_dhdy', compice, mapfcopy , 'unset', 'unset')
          call addmrg_to(compice, 'So_dhdy', mrg_from=compocn, mrg_fld='So_dhdy', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to ice: ocean melt and freeze potential from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'Fioo_q')
       call addfld_to(compice, 'Fioo_q')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'Fioo_q', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'Fioo_q', rc=rc)) then
          call addmap_from(compocn, 'Fioo_q', compice,  mapfcopy, 'unset', 'unset')
          call addmrg_to(compice, 'Fioo_q', mrg_from=compocn, mrg_fld='Fioo_q', mrg_type='copy')
       end if
    end if
    !-----------------------------
    ! to ice: Ratio of ocean surface level abund. H2_16O/H2O/Rstd from ocean
    !-----------------------------
    if (flds_wiso) then
       if (phase == 'advertise') then
          call addfld_from(compocn, 'So_roce_wiso')
          call addfld_to(compice, 'So_roce_wiso')
       else
          if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_roce_wiso', rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compice)         , 'So_roce_wiso', rc=rc)) then
             call addmap_from(compocn, 'So_roce_wiso', compice,  mapfcopy, 'unset', 'unset')
             call addmrg_to(compice, 'So_roce_wiso', mrg_from=compocn, mrg_fld='So_roce_wiso', mrg_type='copy')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: wave elevation spectrum (field with ungridded dimensions)
    ! ---------------------------------------------------------------------
    if (wav_coupling_to_cice) then
       if (phase == 'advertise') then
          call addfld_from(compwav, 'Sw_elevation_spectrum')
          call addfld_to(compice, 'Sw_elevation_spectrum')
       else
          if ( fldchk(is_local%wrap%FBExp(compice)        , 'Sw_elevation_spectrum', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_elevation_spectrum', rc=rc)) then
             call addmap_from(compwav, 'Sw_elevation_spectrum', compice, mapbilnr_nstod, 'one', wav2ocn_map)
             call addmrg_to(compice, 'Sw_elevation_spectrum', &
                  mrg_from=compwav, mrg_fld='Sw_elevation_spectrum', mrg_type='copy')
          end if
       end if
    end if

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    !----------------------------------------------------------
    ! to wav: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compice, 'Si_ifrac')
       call addfld_to(compwav, 'Si_ifrac')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Si_ifrac', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_ifrac', rc=rc)) then
             ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
          call addmap_from(compice, 'Si_ifrac', compwav, mapbilnr, 'one', ice2wav_map)
          call addmrg_to(compwav, 'Si_ifrac', mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')
       end if
    end if
    !----------------------------------------------------------
    ! to wav: ice thickness from ice
    !----------------------------------------------------------
    if (wav_coupling_to_cice) then
       if (phase == 'advertise') then
          call addfld_from(compice, 'Si_thick')
          call addfld_to(compwav, 'Si_thick')
       else
          if (fldchk(is_local%wrap%FBexp(compwav)         , 'Si_thick', rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_thick', rc=rc)) then
             call addmap_from(compice, 'Si_thick', compwav, mapbilnr, 'one', ice2wav_map)
             call addmrg_to(compwav, 'Si_thick', mrg_from=compice, mrg_fld='Si_thick', mrg_type='copy')
          end if
       end if
    end if
    !----------------------------------------------------------
    ! to wav: ice floe diameter from ice
    !----------------------------------------------------------
    if (wav_coupling_to_cice) then
       if (phase == 'advertise') then
          call addfld_from(compice, 'Si_floediam')
          call addfld_to(compwav, 'Si_floediam')
       else
          if (fldchk(is_local%wrap%FBexp(compwav)         , 'Si_floediam', rc=rc) .and. &
              fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_floediam', rc=rc)) then
             call addmap_from(compice, 'Si_floediam', compwav, mapbilnr, 'one', ice2wav_map)
             call addmrg_to(compwav, 'Si_floediam', mrg_from=compice, mrg_fld='Si_floediam', mrg_type='copy')
          end if
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to wav: ocean surface temperature from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_t')
       call addfld_to(compwav, 'So_t')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_t', rc=rc)) then
          call addmap_from(compocn, 'So_t', compwav, mapbilnr, 'one', ocn2wav_map)
          call addmrg_to(compwav, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if
    ! ---------------------------------------------------------------------
    ! to wav: ocean currents from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_u')
       call addfld_to(compwav, 'So_u')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_u', rc=rc)) then
          ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
          call addmap_from(compocn, 'So_u', compwav, mapbilnr, 'one', ocn2wav_map)
          call addmrg_to(compwav, 'So_u', mrg_from=compocn, mrg_fld='So_u', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_v')
       call addfld_to(compwav, 'So_v')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_v', rc=rc)) then
          ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
          call addmap_from(compocn, 'So_v', compwav, mapbilnr, 'one', ocn2wav_map)
          call addmrg_to(compwav, 'So_v', mrg_from=compocn, mrg_fld='So_v', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: ocean boundary layer depth from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compocn, 'So_bldepth')
       call addfld_to(compwav, 'So_bldepth')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_bldepth', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_bldepth', rc=rc)) then
          ! By default will be using a custom map - but if one is not available, use a generated bilinear instead
          call addmap_from(compocn, 'So_bldepth', compwav, mapbilnr, 'one', ocn2wav_map)
          call addmrg_to(compwav, 'So_bldepth', mrg_from=compocn, mrg_fld='So_bldepth', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: zonal and meridional winds at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_u')
       call addfld_to(compwav, 'Sa_u')
       call addfld_from(compatm, 'Sa_u10m')
       call addfld_to(compwav, 'Sa_u10m')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u', rc=rc)) then
          call addmap_from(compatm, 'Sa_u', compwav, mapbilnr, 'one', atm2wav_map)
          call addmrg_to(compwav, 'Sa_u', mrg_from=compatm, mrg_fld='Sa_u', mrg_type='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_u10m', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u10m', rc=rc)) then
          call addmap_from(compatm, 'Sa_u10m', compwav, mapbilnr, 'one', atm2wav_map)
          call addmrg_to(compwav, 'Sa_u10m', mrg_from=compatm, mrg_fld='Sa_u10m', mrg_type='copy')
       end if
    end if
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_v')
       call addfld_to(compwav, 'Sa_v')
       call addfld_from(compatm, 'Sa_v10m')
       call addfld_to(compwav, 'Sa_v10m')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v', rc=rc)) then
          call addmap_from(compatm, 'Sa_v', compwav, mapbilnr, 'one', atm2wav_map)
          call addmrg_to(compwav, 'Sa_v', mrg_from=compatm, mrg_fld='Sa_v', mrg_type='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_v10m', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v10m', rc=rc)) then
          call addmap_from(compatm, 'Sa_v10m', compwav, mapbilnr, 'one', atm2wav_map)
          call addmrg_to(compwav, 'Sa_v10m', mrg_from=compatm, mrg_fld='Sa_v10m', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: temperature at lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(compatm, 'Sa_tbot')
       call addfld_to(compwav, 'Sa_tbot')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_tbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_tbot', rc=rc)) then
          call addmap_from(compatm, 'Sa_tbot', compwav, mapbilnr, 'one', atm2wav_map)
          call addmrg_to(compwav, 'Sa_tbot', mrg_from=compatm, mrg_fld='Sa_tbot', mrg_type='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: zonal and meridional wind stress
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_to(compwav , 'Fwxx_taux')
       call addfld_to(compwav , 'Fwxx_tauy')
    end if

    !=====================================================================
    ! FIELDS TO RIVER (comprof)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to rof: liquid and ice from glc
    ! ---------------------------------------------------------------------
    do ns = 1, is_local%wrap%num_icesheets
       if (phase == 'advertise') then
          call addfld_from(compglc(ns), 'Fgrg_rofl')
          call addfld_from(compglc(ns), 'Fgrg_rofi')
          call addfld_to(comprof, 'Fgrg_rofl')
          call addfld_to(comprof, 'Fgrg_rofi')
       else
          ! Note: we are assuming that the rof mesh has a mask of one everywhere
          if ( fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fgrg_rofl', rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)                 , 'Fgrg_rofl', rc=rc)) then
             call addmap_from(compglc(ns), 'Fgrg_rofl', comprof, mapconsd, 'gfrac' , 'unset')
             ! Custom merge in med_phases_prep_rof
          end if
          if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fgrg_rofi', rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(comprof)                 , 'Fgrg_rofi', rc=rc)) then
             call addmap_from(compglc(ns), 'Fgrg_rofi', comprof, mapconsd, 'gfrac', 'unset')
             ! Custom merge in med_phases_prep_rof
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: liquid and ice from glc water isoptopes
    ! ---------------------------------------------------------------------
    do ns = 1, is_local%wrap%num_icesheets
       if (phase == 'advertise') then
          call addfld_from(compglc(ns), 'Fgrg_rofl_wiso')
          call addfld_from(compglc(ns), 'Fgrg_rofi_wiso')
          call addfld_to(comprof, 'Fgrg_rofl_wiso')
          call addfld_to(comprof, 'Fgrg_rofi_wiso')
       else
          if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fgrg_rofl_wiso' , rc=rc)) then
             call addmap_from(compglc(ns), 'Fgrg_rofl_wiso', comprof, mapconsd, 'one' , 'unset')
             ! TODO: implement custom merge
          end if
          if (fldchk(is_local%wrap%FBImp(compglc(ns), compglc(ns)), 'Fgrg_rofi_wiso' , rc=rc)) then
             call addmap_from(compglc(ns), 'Fgrg_rofi_wiso', comprof, mapconsd, 'one', 'unset')
             ! TODO: implement custom merge
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid surface)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Flrl_rofsur')
       call addfld_to(comprof, 'Flrl_rofsur')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofsur', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofsur', rc=rc)) then
          call addmap_from(complnd, 'Flrl_rofsur', comprof, mapconsf, map_fracname_lnd2rof, 'unset')
          call addmrg_to(comprof, 'Flrl_rofsur', &
               mrg_from=complnd, mrg_fld='Flrl_rofsur', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2rof)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (ice surface)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Flrl_rofi')
       call addfld_to(comprof, 'Flrl_rofi')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofi', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofi', rc=rc)) then
          call addmap_from(complnd, 'Flrl_rofi', comprof, mapconsf, map_fracname_lnd2rof, 'unset')
          call addmrg_to(comprof, 'Flrl_rofi', &
               mrg_from=complnd, mrg_fld='Flrl_rofi', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2rof)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid glacier, wetland, and lake)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Flrl_rofgwl')
       call addfld_to(comprof, 'Flrl_rofgwl')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofgwl', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofgwl', rc=rc)) then
          call addmap_from(complnd, 'Flrl_rofgwl', comprof, mapconsf, map_fracname_lnd2rof, 'unset')
          call addmrg_to(comprof, 'Flrl_rofgwl', &
               mrg_from=complnd, mrg_fld='Flrl_rofgwl', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2rof)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid subsurface)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Flrl_rofsub')
       call addfld_to(comprof, 'Flrl_rofsub')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofsub', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofsub', rc=rc)) then
          call addmap_from(complnd, 'Flrl_rofsub', comprof, mapconsf, map_fracname_lnd2rof, 'unset')
          call addmrg_to(comprof, 'Flrl_rofsub', &
               mrg_from=complnd, mrg_fld='Flrl_rofsub', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2rof)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to rof: irrigation flux from land (withdrawal from rivers)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld_from(complnd, 'Flrl_irrig')
       call addfld_to(comprof, 'Flrl_irrig')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_irrig', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_irrig', rc=rc)) then
          call addmap_from(complnd, 'Flrl_irrig', comprof, mapconsf, map_fracname_lnd2rof, 'unset')
          call addmrg_to(comprof, 'Flrl_irrig', &
               mrg_from=complnd, mrg_fld='Flrl_irrig', mrg_type='copy_with_weights', mrg_fracname=mrg_fracname_lnd2rof)
       end if
    end if

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
       call addfld_from(complnd, 'Sl_tsrf_elev')   ! surface temperature of glacier (1->glc_nec+1)
       call addfld_from(complnd, 'Sl_topo_elev')   ! surface heights of glacier     (1->glc_nec+1)
       call addfld_from(complnd, 'Flgl_qice_elev') ! glacier ice flux               (1->glc_nec+1)
       do ns = 1,is_local%wrap%num_icesheets
          call addfld_to(compglc(ns), 'Sl_tsrf')
          call addfld_to(compglc(ns), 'Flgl_qice')
       end do
    else
       ! custom mapping, accumulation and merging will be done in prep_glc_mod.F90
       do ns = 1,is_local%wrap%num_icesheets
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Flgl_qice_elev', rc=rc)) then
             call addmap_from(complnd, 'Flgl_qice_elev', compglc(ns), mapbilnr, map_fracname_lnd2glc, 'unset')
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Sl_tsrf_elev'  , rc=rc)) then
             call addmap_from(complnd, 'Sl_tsrf_elev', compglc(ns), mapbilnr, map_fracname_lnd2glc, 'unset')
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd) , 'Sl_topo_elev'  , rc=rc)) then
             ! This is needed just for mapping to glc - but is not sent as a field
             call addmap_from(complnd, 'Sl_topo_elev', compglc(ns), mapbilnr, map_fracname_lnd2glc, 'unset')
          end if
       end do
    end if

    !-----------------------------
    ! to glc: from ocn
    !-----------------------------
    if (is_local%wrap%ocn2glc_coupling) then
       if (phase == 'advertise') then
          call addfld_from(compocn, 'So_t_depth')
          call addfld_from(compocn, 'So_s_depth')
          do ns = 1,is_local%wrap%num_icesheets
             call addfld_to(compglc(ns), 'So_t_depth')
             call addfld_to(compglc(ns), 'So_s_depth')
          end do
       else
          ! custom mapping, accumulation and merging will be done in prep_glc_mod.F90
          ! the following is used to create the route handle
          do ns = 1,is_local%wrap%num_icesheets
             if ( fldchk(is_local%wrap%FBImp(compocn,compocn) , 'So_t_depth', rc=rc)) then
                call addmap_from(compocn, 'So_t_depth', compglc(ns), mapbilnr, 'none', 'unset')
             end if
             if ( fldchk(is_local%wrap%FBImp(compocn,compocn) , 'So_s_depth', rc=rc)) then
                call addmap_from(compocn, 'So_s_depth', compglc(ns), mapbilnr, 'none', 'unset')
             end if
          end do
       end if
    end if

  end subroutine esmFldsExchange_cesm

end module esmFldsExchange_cesm_mod
