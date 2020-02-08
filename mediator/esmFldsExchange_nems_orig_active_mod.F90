module esmFldsExchange_nems_orig_active_mod

  ! *************************************************************************
  ! NOTE: in NEMS, atm->ocn and atm->ice some of the fields are mapped
  ! bilinearly or with patch followed by with near neighbor
  ! interpolation - the bilinear and patch nearest neighbor is not available
  ! in CMEPS yet, so conservative nearest neighbor will be done instead
  ! *************************************************************************

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !
  ! Mediator field naming convention
  !  state-prefix:
  !    first 3 characters: Sa_, Si_, So_ (a => atm, i => ice, o => ocn, x=>mediator)
  !  state-name:
  !    what follows state prefix
  !  flux-prefix:
  !    first 5 characters: Flmn__
  !    lm => between components l and m
  !    n  => computed by component n
  !    example: Fioi => ice/ocn flux computed by ice
  !    example: Faxz => atm/mediator flux computed by atm
  !  flux-name:
  !    what follows flux-prefix
  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------
  ! Merging API: (specifies merges to destination field from source fields)
  !---------------------------------------------------------------------
  ! call med_fldList_AddMrg(fld, fldname, &
  !                         mrg_from1, mrg_fld1, mrg_type1, mrg_fracname1, &
  !                         mrg_from2, mrg_fld2, mrg_type2, mrg_fracname2, &
  !                         mrg_from3, mrg_fld3, mrg_type3, mrg_fracname3, &
  !                         mrg_from4, mrg_fld4, mrg_type4, mrg_fracname4)
  ! fld => fldlistTo(destination component name)
  ! fldname => field name in fldListTo(destination component name)
  ! mrg_fromN => source component name
  ! mrg_fldN  => field name in source component name
  ! mrg_typeN => type of merge
  !              [(copy,copy_with_weights,sum,sum_with_weights,merge]
  ! mrge_fracnameN => if merge has weights - use this fractional field name
  !              [ifrac, lfrac, ofrac, etc\
  !
  !---------------------------------------------------------------------
  ! Mapping API:  (specifies mapping from source fields to a destination field)
  !---------------------------------------------------------------------
  ! call med_fldList_AddMap(flds, fldname, destcomp, maptype, mapnorm)
  !   flds => fldlistFr(source component name)
  !   fldname  => field name in fldlistFr(source component name)
  !   destname => destination component name
  !   maptype  => source to destinationmapping type for fldname
  !               [bilnr,consf,consd,patch,fcopy,nstod,nstod_consd,nstod_consf]
  !   mapnorm  => normalization to use in source -> destination mapping
  !               [unset, one, none, custom, lfrin, lfrac, ifrac, ofrac]
  !---------------------------------------------------------------------
  !
  ! ------------
  ! med -> fv3
  ! ------------
  ! So_t      - sea_surface_temperature           from ocn only
  ! ---       - sea_ice_surface_temperature       from ice only (* should be Si_t but cesm/nems have conflict)
  ! Si_ifrac  - ice_fraction                      from ice only
  ! Si_vice   - mean_ice_volume                   from ice only
  ! Si_vsno   - mean_snow_volume                  from ice only
  ! Faii_evap - mean_evap_rate_atm_into_ice       from_ice_only
  ! Faii_lwup - mean_up_lw_flx_ice                from ice only
  ! Faii_lat  - mean_laten_heat_flx_ice           from ice only
  ! Faii_sen  - mean_sensi_heat_flx_ice           from ice only
  ! Faii_taux - mean_zonal_moment_flx_ice         from ice only
  ! Faii_tauy - mean_merid_moment_flx_ice         from ice only
  ! -------------------------------------------------
  !
  ! ------------
  ! fv3 -> med
  ! ------------
  ! Sa_pbot        - inst_pres_height_surface      to ice
  ! Sa_z           - inst_height_lowest            to ice
  ! Sa_tbot        - inst_temp_height_lowest       to ice
  ! Sa_shum        - inst_spec_humid_height_lowest to ice
  ! Sa_u           - inst_zonal_wind_height_lowest to ice
  ! Sa_v           - inst_merid_wind_height_lowest to ice
  ! Sa_pbot        - inst_pres_height_lowest       to ice
  ! Faxa_lwdn      - mean_down_lw_flx              to ice
  ! Faxa_swvdr     - mean_down_sw_vis_dir_flx      to ice and ocn
  ! Faxa_swvdf     - mean_down_sw_vis_dif_flx      to ice and ocn
  ! Faxa_swndr     - mean_down_sw_ir_dir_flx       to ice and ocn
  ! Faxa_swndf     - mean_down_sw_ir_dif_flx       to ice and ocn
  ! Faxa_taux      - mean_zonal_moment_flx         to ocn
  ! Faxa_tauy      - mean_merid_moment_flx         to ocn
  ! Faxa_sen       - mean_sensi_heat_flx           to ocn
  ! Faxa_lat       - mean_laten_heat_flx           to ocn
  ! Faxa_swdn      - mean_down_sw_flx              to ocn
  ! Faxa_lwnet     - mean_net_lw_flx               to ocn
  ! Faxa_rain      - mean_prec_rate                to ocn
  ! Faxa_snow      - mean_fprec_rate               to ocn
  ! -------------------------------------------------
  !
  ! ------------
  ! med -> mom6 (see *note below)
  ! ------------
  ! Fioi_salt      - mean_salt_rate
  ! Fioi_meltw     - mean_fresh_water_to_ocean_rate
  ! Fioi_melth     - net_heat_flx_to_ocn
  ! Foxx_taux      - mean_zonal_moment_flx
  ! Foxx_tauy      - mean_meridional_moment_flx
  ! Foxx_sen       - mean_sensi_heat_flx
  ! Foxx_evap      - mean_evap_rate
  ! ---            - mean_net_sw_vis_dir_flx (* should be Foxx_swnet_vdr but fd.yaml alias confict)
  ! ---            - mean_net_sw_vis_dif_flx (* should be Foxx_swnet_vdf but fd.yaml alias confict)
  ! ---            - mean_net_sw_ir_dir_flx  (* should be Foxx_swnet_idr but fd.yaml alias confict)
  ! ---            - mean_net_sw_ir_dif_flx  (* should be Foxx_swnet_idf but fd.yaml alias confict)
  ! Faxa_lwnet     - mean_net_lw_flx
  ! Faxa_rain      - mean_prec_rate
  ! Faxa_snow      - mean_fprec_rate
  ! Sa_pslv        - inst_pres_height_surface
  !
  ! ------------
  ! mom6 -> med
  ! ------------
  ! So_omask - ocean_mask
  ! So_t     - sea_surface_temperature
  ! So_s     - s_surf
  ! So_u     - ocn_current_zonal
  ! So_v     - ocn_current_merid
  ! So_dhdx  - sea_surface_slope_zonal
  ! So_dhdy  - sea_surface_slope_merid
  ! Fioo_q   - freezing_melting_potential
  !
  ! ------------
  ! cice -> med (see * note below)
  ! ------------
  ! Si_imask       - ice_mask
  ! Si_ifrac       - ice_fraction
  ! ---            - sea_ice_surface_temperature (* should be Si_t but cesm/nems have conflict)
  ! Si_avsdr       - inst_ice_vis_dir_albedo
  ! Si_anidr       - inst_ice_ir_dir_albedo
  ! Si_avsdf       - inst_ice_vis_dif_albedo
  ! Si_anidf       - inst_ice_ir_dif_albedo
  ! Faii_taux      - stress_on_air_ice_zonal
  ! Faii_tauy      - stress_on_air_ice_merid
  ! Fioi_taux      - stress_on_ocn_ice_zonal
  ! Fioi_tauy      - stress_on_ocn_ice_merid
  ! Fioi_swpen     - mean_sw_pen_to_ocn
  ! ---            - mean_net_sw_vis_dir_flx (* should be Fioi_swpen_ but fd.yaml alias confict)
  ! ---            - mean_net_sw_vis_dif_flx (* should be Fioi_swpen_ but fd.yaml alias confict)
  ! ---            - mean_net_sw_ir_dir_flx  (* should be Fioi_swpen_ but fd.yaml alias confict)
  ! ---            - mean_net_sw_ir_dif_flx  (* should be Fioi_swpen_ but fd.yaml alias confict)
  ! Faii_lwup      - mean_up_lw_flx_ice
  ! Faii_sen       - mean_sensi_heat_flx_atm_into_ice
  ! Faii_lat       - mean_laten_heat_flx_atm_into_ice
  ! Faii_evap      - mean_evap_rate_atm_into_ice
  ! Fioi_meltw     - mean_fresh_water_to_ocean_rate
  ! Fioi_salt      - mean_salt_rate
  ! Fioi_melth     - net_heat_flx_to_ocn
  ! Si_vice        - mean_ice_volume
  ! Si_vsno        - mean_snow_volume
  !
  ! ------------
  ! med -> cice
  ! ------------
  ! Sa_z       - inst_height_lowest
  ! Sa_tbot    - inst_temp_height_lowest
  ! Sa_shum    - inst_spec_humid_height_lowest
  ! Sa_u       - inst_zonal_wind_height_lowest
  ! Sa_v       - inst_merid_wind_height_lowest
  ! Sa_pbot    - inst_pres_height_lowest
  ! Faxa_lwdn  - mean_down_lw_flx
  ! Faxa_swndf - mean_down_sw_vis_dir_flx
  ! Faxa_swvdf - mean_down_sw_vis_dif_flx
  ! Faxa_swndr - mean_down_sw_ir_dir_flx
  ! Faxa_swndf - mean_down_sw_ir_dif_flx
  ! Faxa_rain  - mean_prec_rate
  ! Faxa_snow  - mean_fprec_rate
  ! So_t       - sea_surface_temperature
  ! So_s       - s_surf
  ! So_dhdx    - sea_surface_slope_zonal
  ! So_dhdy    - sea_surface_slope_merid
  ! So_u       - ocn_current_zonal
  ! So_v       - ocn_current_merid
  ! Fioo_q     - freezing_melting_potential
  ! Sa_dens    - air_density_height_lowest
  ! -------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_nems_orig_active

  character(*), parameter :: u_FILE_u = &
       __FILE__

  !================================================================================
contains
  !================================================================================

  subroutine esmFldsExchange_nems_orig_active(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_internalstate_mod , only : InternalState
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : compmed, compatm, compocn, compice, ncomps
    use esmflds               , only : mapbilnr, mapfcopy, mapnstod, mapnstod_consf
    use esmflds               , only : fldListTo, fldListFr
    use med_internalstate_mod , only : mastertask, logunit

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: i, n
    character(len=CL)   :: cvalue
    character(len=CS)   :: name
    character(len=*) , parameter   :: subname='(esmFldsExchange_ufs)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

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
    ! FIELDS TO MEDIATOR (for fractions)
    !=====================================================================

    ! to med: masks from components
    call addfld(fldListFr(compocn)%flds, 'So_omask')
    call addfld(fldListFr(compice)%flds, 'Si_imask')
    call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset')

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    ! to atm: Fractions
    ! the following are computed in med_phases_prep_atm
    call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
    call addfld(fldListTo(compatm)%flds, 'So_ofrac')

    ! to atm: zonal surface stress from ice
    call addfld(fldListFr(compice)%flds, 'Faii_taux')
    call addfld(fldListTo(compatm)%flds, 'Faii_taux')
    call addmap(fldListFr(compice)%flds, 'Faii_taux', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_taux', &
         mrg_from1=compice, mrg_fld1='Faii_taux', mrg_type1='copy')

    ! to atm: meridional surface stress from ice
    call addfld(fldListFr(compice)%flds, 'Faii_tauy')
    call addfld(fldListTo(compatm)%flds, 'Faii_tauy')
    call addmap(fldListFr(compice)%flds, 'Faii_tauy', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_tauy', &
         mrg_from1=compice, mrg_fld1='Faii_tauy', mrg_type1='copy')

    ! to atm: surface latent heat flux from ice
    call addfld(fldListFr(compice)%flds, 'Faii_lat')
    call addfld(fldListTo(compatm)%flds, 'Faii_lat')
    call addmap(fldListFr(compice)%flds, 'Faii_lat', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_lat', &
         mrg_from1=compice, mrg_fld1='Faii_lat', mrg_type1='copy')

    ! to atm: surface sensible heat flux from ice
    call addfld(fldListFr(compice)%flds, 'Faii_sen')
    call addfld(fldListTo(compatm)%flds, 'Faii_sen')
    call addmap(fldListFr(compice)%flds, 'Faii_sen', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_sen', &
         mrg_from1=compice, mrg_fld1='Faii_sen', mrg_type1='copy')

    ! to atm: surface upward longwave heat flux from ice
    call addfld(fldListFr(compice)%flds, 'Faii_lwup')
    call addfld(fldListTo(compatm)%flds, 'Faii_lwup')
    call addmap(fldListFr(compice)%flds, 'Faii_lwup', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_lwup', &
         mrg_from1=compice, mrg_fld1='Faii_lwup', mrg_type1='copy')

    ! to atm: evaporation water flux from water from ice
    call addfld(fldListFr(compice)%flds, 'Faii_evap')
    call addfld(fldListTo(compatm)%flds, 'Faii_evap')
    call addmap(fldListFr(compice)%flds, 'Faii_evap', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Faii_evap', &
         mrg_from1=compice, mrg_fld1='Faii_evap', mrg_type1='copy')

    ! to atm: unmerged temperatures from ice
    call addfld(fldListFr(compice)%flds, 'sea_ice_surface_temperature')
    call addfld(fldListTo(compatm)%flds, 'sea_ice_surface_temperature')
    call addmap(fldListFr(compice)%flds, 'sea_ice_surface_temperature', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'sea_ice_surface_temperature', &
         mrg_from1=compice, mrg_fld1='sea_ice_surface_temperature', mrg_type1='copy')

    ! to atm: unmerged temperatures from ocn
    call addfld(fldListFr(compocn)%flds, 'So_t')
    call addfld(fldListTo(compatm)%flds, 'So_t')
    call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapnstod_consf, 'none' )
    call addmrg(fldListTo(compatm)%flds, 'So_t', &
         mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

    ! to atm: mean ice volume per unit area  from ice
    call addfld(fldListFr(compice)%flds, 'Si_vice')
    call addfld(fldListTo(compatm)%flds, 'Si_vice')
    call addmap(fldListFr(compice)%flds, 'Si_vice', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Si_vice', &
         mrg_from1=compice, mrg_fld1='Si_vice', mrg_type1='copy')

    ! to atm: mean snow volume per unit area from ice
    call addfld(fldListFr(compice)%flds, 'Si_vsno')
    call addfld(fldListTo(compatm)%flds, 'Si_vsno')
    call addmap(fldListFr(compice)%flds, 'Si_vsno', compatm, mapnstod_consf, 'ifrac')
    call addmrg(fldListTo(compatm)%flds, 'Si_vsno', &
         mrg_from1=compice, mrg_fld1='Si_vsno', mrg_type1='copy')

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! to ocn: fractional ice coverage wrt ocean from ice
    call addfld(fldListFr(compice)%flds, 'Si_ifrac')
    call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
    call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')

    ! to ocn: merged longwave net heat flux
    call addfld(fldListFr(compatm)%flds, 'Faxa_lwnet')
    call addfld(fldListTo(compocn)%flds, 'Foxx_lwnet')
    call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapnstod_consf, 'one'  )
    call addmrg(fldListTo(compocn)%flds, 'Faxa_lwnet', &
         mrg_from1=compatm, mrg_fld1='Faxa_lwnet', mrg_type1='copy_with_weights', mrg_fracname1='ofrac')

    ! to ocn: net shortwave radiation from med
    ! *** custom merge in med_phases_prep_ocn ***
    call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
    call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
    call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
    call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapnstod_consf, 'one')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapnstod_consf, 'one')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapnstod_consf, 'one')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapnstod_consf, 'one')
    call addfld(fldListFr(compice)%flds, 'mean_net_sw_vis_dir_flx')
    call addfld(fldListFr(compice)%flds, 'mean_net_sw_vis_dif_flx')
    call addfld(fldListFr(compice)%flds, 'mean_net_sw_ir_dir_flx' )
    call addfld(fldListFr(compice)%flds, 'mean_net_sw_ir_dif_flx' )
    call addmap(fldListFr(compice)%flds, 'mean_net_sw_vis_dir_flx', compocn, mapfcopy, 'unset')
    call addmap(fldListFr(compice)%flds, 'mean_net_sw_vis_dif_flx', compocn, mapfcopy, 'unset')
    call addmap(fldListFr(compice)%flds, 'mean_net_sw_ir_dir_flx' , compocn, mapfcopy, 'unset')
    call addmap(fldListFr(compice)%flds, 'mean_net_sw_ir_dif_flx' , compocn, mapfcopy, 'unset')
    call addfld(fldListTo(compocn)%flds, 'mean_net_sw_vis_dir_flx')
    call addfld(fldListTo(compocn)%flds, 'mean_net_sw_vis_dif_flx')
    call addfld(fldListTo(compocn)%flds, 'mean_net_sw_ir_dir_flx' )
    call addfld(fldListTo(compocn)%flds, 'mean_net_sw_ir_dif_flx' )

    !  to ocn: precipitation rate water equivalent from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_rain')
    call addfld(fldListTo(compocn)%flds, 'Faxa_rain')
    call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compocn, mapnstod_consf, 'one')
    call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', &
         mrg_from1=compatm, mrg_fld1='Faxa_rain', mrg_type1='copy_with_weights', mrg_fracname1='ofrac')

    !  to ocn: snow rate water equivalent from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_snow')
    call addfld(fldListTo(compocn)%flds, 'Faxa_snow')
    call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compocn, mapnstod_consf, 'one')
    call addmrg(fldListTo(compocn)%flds, 'Faxa_snow', &
         mrg_from1=compatm, mrg_fld1='Faxa_snow', mrg_type1='copy_with_weights', mrg_fracname1='ofrac')

    ! to ocn: sea level pressure from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
    call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
    call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, mapnstod_consf, 'one') ! bilinear in nems
    call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', &
         mrg_from1=compatm, mrg_fld1='Sa_pslv', mrg_type1='copy')

    ! to ocn: merged sensible heat flux
    call addfld(fldListFr(compatm)%flds, 'Faxa_sen')
    call addfld(fldListTo(compocn)%flds, 'Foxx_sen')
    call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapnstod_consf, 'one')

    ! to ocn: surface latent heat flux and evaporation water flux
    ! *** custom merge in med_phases_prep_ocn ***
    call addfld(fldListFr(compatm)%flds, 'Faxa_lat')
    call addfld(fldListTo(compocn)%flds, 'Foxx_evap')
    call addmap(fldListFr(compatm)%flds, 'Faxa_lat' , compocn, mapnstod_consf, 'one')

    ! to ocn: merge zonal surface stress from ice and atm
    ! *** custom merge in med_phases_prep_ocn ***
    call addfld(fldListFr(compice)%flds, 'Fioi_taux')
    call addfld(fldListFr(compice)%flds, 'Fioi_tauy')
    call addfld(fldListFr(compatm)%flds, 'Faxa_taux')
    call addfld(fldListFr(compatm)%flds, 'Faxa_tauy')
    call addfld(fldListTo(compocn)%flds, 'Foxx_taux')
    call addfld(fldListTo(compocn)%flds, 'Foxx_tauy')
    call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset')
    call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset')
    call addmap(fldListFr(compatm)%flds, 'Faxa_tauy', compocn, mapnstod_consf, 'one'  )
    call addmap(fldListFr(compatm)%flds, 'Faxa_taux', compocn, mapnstod_consf, 'one'  )

    ! to ocn: water flux due to melting ice from ice
    call addfld(fldListFr(compice)%flds, 'Fioi_meltw')
    call addfld(fldListTo(compocn)%flds, 'Fioi_meltw')
    call addmap(fldListFr(compice)%flds, 'Fioi_meltw', compocn,  mapfcopy, 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Fioi_meltw', &
         mrg_from1=compice, mrg_fld1='Fioi_meltw', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')

    ! to ocn: heat flux from melting ice from ice
    call addfld(fldListFr(compice)%flds, 'Fioi_melth')
    call addfld(fldListTo(compocn)%flds, 'Fioi_melth')
    call addmap(fldListFr(compice)%flds, 'Fioi_melth', compocn,  mapfcopy, 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Fioi_melth', &
         mrg_from1=compice, mrg_fld1='Fioi_melth', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')

    ! to ocn: salt flux from ice
    call addfld(fldListFr(compice)%flds, 'Fioi_salt')
    call addfld(fldListTo(compocn)%flds, 'Fioi_salt')
    call addmap(fldListFr(compice)%flds, 'Fioi_salt', compocn,  mapfcopy, 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Fioi_salt', &
         mrg_from1=compice, mrg_fld1='Fioi_salt', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! to ice: downward longwave heat flux from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
    call addfld(fldListTo(compice)%flds, 'Faxa_lwdn')
    call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', mrg_type1='copy')

    ! to ice: downward direct near-infrared incident solar radiation  from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
    call addfld(fldListTo(compice)%flds, 'Faxa_swndr')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', mrg_type1='copy')

    ! to ice: downward direct visible incident solar radiation from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
    call addfld(fldListTo(compice)%flds, 'Faxa_swvdr')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', mrg_type1='copy')

    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')
    call addfld(fldListTo(compice)%flds, 'Faxa_swndf')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', mrg_type1='copy')

    ! to ice: downward Diffuse visible incident solar radiation  from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
    call addfld(fldListTo(compice)%flds, 'Faxa_swvdf')
    call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', mrg_type1='copy')

    ! to ice: rain rate from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_rain')
    call addfld(fldListTo(compice)%flds, 'Faxa_rain')
    call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_rain', mrg_from1=compatm, mrg_fld1='Faxa_rain', mrg_type1='copy')

    ! to ice: snow rate from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_snow')
    call addfld(fldListTo(compice)%flds, 'Faxa_snow')
    call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compice, mapnstod_consf, 'one')
    call addmrg(fldListTo(compice)%flds, 'Faxa_snow', mrg_from1=compatm, mrg_fld1='Faxa_snow', mrg_type1='copy')

    ! to ice: zonal wind at the lowest model level from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_u')
    call addfld(fldListTo(compice)%flds, 'Sa_u')
    call addmap(fldListFr(compatm)%flds, 'Sa_u', compice, mapnstod_consf, 'one') ! patch in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')

    ! to ice: meridional wind at the lowest model level from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_v')
    call addfld(fldListTo(compice)%flds, 'Sa_v')
    call addmap(fldListFr(compatm)%flds, 'Sa_v', compice, mapnstod_consf, 'one') ! patch in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')

    ! to ice: height at the lowest model level from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_z')
    call addfld(fldListTo(compice)%flds, 'Sa_z')
    call addmap(fldListFr(compatm)%flds, 'Sa_z', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_z', mrg_from1=compatm, mrg_fld1='Sa_z', mrg_type1='copy')

    ! to ice: pressure at the lowest model level fromatm
    call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
    call addfld(fldListTo(compice)%flds, 'Sa_pbot')
    call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_pbot', mrg_from1=compatm, mrg_fld1='Sa_pbot', mrg_type1='copy')

    ! to ice: temperature at the lowest model level from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
    call addfld(fldListTo(compice)%flds, 'Sa_tbot')
    call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')

    ! to ice: potential temperature at the lowest model level from atm
    ! this is a derived quantity in med_phases_prep_ice
    call addfld(fldListFr(compatm)%flds, 'Sa_ptem')
    call addfld(fldListTo(compice)%flds, 'Sa_ptem')
    call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_ptem', mrg_from1=compatm, mrg_fld1='Sa_ptem', mrg_type1='copy')

    ! to ice: density at the lowest model level from atm
    ! this is a derived quantity in med_phases_prep_ice
    call addfld(fldListFr(compatm)%flds, 'Sa_dens')
    call addfld(fldListTo(compice)%flds, 'Sa_dens')
    call addmap(fldListFr(compatm)%flds, 'Sa_dens', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_dens', mrg_from1=compatm, mrg_fld1='Sa_dens', mrg_type1='copy')

    ! to ice: specific humidity at the lowest model level from atm
    call addfld(fldListFr(compatm)%flds, 'Sa_shum')
    call addfld(fldListTo(compice)%flds, 'Sa_shum')
    call addmap(fldListFr(compatm)%flds, 'Sa_shum', compice, mapnstod_consf, 'one') ! bilinear in nems
    call addmrg(fldListTo(compice)%flds, 'Sa_shum', mrg_from1=compatm, mrg_fld1='Sa_shum', mrg_type1='copy')

    ! to ice: sea surface temperature from ocn
    call addfld(fldListFr(compocn)%flds, 'So_t')
    call addfld(fldListTo(compice)%flds, 'So_t')
    call addmap(fldListFr(compocn)%flds, 'So_t', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

    ! to ice: sea surface salinity from ocn
    call addfld(fldListFr(compocn)%flds, 'So_s')
    call addfld(fldListTo(compice)%flds, 'So_s')
    call addmap(fldListFr(compocn)%flds, 'So_s', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_s', mrg_from1=compocn, mrg_fld1='So_s', mrg_type1='copy')

    ! to ice: zonal sea water velocity from ocn
    call addfld(fldListFr(compocn)%flds, 'So_u')
    call addfld(fldListTo(compice)%flds, 'So_u')
    call addmap(fldListFr(compocn)%flds, 'So_u', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_u', mrg_from1=compocn, mrg_fld1='So_u', mrg_type1='copy')

    ! to ice: meridional sea water velocity from ocn
    call addfld(fldListFr(compocn)%flds, 'So_v')
    call addfld(fldListTo(compice)%flds, 'So_v')
    call addmap(fldListFr(compocn)%flds, 'So_v', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_v', mrg_from1=compocn, mrg_fld1='So_v', mrg_type1='copy')

    ! to ice: zonal sea surface slope from ocn
    call addfld(fldListFr(compocn)%flds, 'So_dhdx')
    call addfld(fldListTo(compice)%flds, 'So_dhdx')
    call addmap(fldListFr(compocn)%flds, 'So_dhdx', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_dhdx', mrg_from1=compocn, mrg_fld1='So_dhdx', mrg_type1='copy')

    ! to ice: meridional sea surface slope from ocn
    call addfld(fldListFr(compocn)%flds, 'So_dhdy')
    call addfld(fldListTo(compice)%flds, 'So_dhdy')
    call addmap(fldListFr(compocn)%flds, 'So_dhdy', compice, mapfcopy , 'unset')
    call addmrg(fldListTo(compice)%flds, 'So_dhdy', mrg_from1=compocn, mrg_fld1='So_dhdy', mrg_type1='copy')

    ! to ice: ocean melt and freeze potential from ocn
    call addfld(fldListFr(compocn)%flds, 'Fioo_q')
    call addfld(fldListTo(compice)%flds, 'Fioo_q')
    call addmap(fldListFr(compocn)%flds, 'Fioo_q', compice,  mapfcopy, 'unset')
    call addmrg(fldListTo(compice)%flds, 'Fioo_q', mrg_from1=compocn, mrg_fld1='Fioo_q', mrg_type1='copy')

  end subroutine esmFldsExchange_nems_orig_active

end module esmFldsExchange_nems_orig_active_mod
