module esmFldsExchange_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange(gcomp, phase, rc)

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
    use esmflds               , only : coupling_mode, mapuv_with_cart3d
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use med_internalstate_mod , only : mastertask, logunit

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    integer             :: i, n
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=64), allocatable :: flds(:)
    character(len=64), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    if (phase /= 'advertise') then
       return
    end if

    ! Determine supported coupling model
    coupling_mode = 'nems_orig'

    ! Initialize if use 3d cartesian mapping for u,v
    mapuv_with_cart3d = .false.

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
    ! FIELDS TO MEDIATOR component
    !=====================================================================

    ! masks from components
    call addfld(fldListFr(compocn)%flds, 'So_omask')
    call addfld(fldListFr(compice)%flds, 'Si_imask')
    call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')

    ! atm and ocn fields required for atm/ocn flux calculation'
    call addfld(fldListFr(compatm)%flds, 'Sa_u   ')
    call addfld(fldListFr(compatm)%flds, 'Sa_v'   )
    call addfld(fldListFr(compatm)%flds, 'Sa_z'   )
    call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
    call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
    call addfld(fldListFr(compatm)%flds, 'Sa_shum')
    call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, mapnstod_consf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, mapnstod_consf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, mapnstod_consf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapnstod_consf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, mapnstod_consf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Sa_shum', compocn, mapnstod_consf, 'none', 'unset')

    ! unused fields needed by the atm/ocn flux computation
    call addfld(fldListMed_aoflux%flds , 'So_tref'   )
    call addfld(fldListMed_aoflux%flds , 'So_qref'   )
    call addfld(fldListMed_aoflux%flds , 'So_u10'    )
    call addfld(fldListMed_aoflux%flds , 'So_ustar'  )
    call addfld(fldListMed_aoflux%flds , 'So_ssq'    )
    call addfld(fldListMed_aoflux%flds , 'So_re'     )
    call addfld(fldListMed_aoflux%flds , 'So_duu10n' )
    call addfld(fldListMed_aoflux%flds , 'Faox_lwup' )
    call addfld(fldListMed_aoflux%flds , 'Faox_sen'  )
    call addfld(fldListMed_aoflux%flds , 'Faox_lat'  )
    call addfld(fldListMed_aoflux%flds , 'Faox_evap' )
    call addfld(fldListMed_aoflux%flds , 'Faox_taux' )
    call addfld(fldListMed_aoflux%flds , 'Faox_tauy' )

    ! unused fields from ice - but that are needed to be realized by the cice cap
    call addfld(fldListFr(compice)%flds, 'Si_avsdf')
    call addfld(fldListFr(compice)%flds, 'Si_avsdr')
    call addfld(fldListFr(compice)%flds, 'Si_anidf')
    call addfld(fldListFr(compice)%flds, 'Si_anidr')
    call addfld(fldListFr(compice)%flds, 'mean_sw_pen_to_ocn')
    call addfld(fldListFr(compice)%flds, 'Faii_evap')

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    ! TODO: what about land mask?

    ! to atm: Fractions (computed in med_phases_prep_atm)
    call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
    call addfld(fldListTo(compatm)%flds, 'So_ofrac')

    ! to atm: unmerged from ice
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux, 
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux 
    ! - evaporation water flux from water
    allocate(suffix(6))
    suffix = (/'taux     ', 'tauy     ', 'lat      ', 'sen      ', 'lwup     ', 'evap     '/)
    do n = 1,size(suffix)
       call addfld(fldListFr(compice)%flds, 'Faii_'//trim(suffix(n)))
      !call addfld(fldListTo(compatm)%flds, 'Faii_'//trim(suffix(n))) ! TODO: add this in FV3 renaming
       call addfld(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n))) ! TODO: rmeove this in FV3 renameing

       call addmap(fldListFr(compice)%flds , 'Faii_'//trim(suffix(n)), compatm, mapnstod_consf, 'ifrac', 'unset')
       call addmrg(fldListTo(compatm)%flds , 'Faxx_'//trim(suffix(n)), &
            mrg_from1=compice, mrg_fld1='Faii_'//trim(suffix(n)), mrg_type1='copy')
      ! call addmrg(fldListTo(compatm)%flds , 'Faii_'//trim(suffix(n)), &
      !     mrg_from1=compice, mrg_fld1='Faii_'//trim(suffix(n)), mrg_type1='copy')
    end do
    deallocate(suffix)

    ! to atm: unmerged surface temperatures from ice and ocn
    call addfld(fldListFr(compice)%flds, 'Si_t')
    call addfld(fldListTo(compatm)%flds, 'Si_t')
    call addfld(fldListFr(compocn)%flds, 'So_t')
    call addfld(fldListTo(compatm)%flds, 'So_t')
    call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapnstod_consf, 'ifrac', 'unset')
    call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapnstod_consf, 'none' , 'unset')
    call addmrg(fldListTo(compatm)%flds, 'Si_t', mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='copy')
    call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

    ! to atm: surface snow depth             from ice unmerged
    ! to atm: mean ice volume per unit area  from ice unmerged
    ! to atm: mean snow volume per unit area from ice unmerged
    allocate(flds(2))
    flds = (/'Si_vice ', 'Si_vsno '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
          call addmap(fldListFr(compice)%flds, trim(fldname), compatm, mapnstod_consf, 'ifrac', 'unset')
          call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from1=compice, mrg_fld1=trim(fldname), mrg_type1='copy')
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! to ocn: fractional ice coverage wrt ocean from ice
    call addfld(fldListFr(compice)%flds, 'Si_ifrac')
    call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
    call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')

    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    allocate(flds(5))
    flds = (/'Faxa_lwdn ', 'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swndf'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       if (trim(coupling_mode) == 'nems_orig') then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapnstod_consf, 'none', 'unset')
       else
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'none', 'unset')
       end if
       call addmrg(fldListTo(compocn)%flds, trim(fldname), &
            mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end do
    deallocate(flds)

    ! to ocn: net shortwave radiation from med (custom merge in med_phases_prep_ocn)
    call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdr')
    call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdf')
    call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idr')
    call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idf')
    call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idf')
    call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdr')
    call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdf')
    call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idr')
    call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdr' , compocn, mapfcopy, 'unset', 'unset')
    call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdf' , compocn, mapfcopy, 'unset', 'unset')
    call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idr' , compocn, mapfcopy, 'unset', 'unset')
    call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idf' , compocn, mapfcopy, 'unset', 'unset')

    ! to ocn: merged longwave net heat flux
    call addfld(fldListTo(compocn)%flds , 'Foxx_lwnet')
    call addfld(fldListFr(compatm)%flds , 'Faxa_lwdn')
    call addfld(fldListFr(compatm)%flds , 'Faxa_lwnet')
    if (trim(coupling_mode) == 'nems_orig') then
       call addmap(fldListFr(compatm)%flds , 'Faxa_lwdn' , compocn, mapnstod_consf, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds , 'Faxa_lwnet', compocn, mapnstod_consf, 'none', 'unset')
    else
       call addmap(fldListFr(compatm)%flds , 'Faxa_lwdn' , compocn, mapconsf, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds , 'Faxa_lwnet', compocn, mapconsf, 'none', 'unset')
    end if

    !  to ocn: precipitation rate water equivalent and snow rate water equivalent from atm
    call addfld(fldListTo(compocn)%flds, 'Faxa_rain' )
    call addfld(fldListTo(compocn)%flds, 'Faxa_snow' )
    call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
    call addfld(fldListFr(compatm)%flds, 'Faxa_snow' )
    if (trim(coupling_mode) == 'nems_orig') then
       call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compocn, mapnstod_consf, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compocn, mapnstod_consf, 'none', 'unset')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compocn, mapconsf, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compocn, mapconsf, 'none', 'unset')
    end if
    call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', mrg_from1=compatm, mrg_fld1='Faxa_rain', mrg_type1='copy')
    call addmrg(fldListTo(compocn)%flds, 'Faxa_snow', mrg_from1=compatm, mrg_fld1='Faxa_snow', mrg_type1='copy')

    ! to ocn: merged sensible heat flux (custom merge in med_phases_prep_ocn)
    call addfld(fldListTo(compocn)%flds, 'Foxx_sen')
    call addfld(fldListFr(compatm)%flds, 'Faxa_sen')
    call addmap(fldListFr(compatm)%flds, 'Faxa_sen'  , compocn, mapconsf, 'none'  , 'unset')

    ! to ocn: surface latent heat flux and evaporation water flux (custom merge in med_phases_prep_ocn)
    call addfld(fldListTo(compocn)%flds, 'Foxx_evap')
    call addfld(fldListFr(compatm)%flds, 'Faxa_lat' )
    call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, mapconsf, 'none', 'unset')

    ! to ocn: sea level pressure from atm
    call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
    call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
    if (trim(coupling_mode) == 'nems_orig') then
       call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, mapnstod_consf, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compice, mapnstod_consf, 'none', 'unset')
    else 
       call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, mapbilnr, 'none', 'unset')
       call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compice, mapbilnr, 'none', 'unset')
    end if
    call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', mrg_from1=compatm, mrg_fld1='Sa_pslv', mrg_type1='copy')

    ! to ocn: merge zonal surface stress from ice and atm (custom merge calculation in med_phases_prep_ocn)
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)
    do n = 1,size(suffix)
       call addfld(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)))
       call addfld(fldListFr(compice)%flds, 'Fioi_'//trim(suffix(n)))
       call addfld(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)))
       if (trim(coupling_mode) == 'nems_orig') then
          call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)), compocn, mapnstod_consf, 'none', 'unset')
       else
          call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)), compocn, mapconsf, 'none', 'unset')
       end if
       call addmap(fldListFr(compice)%flds, 'Fioi_'//trim(suffix(n)), compocn, mapfcopy, 'unset', 'unset')
    end do
    deallocate(suffix)

    ! to ocn: water flux due to melting ice from ice
    call addfld(fldListFr(compice)%flds, 'Fioi_meltw')
    call addfld(fldListTo(compocn)%flds, 'Fioi_meltw')
    call addmap(fldListFr(compice)%flds, 'Fioi_meltw',    compocn,  mapfcopy, 'unset', 'unset')
    call addmrg(fldListTo(compocn)%flds, 'Fioi_meltw', &
         mrg_from1=compice, mrg_fld1='Fioi_meltw', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')

    ! to ocn: heat flux from melting ice from ice
    ! to ocn: salt flux from ice
    allocate(flds(2))
    flds = (/'Fioi_melth ', 'Fioi_salt  '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compice)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
       call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
       call addmrg(fldListTo(compocn)%flds, trim(fldname), &
            mrg_from1=compice, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! to ice: downward longwave heat flux from atm
    ! to ice: downward direct near-infrared incident solar radiation  from atm
    ! to ice: downward direct visible incident solar radiation        from atm
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! to ice: downward Diffuse visible incident solar radiation       from atm
    allocate(flds(5))
    flds = (/'Faxa_lwdn  '    , 'Faxa_swndr '   , 'Faxa_swvdr '   , 'Faxa_swndf ' , 'Faxa_swvdf '/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapconsf, 'none', 'unset')
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
    end do
    deallocate(flds)

    ! to ice: convective and large scale precipitation rate water equivalent from atm
    ! to ice: rain and snow rate from atm
    call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
    call addfld(fldListFr(compatm)%flds, 'Faxa_snow' )
    call addfld(fldListTo(compice)%flds, 'Faxa_rain' ) 
    call addfld(fldListTo(compice)%flds, 'Faxa_snow' )
    call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compice, mapconsf, 'none', 'unset')
    call addmap(fldListFr(compatm)%flds, 'Faxa_snow', compice, mapconsf, 'none', 'unset')
    call addmrg(fldListTo(compice)%flds, 'Faxa_rain', mrg_from1=compatm, mrg_fld1='Faxa_rain', mrg_type1='copy')
    call addmrg(fldListTo(compice)%flds, 'Faxa_snow', mrg_from1=compatm, mrg_fld1='Faxa_snow', mrg_type1='copy')

    ! to ice: height at the lowest model level from atm
    ! to ice: pressure at the lowest model level from atm
    ! to ice: temperature at the lowest model level from atm
    ! to ice: potential temperature at the lowest model level from atm
    ! to ice: density at the lowest model level from atm
    ! to ice: zonal wind at the lowest model level from atm
    ! to ice: meridional wind at the lowest model level from atm
    ! to ice: specific humidity at the lowest model level from atm
    allocate(flds(7))
    flds = (/'Sa_z        ', 'Sa_pbot     ', 'Sa_tbot     ', 'Sa_dens     ', 'Sa_u        ', 'Sa_v        ', 'Sa_shum     '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       if (trim(coupling_mode) == 'nems_orig') then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapnstod_consf, 'none', 'unset')
       else
          if (trim(fldname) == 'Sa_u' .or. trim(fldname) == 'Sa_v') then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mappatch, 'none', 'unset')
          else
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapbilnr, 'none', 'unset')
          end if
       end if
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
    end do
    deallocate(flds)

    ! to ice: sea surface temperature from ocn
    ! to ice: sea surface salinity from ocn
    ! to ice: zonal sea water velocity from ocn
    ! to ice: meridional sea water velocity from ocn
    ! to ice: zonal sea surface slope from ocean
    ! to ice: meridional sea surface slope from ocn
    allocate(flds(6))
    flds = (/'So_t   ', 'So_s   ', 'So_u   ', 'So_v   ', 'So_dhdx', 'So_dhdy'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       call addfld(fldListTo(compice)%flds, trim(fldname))
       call addfld(fldListFr(compocn)%flds, trim(fldname))
       call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy , 'unset', 'unset')
       call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from1=compocn, mrg_fld1=trim(fldname), mrg_type1='copy')
    end do
    deallocate(flds)

    ! to ice: ocean melt and freeze potential from ocn
    call addfld(fldListTo(compice)%flds, 'Fioo_q')
    call addfld(fldListFr(compocn)%flds, 'Fioo_q')
    call addmap(fldListFr(compocn)%flds, 'Fioo_q', compice,  mapfcopy, 'unset', 'unset')
    call addmrg(fldListTo(compice)%flds, 'Fioo_q', mrg_from1=compocn, mrg_fld1='Fioo_q', mrg_type1='copy')

  end subroutine esmFldsExchange

end module esmFldsExchange_mod
