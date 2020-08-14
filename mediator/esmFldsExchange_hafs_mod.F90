module esmFldsExchange_hafs_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_hafs

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange_hafs(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : compmed, compatm, compocn, compice, ncomps
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use esmflds               , only : mapuv_with_cart3d
    use esmflds               , only : fldListTo, fldListFr
    use esmFlds               , only : coupling_mode

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: dbrc
    integer             :: num, i, n
    integer             :: n1, n2, n3, n4
    logical             :: isPresent
    !character(len=5)    :: iso(2)
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CX)   :: atm2ice_fmap='unset', atm2ice_smap='unset', atm2ice_vmap='unset'
    character(len=CX)   :: atm2ocn_fmap='unset', atm2ocn_smap='unset', atm2ocn_vmap='unset'
    character(len=CX)   :: ice2atm_fmap='unset', ice2atm_smap='unset'
    character(len=CX)   :: ocn2atm_fmap='unset', ocn2atm_smap='unset'
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    if (phase /= 'advertise') then
       nullify(is_local%wrap)
       call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------------
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

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmapname', value=ice2atm_fmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('ice2atm_fmapname = '// trim(ice2atm_fmap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smapname', value=ice2atm_smap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('ice2atm_smapname = '// trim(ice2atm_smap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', value=ocn2atm_smap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('ocn2atm_smapname = '// trim(ocn2atm_smap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', value=ocn2atm_fmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('ocn2atm_fmapname = '// trim(ocn2atm_fmap), ESMF_LOGMSG_INFO)
    end if

    ! to ice

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmapname', value=atm2ice_fmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ice_fmapname = '// trim(atm2ice_fmap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', value=atm2ice_smap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ice_smapname = '// trim(atm2ice_smap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmapname', value=atm2ice_vmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ice_vmapname = '// trim(atm2ice_vmap), ESMF_LOGMSG_INFO)
    end if

    ! to ocn

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', value=atm2ocn_fmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ocn_fmapname = '// trim(atm2ocn_fmap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', value=atm2ocn_smap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ocn_smapname = '// trim(atm2ocn_smap), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', value=atm2ocn_vmap, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call ESMF_LogWrite('atm2ocn_vmapname = '// trim(atm2ocn_vmap), ESMF_LOGMSG_INFO)
    end if

    !----------------------------------------------------------
    ! Initialize if use 3d cartesian mapping for u,v
    !----------------------------------------------------------

    mapuv_with_cart3d = .false.

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
       call addfld(fldListFr(compocn)%flds, 'So_omask')
       call addfld(fldListFr(compice)%flds, 'Si_imask')
    else
       call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
    end if

    ! ---------------------------------------------------------------------
    ! to med: atm and ocn fields required for atm/ocn flux calculation
    ! ---------------------------------------------------------------------
    if (phase /= 'advertise') then
      allocate(flds(6))
      flds = (/'Sa_u   ', 'Sa_v   ', 'Sa_z   ', 'Sa_tbot', 'Sa_pbot', 'Sa_shum'/)
      do n = 1,size(flds)
         fldname = trim(flds(n))
         call addfld(fldListFr(compatm)%flds, trim(fldname))
         if (trim(fldname) == 'Sa_u' .or. trim(fldname) == 'Sa_v') then
            !call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mappatch, 'one', atm2ocn_vmap)
            call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', atm2ocn_smap)
         else
            call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', atm2ocn_smap)
         end if
      end do
      deallocate(flds)
    end if

    ! ---------------------------------------------------------------------
    ! to med: unused fields needed by the atm/ocn flux computation
    ! ---------------------------------------------------------------------

    if (phase /= 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_u')
       !call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, mappatch, 'one', atm2ocn_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, mapbilnr, 'one', atm2ocn_vmap)

       call addfld(fldListFr(compatm)%flds, 'Sa_v')
       !call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, mappatch, 'one', atm2ocn_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, mapbilnr, 'one', atm2ocn_vmap)

       call addfld(fldListFr(compatm)%flds, 'Sa_z')
       call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, mapbilnr, 'one', atm2ocn_smap)

       call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
       call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
       call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addfld(fldListFr(compatm)%flds, 'Sa_shum')
       call addmap(fldListFr(compatm)%flds, 'Sa_shum', compocn, mapbilnr, 'one', atm2ocn_smap)

       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_ptem')
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, mapbilnr, 'one', atm2ocn_smap)
       end if

       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_dens')
          call addmap(fldListFr(compatm)%flds, 'Sa_dens', compocn, mapbilnr, 'one', atm2ocn_smap)
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one', atm2ocn_fmap)
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    !----------------------------------------------------------
    ! to atm: Fractions
    !----------------------------------------------------------
    if (phase == 'advertise') then
       ! the following are computed in med_phases_prep_atm
       call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
       call addfld(fldListTo(compatm)%flds, 'So_ofrac')
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: req. fields to satisfy mediator (can be removed later)
    !----------------------------------------------------------
    allocate(flds(2))
    flds = (/'Faxa_snowc', 'Faxa_snowl'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    allocate(flds(4))
    flds = (/'Sa_topo', 'Sa_z   ', 'Sa_ptem', 'Sa_pbot'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', atm2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
    else
       call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')
       call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_lwdn ', 'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swvdf'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from1=compatm, mrg_fld1=trim(fldname), &
                  mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: longwave net heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwnet')
       call addfld(fldListTo(compocn)%flds, 'Foxx_lwnet')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapconsf, 'one', atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
            mrg_from1=compatm, mrg_fld1='Faxa_lwnet', mrg_type1='copy')
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
          call addmap(fldListFr(compatm)%flds, 'Faxa_swdn', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swdn', &
               mrg_from1=compatm, mrg_fld1='Faxa_swdn', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one', atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Foxx_swnet', &
            mrg_from1=compatm, mrg_fld1='Faxa_swnet', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rainl')
       call addfld(fldListFr(compatm)%flds, 'Faxa_rain' )
       call addfld(fldListTo(compocn)%flds, 'Faxa_rain' )
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain' , rc=rc)) then
           call addmap(fldListFr(compatm)%flds, 'Faxa_rainl', compocn, mapconsf, 'one', atm2ocn_fmap)
           call addmap(fldListFr(compatm)%flds, 'Faxa_rainc', compocn, mapconsf, 'one', atm2ocn_fmap)
           call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', &
                mrg_from1=compatm, mrg_fld1='Faxa_rainc:Faxa_rainl', &
                mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
       else if (fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain', rc=rc)) then
           call addmap(fldListFr(compatm)%flds, 'Faxa_rain', compocn, mapconsf, 'one', atm2ocn_fmap)
           call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', mrg_from1=compatm, mrg_fld1='Faxa_rain', &
                mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sensible heat flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_sen')
       call addfld(fldListTo(compocn)%flds , 'Foxx_sen')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapconsf, 'one', atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
            mrg_from1=compatm, mrg_fld1='Faxa_sen', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_lat')
       call addfld(fldListTo(compocn)%flds , 'Foxx_lat')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, mapconsf, 'one', atm2ocn_fmap)
       call addmrg(fldListTo(compocn)%flds, 'Foxx_lat', &
            mrg_from1=compatm, mrg_fld1='Faxa_lat', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! to ocn: zonal wind at the lowest model level from atm
    ! to ocn: meridional wind at the lowest model level from atm
    ! to ocn: wind speed at the lowest model level from atm
    ! to ocn: temperature at the lowest model level from atm
    ! to ocn: sea surface skin temperature
    ! to ocn: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(7))
    flds = (/'Sa_pslv', 'Sa_u   ', 'Sa_v   ', 'Sa_wspd', 'Sa_tbot', 'Sa_tskn', 'Sa_shum'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if (fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname), rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(compocn)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', atm2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: zonal and meridional surface stress from atm
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds , 'Faxa_'//trim(suffix(n)))
          call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(suffix(n)))
       else
          call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)), compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)), &
               mrg_from1=compatm, mrg_fld1='Faxa_'//trim(suffix(n)), mrg_type1='copy')
       end if
    end do
    deallocate(suffix)

    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: density at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(1))
    flds = (/'Sa_dens'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapbilnr, 'one', atm2ice_smap)
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: zonal sea water velocity from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(2))
    flds = (/'So_u   ', 'So_v   '/)

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
                  mrg_from1=compocn, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

  end subroutine esmFldsExchange_hafs

end module esmFldsExchange_hafs_mod
