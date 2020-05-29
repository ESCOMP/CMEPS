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
    use esmflds               , only : compmed, compatm, complnd, compocn
    use esmflds               , only : compice, comprof, compwav, compglc, ncomps
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use esmflds               , only : mapuv_with_cart3d
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
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
    character(len=5)    :: iso(2)
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CX)   :: atm2ice_fmap='unset', atm2ice_smap='unset', atm2ice_vmap='unset'
    character(len=CX)   :: atm2ocn_fmap='unset', atm2ocn_smap='unset', atm2ocn_vmap='unset'
    character(len=CX)   :: atm2lnd_fmap='unset', atm2lnd_smap='unset'
    character(len=CX)   :: glc2lnd_smap='unset', glc2lnd_fmap='unset'
    character(len=CX)   :: glc2ice_rmap='unset'
    character(len=CX)   :: glc2ocn_liq_rmap='unset', glc2ocn_ice_rmap='unset'
    character(len=CX)   :: ice2atm_fmap='unset', ice2atm_smap='unset'
    character(len=CX)   :: ocn2atm_fmap='unset', ocn2atm_smap='unset'
    character(len=CX)   :: lnd2atm_fmap='unset', lnd2atm_smap='unset'
    character(len=CX)   :: lnd2glc_fmap='unset', lnd2glc_smap='unset'
    character(len=CX)   :: lnd2rof_fmap='unset'
    character(len=CX)   :: rof2lnd_fmap='unset'
    character(len=CX)   :: rof2ocn_fmap='unset', rof2ocn_ice_rmap='unset', rof2ocn_liq_rmap='unset'
    character(len=CX)   :: atm2wav_smap='unset', ice2wav_smap='unset', ocn2wav_smap='unset'
    character(len=CX)   :: wav2ocn_smap='unset'
    logical             :: flds_co2a  ! use case
    logical             :: flds_co2b  ! use case
    logical             :: flds_co2c  ! use case
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !iso(1) = ''
    !iso(2) = '_wiso'

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
      flds = (/'Sa_u', 'Sa_v', 'Sa_z', 'Sa_tbot', 'Sa_pbot', 'Sa_shum'/)
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
       allocate(flds(13))
       flds = (/'So_tref  ', 'So_qref  ', 'So_u10   ', 'So_ustar ', 'So_ssq   ', &
                'So_re    ', 'So_duu10n', 'Faox_lwup', 'Faox_sen ', 'Faox_lat ', &
                'Faox_evap', 'Faox_taux', 'Faox_tauy'/)
       do n = 1, size(flds)
          fldname = trim(flds(n))
          call addfld(fldListMed_aoflux%flds, trim(fldname))
       end do
       deallocate(flds)
    end if

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

       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum_wiso', rc=rc)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_shum_wiso')
          call addmap(fldListFr(compatm)%flds, 'Sa_shum_wiso', compocn, mapbilnr, 'one', atm2ocn_smap)
       end if

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
    ! TODO (mvertens, 2019-01-11): budget implemention needs to be done in CMEPS
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Fall_swnet')
       call addfld(fldListFr(compice)%flds, 'Faii_swnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compice, mapconsf, 'one'  , atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
       end if
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_swnet', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Faii_swnet', compocn, mapfcopy, 'unset', 'unset')
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
    allocate(suffix(4))
    suffix = (/'avsdr', 'avsdf', 'anidr', 'anidf'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Sl_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds, 'Si_'//trim(suffix(n)))
          call addfld(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)))
       else
          ! (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Si_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_ocnalb_a        , 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_smap)
             call addmap(fldListFr(compice)%flds, 'Si_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_smap)
             call addmap(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_smap)
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Sl_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBexp(compatm), 'Sx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_smap)
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='So_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged reference temperature at 2 meters
    ! to atm: merged reference specific humidity at 2 meters
    ! ---------------------------------------------------------------------
    allocate(suffix(3))
    suffix = (/'tref     ', 'u10      ', 'qref     '/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds , 'Sl_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds , 'Si_'//trim(suffix(n)))
          call addfld(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds , 'Sx_'//trim(suffix(n)))
       else
          ! (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Si_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmap(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)), compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds , 'Sx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Sl_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='So_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged zonal surface stress
    ! to atm: merged meridional surface stress
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)))
          call addfld(fldListFr(complnd)%flds, 'Fall_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds, 'Faii_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n)))
       else
          ! (non aqua-planet)
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap)
             call addmap(fldListFr(complnd)%flds , 'Fall_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Faii_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Fall_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Faii_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='Faox_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBexp(compatm), 'Faxx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='Faox_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

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
       ! merged ocn/ice/lnd temp and unmerged ocn temp
       if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_t', compatm, mapconsf , 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=complnd, mrg_fld1='Sl_t', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_t', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compocn, mrg_fld3='So_t', mrg_type3='merge', mrg_fracname3='ofrac')
          call addmrg(fldListTo(compatm)%flds, 'So_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

       ! aqua-planet - merged and unmerged ocn temp are the same
       else if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf, 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='merge', mrg_fracname1='ofrac')
          call addmrg(fldListTo(compatm)%flds, 'So_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface saturation specific humidity in ocean from med aoflux
    ! to atm: square of exch. coeff (tracers)               from med aoflux
    ! to atm: surface fraction velocity                     from med aoflux
    ! ---------------------------------------------------------------------
    allocate(flds(3))
    flds = (/'So_ssq', 'So_ustar', 'So_re'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds  , trim(fldname))
          call addfld(fldListTo(compatm)%flds , trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm) , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o , trim(fldname), rc=rc)) then
             call addmap(fldListMed_aoflux%flds    , trim(fldname), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds   , trim(fldname), &
                  mrg_from1=compmed, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: req. fields to satisfy mediator (can be removed later)
    !----------------------------------------------------------
    allocate(flds(2))
    flds = (/'Faxa_snowc', 'Faxa_snowl'/)
!
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
    flds = (/'Sa_topo', 'Sa_z', 'Sa_ptem', 'Sa_pbot'/)

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
    allocate(flds(6))
    flds = (/'Faxa_lwdn', 'Faxa_swdn', 'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swvdf'/)

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
    ! to ocn: surface upward longwave heat flux from mediator
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup')
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwup')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lwup', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwup', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac')
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
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compatm, mrg_fld2='Faxa_lwdn', mrg_type2='merge', mrg_fracname2='ofrac')
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
          call addmap(fldListFr(compatm)%flds, 'Faxa_swdn', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swdn', &
               mrg_from1=compatm, mrg_fld1='Faxa_swdn', mrg_type1='copy')
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
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_fmap)
       end if
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
               mrg_from1=compmed, mrg_fld1='Faox_sen', mrg_type1='merge', mrg_fracname1='ofrac')
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
               mrg_from1=compmed, mrg_fld1='Faox_lat', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_evap', &
               mrg_from1=compmed, mrg_fld1='Faox_evap', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds , 'Faox_lat_wiso' )
       call addfld(fldListTo(compocn)%flds, 'Foxx_lat_wiso' )
    else
       if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat_wiso', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat_wiso', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lat_wiso', &
               mrg_from1=compmed, mrg_fld1='Faox_lat_wiso', mrg_type1='merge', mrg_fracname1='ofrac')
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
          call addmap(fldListMed_aoflux%flds , 'So_duu10n', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compocn)%flds, 'So_duu10n', &
               mrg_from1=compmed, mrg_fld1='So_duu10n', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! to ocn: zonal wind at the lowest model level from atm
    ! to ocn: meridional wind at the lowest model level from atm
    ! to ocn: temperature at the lowest model level from atm
    ! to ocn: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Sa_pslv', 'Sa_u', 'Sa_v', 'Sa_tbot', 'Sa_shum'/)

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
    ! to ocn: merge zonal surface stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds , 'Fioi_'//trim(suffix(n)))
          call addfld(fldListFr(compatm)%flds , 'Faxa_'//trim(suffix(n)))
          call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(suffix(n)))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_'//trim(suffix(n)), compocn, mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='Faox_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac', &
                  mrg_from2=compice, mrg_fld2='Fioi_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac')
          end if
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
