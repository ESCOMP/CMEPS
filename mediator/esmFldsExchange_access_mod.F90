module esmFldsExchange_access_mod

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
    use med_internalstate_mod , only : ncomps
    use med_internalstate_mod , only : coupling_mode
    use esmflds               , only : fldListTo
    use esmflds               , only : fldListFr

    !---------------------------------------------------------------------
    ! This is a mediator specific routine that determines ALL possible
    ! fields exchanged between components and their associated routing,
    ! mapping and merging
    !---------------------------------------------------------------------

    implicit none
    public

    public :: esmFldsExchange_access

    character(*), parameter :: u_FILE_u = &
         __FILE__

  !===============================================================================
  contains
  !===============================================================================

    subroutine esmFldsExchange_access(gcomp, phase, rc)

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      character(len=*) , parameter   :: subname='(esmFldsExchange_access)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      if (phase == 'advertise') then
        call esmFldsExchange_access_advt(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (phase == 'fieldcheck') then
        call esmFldsExchange_access_fchk(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (phase == 'initialize') then
        call esmFldsExchange_access_init(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      else
        call ESMF_LogSetError(ESMF_FAILURE, &
           msg=trim(subname)//": Phase is set to "//trim(phase), &
           line=__LINE__, file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_advt(gcomp, phase, rc)

      use esmFlds, only : addfld => med_fldList_AddFld

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      integer             :: num, i, n
      logical             :: isPresent
      character(len=CL)   :: cvalue
      character(len=CS)   :: name, fldname
      character(len=CS)   :: fldname1, fldname2
      character(len=CS), allocatable :: flds(:)
      character(len=CS), allocatable :: S_flds(:)
      character(len=CS), allocatable :: F_flds(:,:)
      character(len=CS), allocatable :: suffix(:)
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_advt)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      !=====================================================================
      ! scalar information
      !=====================================================================

      call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', &
         isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
         call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", &
            value=cvalue, rc=rc)
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
      call addfld(fldListFr(compocn)%flds, 'So_omask')
      call addfld(fldListFr(compocn)%flds, 'So_imask')

      !=====================================================================
      ! FIELDS TO ATMOSPHERE
      !=====================================================================

      call addfld(fldListTo(compatm)%flds, 'So_ofrac')
      call addfld(fldListTo(compatm)%flds, 'So_ifrac')

      ! ---------------------------------------------------------------------
      ! to atm: from ocn
      ! ---------------------------------------------------------------------
      allocate(S_flds(1))
      S_flds = (/'So_t'/) ! sea_surface_temperature
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        call addfld(fldListFr(compocn)%flds, trim(fldname))
        call addfld(fldListTo(compatm)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to atm: from ice
      ! ---------------------------------------------------------------------
      allocate(S_flds(1))
      S_flds = (/'Si_t'/) ! sea_surface_temperature
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        call addfld(fldListFr(compice)%flds, trim(fldname))
        call addfld(fldListTo(compatm)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'inst_pres_height_surface', & ! inst_zonal_wind_height10m
                  'So_duu10n' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld(fldListFr(compatm)%flds, trim(fldname))
         call addfld(fldListTo(compocn)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(15, 2))
      F_flds(1,:) = (/'Faxa_taux ', 'Foxx_taux'/)
      F_flds(2,:) = (/'Faxa_tauy ', 'Foxx_tauy'/)
      F_flds(3,:) = (/'mean_sensi_heat_flx', 'mean_sensi_heat_flx'/)
      F_flds(4,:) = (/'mean_evap_rate', 'mean_evap_rate'/)
      F_flds(5,:) = (/'mean_net_lw_flx', 'mean_net_lw_flx'/)
      F_flds(6,:) = (/'mean_net_sw_vis_dir_flx', 'mean_net_sw_vis_dir_flx'/)
      F_flds(7,:) = (/'mean_net_sw_vis_dif_flx', 'mean_net_sw_vis_dif_flx'/)
      F_flds(8,:) = (/'mean_net_sw_ir_dir_flx', 'mean_net_sw_ir_dir_flx'/)
      F_flds(9,:) = (/'mean_net_sw_ir_dif_flx', 'mean_net_sw_ir_dif_flx'/)
      F_flds(10,:) = (/'Faxa_rainc', 'mean_prec_rate'/)
      F_flds(11,:) = (/'Faxa_rainl', 'mean_prec_rate'/)
      F_flds(12,:) = (/'Faxa_snowc', 'mean_fprec_rate'/)
      F_flds(13,:) = (/'Faxa_snowl', 'mean_fprec_rate'/)
      F_flds(14,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(15,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld(fldListFr(compatm)%flds, trim(fldname1))
         call addfld(fldListTo(compocn)%flds, trim(fldname2))
      end do
      deallocate(F_flds)

      ! from ice
      allocate(F_flds(6, 2))
      F_flds(1,:) = (/'mean_salt_rate', 'mean_salt_rate'/)
      F_flds(2,:) = (/'Si_ifrac', 'Si_ifrac'/) ! ice_fraction
      F_flds(3,:) = (/'mean_fresh_water_to_ocean_rate', 'mean_fresh_water_to_ocean_rate'/)
      F_flds(4,:) = (/'net_heat_flx_to_ocn', 'net_heat_flx_to_ocn'/) ! heat flux sea-ice to ocean
      F_flds(5,:) = (/'Fioi_taux', 'Foxx_taux'/)
      F_flds(6,:) = (/'Fioi_tauy', 'Foxx_tauy'/) ! heat flux sea-ice to ocean
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld(fldListFr(compice)%flds, trim(fldname1))
         call addfld(fldListTo(compocn)%flds, trim(fldname2))
      end do
      deallocate(F_flds)

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: state fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(S_flds(2))
      S_flds = (/'inst_height_lowest', &
                  'inst_zonal_wind_height_lowest', &
                  'inst_merid_wind_height_lowest', &
                  'inst_spec_humid_height_lowest', &
                  'inst_temp_height_lowest', &
                  'inst_pres_height_lowest', &
                  'air_density_height_lowest', &
                  'Sa_ptem'/)
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld(fldListFr(compatm)%flds, trim(fldname))
         call addfld(fldListTo(compice)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      ! from ocn
      allocate(S_flds(7))
      S_flds = (/'sea_surface_slope_zonal', & ! inst_zonal_wind_height10m
                  'sea_surface_slope_merid', & ! inst_merid_wind_height10m
                  'sea_surface_temperature ', & ! inst_temp_height2m
                  's_surf ', & ! inst_spec_humid_height2m
                  'ocn_current_zonal', & ! inst_pres_height_surface
                  'ocn_current_merid', & ! inst_pres_height_surface
                  'freezing_melting_potential' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld(fldListFr(compocn)%flds, trim(fldname))
         call addfld(fldListTo(compice)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ice: flux fields
      ! ---------------------------------------------------------------------
      allocate(F_flds(7, 2))
      F_flds(1,:) = (/'mean_down_sw_vis_dir_flx ', 'mean_down_sw_vis_dir_flx '/)
      F_flds(2,:) = (/'mean_down_sw_ir_dir_flx ', 'mean_down_sw_ir_dir_flx '/)
      F_flds(3,:) = (/'mean_down_sw_vis_dif_flx', 'mean_down_sw_vis_dif_flx'/)
      F_flds(4,:) = (/'mean_down_sw_ir_dif_flx', 'mean_down_sw_ir_dif_flx'/)
      F_flds(5,:) = (/'mean_down_lw_flx', 'mean_down_lw_flx'/)
      F_flds(6,:) = (/'mean_prec_rate', 'mean_prec_rate'/)
      F_flds(7,:) = (/'mean_fprec_rate', 'mean_fprec_rate'/)
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld(fldListFr(compatm)%flds, trim(fldname1))
         call addfld(fldListTo(compice)%flds, trim(fldname2))
      end do
      deallocate(F_flds)

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_advt

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_fchk(gcomp, phase, rc)

      use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
      use med_internalstate_mod , only : InternalState

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      type(InternalState) :: is_local
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_fchk)'
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

    end subroutine esmFldsExchange_access_fchk

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_init(gcomp, phase, rc)

      use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
      use med_internalstate_mod , only : InternalState
      use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
      use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd
      use med_internalstate_mod , only : mapfillv_bilnr
      use med_internalstate_mod , only : mapnstod_consf
      use esmFlds               , only : med_fldList_type
      use esmFlds               , only : addmap => med_fldList_AddMap
      use esmFlds               , only : addmrg => med_fldList_AddMrg

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      type(InternalState) :: is_local
      integer             :: num, i, n
      integer             :: n1, n2, n3, n4
      character(len=CL)   :: cvalue
      character(len=CS)   :: name, fldname
      character(len=CS)   :: fldname1, fldname2
      character(len=CS), allocatable :: flds(:)
      character(len=CS), allocatable :: S_flds(:)
      character(len=CS), allocatable :: F_flds(:,:)
      character(len=CS), allocatable :: suffix(:)
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_init)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      !---------------------------------------
      ! Get the internal state
      !---------------------------------------
      nullify(is_local%wrap)
      call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      !=====================================================================
      ! FIELDS TO ATMOSPHERE
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to atm: sea surface temperature
      ! ---------------------------------------------------------------------
      allocate(S_flds(1))
      S_flds = (/'So_t'/) ! sea_surface_temperature
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compatm), trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname),rc=rc) &
            ) then

            call addmap(fldListFr(compocn)%flds, trim(fldname), compatm, mapconsf, 'ofrac', 'unset')
            call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      allocate(S_flds(1))
      S_flds = (/'Si_t'/) ! sea ice temperature
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compatm), trim(fldname),rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname),rc=rc) &
            ) then

            call addmap(fldListFr(compice)%flds, trim(fldname), compatm, mapconsf, 'ifrac', 'unset')
            call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')

            end if
      end do
      deallocate(S_flds)

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'inst_pres_height_surface', & ! inst_zonal_wind_height10m
                  'So_duu10n' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname), rc=rc) &
            ) then

            call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', atm2ice_map)
            call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(11, 2))
      F_flds(1,:) = (/'mean_zonal_moment_flx', 'mean_zonal_moment_flx'/)
      F_flds(2,:) = (/'mean_merid_moment_flx', 'mean_merid_moment_flx'/)
      F_flds(3,:) = (/'mean_sensi_heat_flx', 'mean_sensi_heat_flx'/)
      F_flds(4,:) = (/'mean_evap_rate', 'mean_evap_rate'/)
      F_flds(5,:) = (/'mean_net_lw_flx', 'mean_net_lw_flx'/)
      F_flds(6,:) = (/'mean_net_sw_vis_dir_flx', 'mean_net_sw_vis_dir_flx'/)
      F_flds(7,:) = (/'mean_net_sw_vis_dif_flx', 'mean_net_sw_vis_dif_flx'/)
      F_flds(8,:) = (/'mean_net_sw_ir_dir_flx', 'mean_net_sw_ir_dir_flx'/)
      F_flds(9,:) = (/'mean_net_sw_ir_dif_flx', 'mean_net_sw_ir_dif_flx'/)
      F_flds(10,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(11,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname1), rc=rc) &
            ) then
            call addmap(fldListFr(compatm)%flds, trim(fldname1), compocn, mapconsf, 'one', 'unset')
            call addmrg(fldListTo(compocn)%flds, trim(fldname2), mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')
         end if
      end do
      deallocate(F_flds)

      ! precip
      call addmap_from(compatm, 'Faxa_rainc', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'mean_prec_rate', mrg_from=compatm, mrg_fld='Faxa_rainc', mrg_type='merge')
      call addmap_from(compatm, 'Faxa_rainl', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'mean_prec_rate', mrg_from=compatm, mrg_fld='Faxa_rainl', mrg_type='merge')

      call addmap_from(compatm, 'Faxa_snowc', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'mean_fprec_rate', mrg_from=compatm, mrg_fld='Faxa_snowc', mrg_type='merge')
      call addmap_from(compatm, 'Faxa_snowl', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'mean_fprec_rate', mrg_from=compatm, mrg_fld='Faxa_snowl', mrg_type='merge')

      ! from ice
      allocate(F_flds(4, 2))
      F_flds(1,:) = (/'mean_salt_rate', 'mean_salt_rate'/)
      F_flds(2,:) = (/'Si_ifrac', 'Si_ifrac'/) ! ice_fraction
      F_flds(3,:) = (/'mean_fresh_water_to_ocean_rate', 'mean_fresh_water_to_ocean_rate'/)
      F_flds(4,:) = (/'net_heat_flx_to_ocn', 'net_heat_flx_to_ocn'/) ! heat flux sea-ice to ocean
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname1),rc=rc) &
            ) then
            call addmap(fldListFr(compice)%flds, trim(fldname1), compocn, mapfcopy, 'unset', 'unset')
            call addmrg(fldListTo(compocn)%flds, trim(fldname2), mrg_from=compice, mrg_fld=trim(fldname1), mrg_type='copy')
         end if
      end do
      deallocate(F_flds)

      ! momentum transfer
      call addmap_from(compice, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
      call addmrg_to(compocn, 'Foxx_taux', mrg_from=compice, mrg_fld='Fioi_taux', mrg_type='merge', mrg_fracname='ifrac')
      call addmap_from(compatm, 'Faxa_taux', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Foxx_taux', mrg_from=compatm, mrg_fld='Faxa_taux', mrg_type='merge', mrg_fracname='ofrac')

      call addmap_from(compice, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
      call addmrg_to(compocn, 'Foxx_tauy', mrg_from=compice, mrg_fld='Fioi_tauy', mrg_type='merge', mrg_fracname='ifrac')
      call addmap_from(compatm, 'Faxa_tauy', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Foxx_tauy', mrg_from=compatm, mrg_fld='Faxa_tauyx', mrg_type='merge', mrg_fracname='ofrac')

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: state fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(S_flds(8))
      S_flds = (/'inst_height_lowest', & ! inst_zonal_wind_height10m
                  'inst_zonal_wind_height_lowest', & ! inst_merid_wind_height10m
                  'inst_merid_wind_height_lowest ', & ! inst_temp_height2m
                  'inst_spec_humid_height_lowest ', & ! inst_spec_humid_height2m
                  'inst_temp_height_lowest', & ! inst_pres_height_surface
                  'inst_pres_height_lowest', &
                  'air_density_height_lowest', &
                  'Sa_ptem' /) ! inst_temp_height_surface

      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice), trim(fldname),rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname),rc=rc) &
            ) then

            call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapbilnr, 'one', 'unset')
            call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! from ocn
      allocate(S_flds(6))
      S_flds = (/'sea_surface_slope_zonal', & ! inst_zonal_wind_height10m
                  'sea_surface_slope_merid', & ! inst_merid_wind_height10m
                  'sea_surface_temperature ', & ! inst_temp_height2m
                  's_surf ', & ! inst_spec_humid_height2m
                  'ocn_current_zonal', & ! inst_pres_height_surface
                  'ocn_current_merid', & ! inst_pres_height_surface
                  'freezing_melting_potential' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname),rc=rc) &
            ) then

            call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy, 'unset', 'unset')
            call addmrg(fldListTo(compice)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ice: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(5, 2))
      F_flds(1,:) = (/'mean_down_sw_vis_dir_flx ', 'mean_down_sw_vis_dir_flx '/)
      F_flds(2,:) = (/'mean_down_sw_ir_dir_flx ', 'mean_down_sw_ir_dir_flx '/)
      F_flds(3,:) = (/'mean_down_sw_vis_dif_flx', 'mean_down_sw_vis_dif_flx'/)
      F_flds(4,:) = (/'mean_down_sw_ir_dif_flx', 'mean_down_sw_ir_dif_flx'/)
      F_flds(5,:) = (/'mean_down_lw_flx', 'mean_down_lw_flx'/)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compice), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname1), rc=rc) &
            ) then

            call addmap(fldListFr(compatm)%flds, trim(fldname1), compice, mapconsf, 'one', 'unset')
            call addmrg(fldListTo(compice)%flds, trim(fldname2), mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')

         end if
      end do
      deallocate(F_flds)

      ! precip
      call addmap_from(compatm, 'Faxa_rainc', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'mean_prec_rate', mrg_from=compatm, mrg_fld='Faxa_rainc', mrg_type='merge')
      call addmap_from(compatm, 'Faxa_rainl', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'mean_prec_rate', mrg_from=compatm, mrg_fld='Faxa_rainl', mrg_type='merge')

      call addmap_from(compatm, 'Faxa_snowc', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'mean_fprec_rate', mrg_from=compatm, mrg_fld='Faxa_snowc', mrg_type='merge')
      call addmap_from(compatm, 'Faxa_snowl', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'mean_fprec_rate', mrg_from=compatm, mrg_fld='Faxa_snowl', mrg_type='merge')

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_init

    !-----------------------------------------------------------------------------

  end module esmFldsExchange_access_mod
