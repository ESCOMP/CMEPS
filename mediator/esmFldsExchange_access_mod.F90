module esmFldsExchange_access_mod

    use ESMF
    use NUOPC
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_kind_mod          , only : CX=>SHR_KIND_CX
    use med_kind_mod          , only : CS=>SHR_KIND_CS
    use med_kind_mod          , only : CL=>SHR_KIND_CL
    use med_kind_mod          , only : R8=>SHR_KIND_R8
    use med_internalstate_mod , only : compmed, compatm, compocn, compwav, compice
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
      call addfld(fldListFr(compice)%flds, 'Si_imask')

      !=====================================================================
      ! FIELDS TO ATMOSPHERE
      !=====================================================================

      call addfld(fldListTo(compatm)%flds, 'So_ofrac')
      call addfld(fldListTo(compatm)%flds, 'Si_ifrac')

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
      call addfld(fldListFr(compice)%flds, 'Si_t')
      call addfld(fldListTo(compatm)%flds, 'Si_t')

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'Sa_pslv', & ! inst_zonal_wind_height10m
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
      allocate(F_flds(13, 2))
      F_flds(1,:) = (/'Faxa_taux ', 'Foxx_taux'/)
      F_flds(2,:) = (/'Faxa_tauy ', 'Foxx_tauy'/)
      F_flds(3,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(4,:) = (/'Foxx_evap', 'Foxx_evap'/)
      F_flds(5,:) = (/'Foxx_lwnet', 'Foxx_lwnet'/)
      F_flds(6,:) = (/'Foxx_swnet_vdr', 'Foxx_swnet_vdr'/)
      F_flds(7,:) = (/'Foxx_swnet_vdf', 'Foxx_swnet_vdf'/)
      F_flds(8,:) = (/'Foxx_swnet_idr', 'Foxx_swnet_idr'/)
      F_flds(9,:) = (/'Foxx_swnet_idf', 'Foxx_swnet_idf'/)
      F_flds(10,:) = (/'Faxa_rainc', 'Faxa_rain'/)
      F_flds(11,:) = (/'Faxa_snowc', 'Faxa_snow'/)
      F_flds(12,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(13,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld(fldListFr(compatm)%flds, trim(fldname1))
         call addfld(fldListTo(compocn)%flds, trim(fldname2))
      end do
      deallocate(F_flds)

      call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
      call addfld(fldListFr(compatm)%flds, 'Faxa_snowc')

      ! from ice
      allocate(F_flds(6, 2))
      F_flds(1,:) = (/'Fioi_salt', 'Fioi_salt'/)
      F_flds(2,:) = (/'Si_ifrac', 'Si_ifrac'/) ! ice_fraction
      F_flds(3,:) = (/'Fioi_meltw', 'Fioi_meltw'/)
      F_flds(4,:) = (/'Fioi_melth', 'Fioi_melth'/) ! heat flux sea-ice to ocean
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
      S_flds = (/'Sa_z', &
                  'Sa_u', &
                  'Sa_v', &
                  'Sa_shum', &
                  'Sa_tbot', &
                  'Sa_pbot', &
                  'Sa_dens', &
                  'Sa_ptem'/)
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld(fldListFr(compatm)%flds, trim(fldname))
         call addfld(fldListTo(compice)%flds, trim(fldname))
      end do
      deallocate(S_flds)

      ! from ocn
      allocate(S_flds(7))
      S_flds = (/'So_dhdx', & 
                 'So_dhdy', &
                 'So_t', & 
                 'So_s', & 
                 'So_u', & 
                 'So_v', & 
                 'Fioo_q' /) 
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
      F_flds(1,:) = (/'Faxa_swvdr ', 'Faxa_swvdr '/)
      F_flds(2,:) = (/'Faxa_swndr ', 'Faxa_swndr '/)
      F_flds(3,:) = (/'Faxa_swvdf', 'Faxa_swvdf'/)
      F_flds(4,:) = (/'Faxa_swndf', 'Faxa_swndf'/)
      F_flds(5,:) = (/'Faxa_lwdn', 'Faxa_lwdn'/)
      F_flds(6,:) = (/'Faxa_rainl', 'Faxa_rain'/)
      F_flds(7,:) = (/'Faxa_snowl', 'Faxa_snow'/)
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld(fldListFr(compatm)%flds, trim(fldname1))
         call addfld(fldListTo(compice)%flds, trim(fldname2))
      end do
      deallocate(F_flds)

      call addfld(fldListFr(compatm)%flds, 'Faxa_rainc')
      call addfld(fldListFr(compatm)%flds, 'Faxa_snowc')
      
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
      call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf, 'ofrac', 'unset')
      call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')

      call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf, 'ifrac', 'unset')
      call addmrg(fldListTo(compatm)%flds, 'Si_t', mrg_from=compice, mrg_fld='Si_t', mrg_type='copy')

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'Sa_pslv', & ! inst_zonal_wind_height10m
                  'So_duu10n' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname), rc=rc) &
            ) then

            call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapbilnr, 'one', 'unset')
            call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(11, 2))
      F_flds(1,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(2,:) = (/'Foxx_evap', 'Foxx_evap'/)
      F_flds(3,:) = (/'Foxx_lwnet', 'Foxx_lwnet'/)
      F_flds(4,:) = (/'Foxx_swnet_vdr', 'Foxx_swnet_vdr'/)
      F_flds(5,:) = (/'Foxx_swnet_vdf', 'Foxx_swnet_vdf'/)
      F_flds(6,:) = (/'Foxx_swnet_idr', 'Foxx_swnet_idr'/)
      F_flds(7,:) = (/'Foxx_swnet_idf', 'Foxx_swnet_idf'/)
      F_flds(8,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(9,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

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
      call addmap(fldListFr(compatm)%flds, 'Faxa_rainc', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', mrg_from=compatm, mrg_fld='Faxa_rainc', mrg_type='sum')
      call addmap(fldListFr(compatm)%flds, 'Faxa_rainl', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Faxa_rain', mrg_from=compatm, mrg_fld='Faxa_rainl', mrg_type='sum')

      call addmap(fldListFr(compatm)%flds, 'Faxa_snowc', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Faxa_snow', mrg_from=compatm, mrg_fld='Faxa_snowc', mrg_type='sum')
      call addmap(fldListFr(compatm)%flds, 'Faxa_snowl', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Faxa_snow', mrg_from=compatm, mrg_fld='Faxa_snowl', mrg_type='sum')

      ! from ice
      allocate(F_flds(4, 2))
      F_flds(1,:) = (/'Fioi_salt', 'Fioi_salt'/)
      F_flds(2,:) = (/'Si_ifrac', 'Si_ifrac'/) ! ice_fraction
      F_flds(3,:) = (/'Fioi_meltw', 'Fioi_meltw'/)
      F_flds(4,:) = (/'Fioi_melth', 'Fioi_melth'/) ! heat flux sea-ice to ocean
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
      call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', mrg_from=compice, mrg_fld='Fioi_taux', mrg_type='merge', mrg_fracname='ifrac')
      call addmap(fldListFr(compatm)%flds, 'Faxa_taux', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', mrg_from=compatm, mrg_fld='Faxa_taux', mrg_type='merge', mrg_fracname='ofrac')

      call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', mrg_from=compice, mrg_fld='Fioi_tauy', mrg_type='merge', mrg_fracname='ifrac')
      call addmap(fldListFr(compatm)%flds, 'Faxa_tauy', compocn, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', mrg_from=compatm, mrg_fld='Faxa_tauy', mrg_type='merge', mrg_fracname='ofrac')

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: state fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(S_flds(8))
      S_flds = (/'Sa_z', & ! inst_zonal_wind_height10m
                  'Sa_u', & ! inst_merid_wind_height10m
                  'Sa_v ', & ! inst_temp_height2m
                  'Sa_shum ', & ! inst_spec_humid_height2m
                  'Sa_tbot', & ! Sa_pslv
                  'Sa_pbot', &
                  'Sa_dens', &
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
      S_flds = (/'So_dhdx', & ! inst_zonal_wind_height10m
                 'So_dhdy', & ! inst_merid_wind_height10m
                 'So_t ', & ! inst_temp_height2m
                 'So_s ', & ! inst_spec_humid_height2m
                 'So_u', & ! Sa_pslv
                 'So_v', & ! Sa_pslv
                 'Fioo_q' /) ! inst_temp_height_surface
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
      F_flds(1,:) = (/'Faxa_swvdr ', 'Faxa_swvdr '/)
      F_flds(2,:) = (/'Faxa_swndr ', 'Faxa_swndr '/)
      F_flds(3,:) = (/'Faxa_swvdf', 'Faxa_swvdf'/)
      F_flds(4,:) = (/'Faxa_swndf', 'Faxa_swndf'/)
      F_flds(5,:) = (/'Faxa_lwdn', 'Faxa_lwdn'/)

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
      call addmap(fldListFr(compatm)%flds, 'Faxa_rainc', compice, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compice)%flds, 'Faxa_rain', mrg_from=compatm, mrg_fld='Faxa_rainc', mrg_type='sum')
      call addmap(fldListFr(compatm)%flds, 'Faxa_rainl', compice, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compice)%flds, 'Faxa_rain', mrg_from=compatm, mrg_fld='Faxa_rainl', mrg_type='sum')

      call addmap(fldListFr(compatm)%flds, 'Faxa_snowc', compice, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compice)%flds, 'Faxa_snow', mrg_from=compatm, mrg_fld='Faxa_snowc', mrg_type='sum')
      call addmap(fldListFr(compatm)%flds, 'Faxa_snowl', compice, mapconsf, 'one', 'unset')
      call addmrg(fldListTo(compice)%flds, 'Faxa_snow', mrg_from=compatm, mrg_fld='Faxa_snowl', mrg_type='sum')

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_init

    !-----------------------------------------------------------------------------

  end module esmFldsExchange_access_mod
