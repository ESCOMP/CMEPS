module esmFldsExchange_hafs_mod

  use ESMF
  use NUOPC
  use med_utils_mod, only : chkerr => med_utils_chkerr
  use med_kind_mod,  only : CX=>SHR_KIND_CX
  use med_kind_mod,  only : CS=>SHR_KIND_CS
  use med_kind_mod,  only : CL=>SHR_KIND_CL
  use med_kind_mod,  only : R8=>SHR_KIND_R8
  use esmflds,       only : compmed
  use esmflds,       only : compatm
  use esmflds,       only : compocn
  use esmflds,       only : ncomps
  use esmflds,       only : fldListTo
  use esmflds,       only : fldListFr
  use esmFlds,       only : coupling_mode

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

  type gcomp_attr
    character(len=CX)   :: atm2ocn_fmap='unset'
    character(len=CX)   :: atm2ocn_smap='unset'
    character(len=CX)   :: atm2ocn_vmap='unset'
    character(len=CX)   :: ocn2atm_fmap='unset'
    character(len=CX)   :: ocn2atm_smap='unset'
    character(len=CS)   :: mapnorm     ='one'
  end type

!===============================================================================
contains
!===============================================================================

  subroutine esmFldsExchange_hafs(gcomp, phase, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    if (phase == 'advertise') then
      call esmFldsExchange_hafs_advt(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'fieldcheck') then
      call esmFldsExchange_hafs_fchk(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'initialize') then
      call esmFldsExchange_hafs_init(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogSetError(ESMF_FAILURE, &
         msg=trim(subname)//": Phase is set to "//trim(phase), &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_advt(gcomp, phase, rc)

    use esmFlds               , only : addfld => med_fldList_AddFld

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
    type(gcomp_attr)    :: hafs_attr
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_advt)'
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
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    call addfld(fldListFr(compocn)%flds, 'So_omask')

    !----------------------------------------------------------
    ! to med: frac from components
    !----------------------------------------------------------
    call addfld(fldListTo(compatm)%flds, 'So_ofrac')

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to atm: surface temperatures from ocn
    ! ---------------------------------------------------------------------
    allocate(S_flds(1))
    S_flds = (/'So_t'/) ! sea_surface_temperature
    do n = 1,size(S_flds)
       fldname = trim(S_flds(n))
       call addfld(fldListFr(compocn)%flds, trim(fldname))
       call addfld(fldListTo(compatm)%flds, trim(fldname))
    end do
    deallocate(S_flds)

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    allocate(S_flds(6))
    S_flds = (/'Sa_u10m', & ! inst_zonal_wind_height10m
               'Sa_v10m', & ! inst_merid_wind_height10m
               'Sa_t2m ', & ! inst_temp_height2m
               'Sa_q2m ', & ! inst_spec_humid_height2m
               'Sa_pslv', & ! inst_pres_height_surface
               'Sa_tskn' /) ! inst_temp_height_surface
    do n = 1,size(S_flds)
       fldname = trim(S_flds(n))
       call addfld(fldListFr(compatm)%flds, trim(fldname))
       call addfld(fldListTo(compocn)%flds, trim(fldname))
    end do
    deallocate(S_flds)

    ! ---------------------------------------------------------------------
    ! to ocn: flux fields
    ! ---------------------------------------------------------------------
    allocate(F_flds(7,2))
    F_flds(1,:) = (/'Faxa_taux ','Faxa_taux '/) ! mean_zonal_moment_flx_atm
    F_flds(2,:) = (/'Faxa_tauy ','Faxa_tauy '/) ! mean_merid_moment_flx_atm
    F_flds(3,:) = (/'Faxa_rain ','Faxa_rain '/) ! mean_prec_rate
    F_flds(4,:) = (/'Faxa_swnet','Faxa_swnet'/) ! mean_net_sw_flx
    F_flds(5,:) = (/'Faxa_lwnet','Faxa_lwnet'/) ! mean_net_lw_flx
    F_flds(6,:) = (/'Faxa_sen  ','Faxa_sen  '/) ! mean_sensi_heat_flx
    F_flds(7,:) = (/'Faxa_lat  ','Faxa_lat  '/) ! mean_laten_heat_flx
    do n = 1,size(F_flds,1)
       fldname1 = trim(F_flds(n,1))
       fldname2 = trim(F_flds(n,2))
       call addfld(fldListFr(compatm)%flds, trim(fldname1))
       call addfld(fldListTo(compocn)%flds, trim(fldname2))
    end do
    deallocate(F_flds)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_advt

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_fchk(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_fchk)'
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

  end subroutine esmFldsExchange_hafs_fchk

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_init(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapnstod, mapnstod_consd
    use esmflds               , only : mapfillv_bilnr
    use esmflds               , only : mapnstod_consf

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
    type(gcomp_attr)    :: hafs_attr
    character(len=CS), allocatable :: flds(:)
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=CS), allocatable :: suffix(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_init)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------
    ! Merging arguments:
    ! mrg_fromN = source component index that for the field to be merged
    ! mrg_fldN  = souce field name to be merged
    ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum',
    !                         'sum_with_weights', 'merge')
    ! NOTE:
    ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn
    ! fluxes or for ocn albedos
    !
    ! NOTE:
    ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in
    ! the atm/ocn flux calculations. Input fields required from either the atm
    ! or the ocn for these computation will use the logical 'use_med_aoflux'
    ! below. This is used to determine mappings between the atm and ocn needed
    ! for these computations.
    !--------------------------------------

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)
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
       if (fldchk(is_local%wrap%FBExp(compatm),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compocn)%flds, trim(fldname), compatm, &
               mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%ocn2atm_smap)
          call addmrg(fldListTo(compatm)%flds, trim(fldname), &
               mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(S_flds)

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    allocate(S_flds(6))
    S_flds = (/'Sa_u10m', & ! inst_zonal_wind_height10m
               'Sa_v10m', & ! inst_merid_wind_height10m
               'Sa_t2m ', & ! inst_temp_height2m
               'Sa_q2m ', & ! inst_spec_humid_height2m
               'Sa_pslv', & ! inst_pres_height_surface
               'Sa_tskn' /) ! inst_temp_height_surface
    do n = 1,size(S_flds)
       fldname = trim(S_flds(n))
       if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
          ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, &
               mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
       end if
    end do
    deallocate(S_flds)

    ! ---------------------------------------------------------------------
    ! to ocn: flux fields
    ! ---------------------------------------------------------------------
    allocate(F_flds(7,2))
    F_flds(1,:) = (/'Faxa_taux ','Faxa_taux '/) ! mean_zonal_moment_flx_atm
    F_flds(2,:) = (/'Faxa_tauy ','Faxa_tauy '/) ! mean_merid_moment_flx_atm
    F_flds(3,:) = (/'Faxa_rain ','Faxa_rain '/) ! mean_prec_rate
    F_flds(4,:) = (/'Faxa_swnet','Faxa_swnet'/) ! mean_net_sw_flx
    F_flds(5,:) = (/'Faxa_lwnet','Faxa_lwnet'/) ! mean_net_lw_flx
    F_flds(6,:) = (/'Faxa_sen  ','Faxa_sen  '/) ! mean_sensi_heat_flx
    F_flds(7,:) = (/'Faxa_lat  ','Faxa_lat  '/) ! mean_laten_heat_flx
    do n = 1,size(F_flds,1)
       fldname1 = trim(F_flds(n,1))
       fldname2 = trim(F_flds(n,2))
       if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname2),rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname1),rc=rc) &
         ) then
          call addmap(fldListFr(compatm)%flds, trim(fldname1), compocn, &
               mapfillv_bilnr, hafs_attr%mapnorm, hafs_attr%atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, trim(fldname2), &
               mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')
       end if
    end do
    deallocate(F_flds)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_init

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_hafs_attr(gcomp, hafs_attr, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    type(gcomp_attr) , intent(inout) :: hafs_attr
    integer          , intent(inout) :: rc

    ! local variables:
    character(32)       :: cname
    integer             :: verbosity, diagnostic
    character(len=CL)   :: cvalue
    logical             :: isPresent
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_attr)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call NUOPC_CompGet(gcomp, name=cname, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Normalization type
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='normalization', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='normalization', &
          value=hafs_attr%mapnorm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
          value=hafs_attr%ocn2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
          value=hafs_attr%ocn2atm_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ocn
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
          value=hafs_attr%atm2ocn_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       value=hafs_attr%atm2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
          value=hafs_attr%atm2ocn_vmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Log Attribute Settings
    if (btest(verbosity,16)) then
       write(cvalue,"(I0)") verbosity
       call ESMF_LogWrite(trim(subname)//': Verbosity        = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       write(cvalue,"(I0)") diagnostic
       call ESMF_LogWrite(trim(subname)//': Diagnostic       = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': normalization    = '// &
          trim(hafs_attr%mapnorm), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_smapname = '// &
          trim(hafs_attr%ocn2atm_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_fmapname = '// &
          trim(hafs_attr%ocn2atm_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_fmapname = '// &
          trim(hafs_attr%atm2ocn_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_smapname = '// &
          trim(hafs_attr%atm2ocn_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_vmapname = '// &
          trim(hafs_attr%atm2ocn_vmap), ESMF_LOGMSG_INFO)
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_hafs_attr

end module esmFldsExchange_hafs_mod
