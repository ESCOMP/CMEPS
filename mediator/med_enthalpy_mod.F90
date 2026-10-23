module med_enthalpy_mod

  !-----------------------------------------------------------------------------
  ! Enthalpy flux calculations shared by the ocn and atm prep phases.
  !
  ! Two paths are supported, selected by the component_computes_enthalpy_flux
  ! attribute:
  !   'med' - the mediator computes the enthalpy flux (coupling_mode 'cesm')
  !   'atm' - the enthalpy flux comes from the prognostic atm (coupling_mode 'noresm')
  !
  ! The global corrections computed here (global_htot_corr, global_hrof_corr) are
  ! applied to the atm export state in med_phases_prep_atm.
  !-----------------------------------------------------------------------------

  use shr_log_mod           , only : shr_log_error
  use med_kind_mod          , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use med_utils_mod         , only : chkerr       => med_utils_ChkErr
  use med_global_sums_mod   , only : med_global_sums
  use med_methods_mod       , only : FB_fldchk    => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr => med_methods_FB_GetFldPtr
  use med_internalstate_mod , only : InternalState
  use med_internalstate_mod , only : compocn, compatm

  implicit none
  private

  public :: med_enthalpy_init                ! called from med_phases_prep_ocn_init
  public :: med_enthalpy_med_computation     ! mediator computes the enthalpy flux
  public :: med_enthalpy_atm_computation     ! enthalpy flux comes from the atm
  public :: med_enthalpy_correction          ! global enthalpy correction
  public :: med_enthalpy_runoff              ! global enthalpy of runoff

  ! Which component computes the enthalpy flux - 'med', 'atm' or 'none'
  character(len=CS), public :: component_computes_enthalpy_flux = 'unset'

  ! Global corrections applied to the atm export state in med_phases_prep_atm
  real(r8), public :: global_htot_corr = 0._r8  ! enthalpy correction
  real(r8), public :: global_hrof_corr = 0._r8  ! enthalpy of run-off

  character(*), parameter :: u_FILE_u = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_enthalpy_init(gcomp, rc)

    ! Determine which component computes the enthalpy flux.

    use ESMF , only : ESMF_GridComp, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: cvalue
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname='(med_enthalpy_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name="component_computes_enthalpy_flux", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       component_computes_enthalpy_flux = trim(cvalue)
    else
       component_computes_enthalpy_flux = 'none'
    end if

  end subroutine med_enthalpy_init

  !-----------------------------------------------------------------------------
  subroutine med_enthalpy_med_computation(gcomp, rc)

    ! Compute enthalpy associated with rain, snow, condensation and liquid river & glc runoff.
    ! The sea-ice model already accounts for the enthalpy flux (as part of melth), so
    ! enthalpy from meltw **is not** included below.
    ! This is the coupling_mode='cesm' path, where the mediator computes the enthalpy flux.

    use ESMF                    , only : ESMF_GridComp, ESMF_SUCCESS
    use med_constants_mod       , only : shr_const_cpsw, shr_const_tkfrz, shr_const_pi
    use med_constants_mod       , only : shr_const_cpice

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)   :: is_local
    integer               :: n
    real(r8)              :: glob_area_inv
    real(r8), pointer     :: tocn(:)
    real(r8), pointer     :: rain(:), hrain(:)
    real(r8), pointer     :: snow(:), hsnow(:)
    real(r8), pointer     :: evap(:), hevap(:)
    real(r8), pointer     :: hcond(:)
    real(r8), pointer     :: rofl(:), hrofl(:)
    real(r8), pointer     :: rofi(:), hrofi(:)
    real(r8), pointer     :: rofl_glc(:), hrofl_glc(:)
    real(r8), pointer     :: rofi_glc(:), hrofi_glc(:)
    real(r8), pointer     :: areas(:)
    real(r8), allocatable :: hcorr(:)
    character(len=*), parameter :: subname='(med_enthalpy_med_computation)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_rain'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrain'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_snow'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hsnow'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hevap'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hcond'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi'      , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi'     , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc'  , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl_glc' , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc'  , rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi_glc' , rc=rc)       &
         ) then
       ! Error check
       if (trim(component_computes_enthalpy_flux) /= 'med') then
          call shr_log_error(trim(subname)//' ERROR: component_computes_enthalpy_flux must be set to med', rc=rc)
          return
       end if
       call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', tocn, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_rain' , rain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrain', hrain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap' , evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hevap', hevap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hcond', hcond, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_snow' , snow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', hsnow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', hrofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , rofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', hrofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Forr_rofl_glc' , rofl_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofl_glc', hrofl_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Forr_rofi_glc' , rofi_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofi_glc', hrofi_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do n = 1,size(tocn)
          ! Need max to ensure that will not have an enthalpy contribution if the water is below 0C
          hrain(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * rain(n)  * shr_const_cpsw
          hsnow(n)  = min((tocn(n) - shr_const_tkfrz), 0._r8) * snow(n)  * shr_const_cpsw
          hevap(n)  = (tocn(n) - shr_const_tkfrz) * min(evap(n), 0._r8)   * shr_const_cpsw
          hcond(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * max(evap(n), 0._r8)  * shr_const_cpsw
          hrofl(n)  = max((tocn(n) - shr_const_tkfrz), 0._r8) * rofl(n)  * shr_const_cpsw
          hrofl_glc(n) = max((tocn(n) - shr_const_tkfrz), 0._r8) * rofl_glc(n)  * shr_const_cpsw
          ! −10 C is a reasonable bulk temperature assumption for iceberg/land-ice runoff
          hrofi(n)  = -10._r8 * rofi(n)  * shr_const_cpice
          hrofi_glc(n) = -10._r8 * rofi_glc(n)  * shr_const_cpice
       end do

       if (FB_fldchk(is_local%wrap%FBExp(compatm), 'Faxx_sen', rc=rc)) then
          ! Determine enthalpy correction factor that will be added to the sensible heat flux sent to the atm
          ! Areas here in radians**2 - this is an instantaneous snapshot that will be sent to the atm - only
          ! need to calculate this if data is sent back to the atm
          allocate(hcorr(size(tocn)))
          glob_area_inv = 1._r8 / (4._r8 * shr_const_pi)
          areas => is_local%wrap%mesh_info(compocn)%areas
          do n = 1,size(tocn)
             hcorr(n) = (hrain(n) + hsnow(n) + hcond(n) + hevap(n) + hrofl(n) + hrofi(n) + hrofl_glc(n) + hrofi_glc(n)) * &
                  areas(n) * glob_area_inv
          end do

          call med_enthalpy_correction(gcomp, hcorr, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          deallocate(hcorr)

       end if
    end if ! condition for using global energy fixer

  end subroutine med_enthalpy_med_computation

  !-----------------------------------------------------------------------------
  subroutine med_enthalpy_atm_computation(gcomp, rc)

    ! Apply the enthalpy flux obtained from the prognostic atm.
    ! This is the coupling_mode='noresm' path, where the atm computes the enthalpy flux.

    use ESMF                    , only : ESMF_GridComp, ESMF_SUCCESS
    use med_constants_mod       , only : shr_const_tkfrz, shr_const_pi
    use med_constants_mod       , only : shr_const_cpice, shr_const_cpfw

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)   :: is_local
    integer               :: n
    real(r8), pointer     :: tocn(:)
    real(r8), pointer     :: rofl(:)
    real(r8), pointer     :: rofi(:)
    real(r8), pointer     :: rofl_glc(:)
    real(r8), pointer     :: rofi_glc(:)
    real(r8), pointer     :: areas(:)
    real(r8), allocatable :: hcorr(:)
    real(r8), pointer     :: dataptr(:)
    real(r8), pointer     :: Faxa_hmat(:)
    real(r8), pointer     :: Faxa_hlat(:)
    real(r8), allocatable :: hrof2atm(:)
    real(r8)              :: ocean_htot_corr
    real(r8)              :: ocean_atot_corr
    real(r8), allocatable :: hrof(:)

    ! if separate_varlat is true then do global ocean average for
    ! hmat_oa only for the net-mass part, and pass in hmat only local
    ! variable latent heat correction part
    logical, parameter    :: separate_varlat=.true.
    real(r8), allocatable :: acorr(:)
    character(len=*), parameter :: subname='(med_enthalpy_atm_computation)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if( FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_hmat', rc=rc) .and.  &
         FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_hlat', rc=rc)) then
       if (trim(component_computes_enthalpy_flux) /= 'atm') then
          call shr_log_error(trim(subname)//' ERROR: component_computes_enthalpy_flux must be set to atm', rc=rc)
          return
       end if
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_hmat', Faxa_hmat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_hlat', Faxa_hlat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', tocn, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       areas => is_local%wrap%mesh_info(compocn)%areas

       ! if separate_varlat is true then do global ocean average only for the
       ! net-mass part, and pass in as hmat only the local variable latent heat correction part
       if (separate_varlat) then
          !-----------------------
          ! Determine enthalpy due to ocean river input
          !-----------------------
          allocate(hrof (size(tocn)))
          if ( FB_fldchk(is_local%wrap%FBExp(compocn),'Foxx_rofl'    ,rc=rc) .and. &
               FB_fldchk(is_local%wrap%FBExp(compocn),'Foxx_rofi'    ,rc=rc) .and. &
               FB_fldchk(is_local%wrap%FBExp(compocn),'Forr_rofl_glc',rc=rc) .and. &
               FB_fldchk(is_local%wrap%FBExp(compocn),'Forr_rofi_glc',rc=rc) ) then

             call FB_GetFldPtr(is_local%wrap%FBExp(compocn),'Foxx_rofl'    , rofl    ,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBExp(compocn),'Foxx_rofi'    , rofi    ,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBExp(compocn),'Forr_rofl_glc', rofl_glc,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBExp(compocn),'Forr_rofi_glc', rofi_glc,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(tocn)
                hrof(n) = shr_const_cpfw  * (tocn(n) - shr_const_tkfrz) * rofl(n) &
                     + shr_const_cpice * (tocn(n) - shr_const_tkfrz) * rofi(n) &
                     + shr_const_cpfw  * (tocn(n) - shr_const_tkfrz) * rofl_glc(n) &
                     + shr_const_cpice * (tocn(n) - shr_const_tkfrz) * rofi_glc(n)
             enddo
          else
             do n = 1,size(tocn)
                hrof(n) = 0._r8
             enddo
          endif

          ! send back to atm if requested by atm
          if (FB_fldchk(is_local%wrap%FBExp(compatm), 'Faxx_hrof', rc=rc)) then
             allocate(hrof2atm(size(tocn)))
             hrof2atm(:) = hrof(:)*areas(:) / (4._r8 * shr_const_pi)

             ! determine module variable global_hrof_corr in med_enthalpy_mod
             call med_enthalpy_runoff(gcomp, hrof2atm, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       !-----------------------
       ! Compute Faxa_hmat_oa
       !-----------------------
       ! Determine hcorr and acorr
       allocate(hcorr(size(tocn)))
       allocate(acorr(size(tocn)))
       if (separate_varlat) then
          do n = 1,size(tocn)
             hcorr(n) = areas(n) *(Faxa_hmat(n) - Faxa_hlat(n) + hrof(n))
             acorr(n) = areas(n)
          end do
       else
          do n = 1,size(tocn)
             hcorr(n) = areas(n) * Faxa_hmat(n)
             acorr(n) = areas(n)
          end do
       endif
       deallocate(hrof)

       ! Compute global integral of hcorr - ocean_oa_htot
       call med_global_sums(gcomp, hcorr, ocean_htot_corr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(hcorr)

       ! Compute global integral of acorr - ocean_atot_corr
       call med_global_sums(gcomp, acorr, ocean_atot_corr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(acorr)

       ! Set value of Faxa_hmat_oa to ratio of ocean_htot_corr and ocean_atot_corr
       call FB_getfldptr(is_local%wrap%FBExp(compocn), 'Faxa_hmat_oa', dataptr, rc=rc)
       if (ocean_atot_corr > 0._r8) then
          dataptr(:) = ocean_htot_corr/ocean_atot_corr
       end if

       !-----------------------
       ! replace full material enthalpy flux with variable latent
       ! heats part only in pointer to ocean export
       !-----------------------

       ! might add another coupling field later but may not be strictly necessary
       if (separate_varlat) then
          do n = 1,size(tocn)
             Faxa_hmat(n) = Faxa_hlat(n)
          end do
       else
          do n = 1,size(tocn)
             Faxa_hmat(n) = 0._r8 ! avoid applying twice for some ocean components such as BLOM
          end do
       endif
    endif

  end subroutine med_enthalpy_atm_computation


  !-----------------------------------------------------------------------------
  subroutine med_enthalpy_correction (gcomp, hcorr, rc)

    use ESMF , only : ESMF_GridComp, ESMF_SUCCESS

    ! Enthalpy correction term calculation called by med_enthalpy_med_computation
    ! Note that this is only called if the following fields are in FBExp(compocn)
    ! 'Faxa_rain','Foxx_hrain','Faxa_snow' ,'Foxx_hsnow',
    ! 'Foxx_evap','Foxx_hevap','Foxx_hcond','Foxx_rofl',
    ! 'Foxx_hrofl','Foxx_rofi','Foxx_hrofi','Foxx_rofl_glc',
    ! 'Foxx_hrofl_glc','Foxx_rofi_glc','Foxx_hrofi_glc'

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    real(r8)            , intent(in)  :: hcorr(:)
    integer             , intent(out) :: rc
    !---------------------------------------

    rc = ESMF_SUCCESS

    call med_global_sums(gcomp, hcorr, global_htot_corr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_enthalpy_correction

  !-----------------------------------------------------------------------------
  subroutine med_enthalpy_runoff(gcomp, hcorr, rc)

    use ESMF , only : ESMF_GridComp, ESMF_SUCCESS

    ! Enthalpy of runoff, called by med_enthalpy_atm_computation
    ! Note that this is only called if the following fields are in FBExp(compocn)
    ! - Faxa_hmat, Faxa_hlat
    ! The result (Faxx_hrof) is sent back to the atm in subroutine med_phases_prep_atm

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    real(r8)            , intent(in)  :: hcorr(:)
    integer             , intent(out) :: rc
    !---------------------------------------

    rc = ESMF_SUCCESS

    call med_global_sums(gcomp, hcorr, global_hrof_corr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_enthalpy_runoff

end module med_enthalpy_mod
