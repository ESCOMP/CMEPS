module med_enthalpy_mod
  use ESMF, only : ESMF_SUCCESS, ESMF_GridComp, ESMF_VMAllreduce, ESMF_REDUCE_SUM
  use shr_kind_mod, only : R8=>shr_kind_r8
  use shr_const_mod, only : tkfrz=>shr_const_tkfrz, cpfw=>shr_const_cpfw, cpice=>shr_const_cpice,&
       cpwv=>shr_const_cpwv, cpsw=>shr_const_cpsw, pi=>shr_const_pi
  use med_utils_mod         , only : chkerr        => med_utils_ChkErr
  use med_methods_mod       , only : FB_fldchk     => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_internalstate_mod, only : compocn, compatm, comprof, InternalState
  use med_internalstate_mod , only : logunit, maintask
  use perf_mod, only : t_startf, t_stopf
    

  implicit none
  public :: med_compute_enthalpy
  logical, public :: mediator_compute_enthalpy = .false.
  
  real(r8) :: global_htot_corr(1) = 0._r8
  character(*), parameter :: u_FILE_u  = &
       __FILE__
contains
  real(r8) function med_enthalpy_get_global_htot_corr()
    ! Just return the latest global_htot_corr
    med_enthalpy_get_global_htot_corr = global_htot_corr(1)
  end function med_enthalpy_get_global_htot_corr

  subroutine med_compute_enthalpy(is_local, rc)
    type(InternalState), intent(in) :: is_local
    integer, intent(out) :: rc

    ! local variables
    
    real(r8), pointer :: tocn(:), rain(:), snow(:), rofl(:), rofi(:), evap(:)
    real(r8), pointer :: rainl(:), rainc(:), tbot(:)
    real(r8), pointer :: snowl(:), snowc(:), ofrac(:)
    real(r8), pointer :: hrain(:), hsnow(:), hevap(:), hcond(:), hrofl(:), hrofi(:)
    real(r8), allocatable :: hcorr(:)
    real(r8), pointer :: areas(:)
    real(r8), parameter :: glob_area_inv = 1._r8 / (4._r8 * pi)
    real(r8) :: local_htot_corr(1)
    
    integer :: n, nmax
    character(len=*), parameter:: subname = "med_compute_enthalpy"

    call t_startf(subname)
    rc = ESMF_SUCCESS

    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', tocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nmax = size(tocn)
       
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm, compocn), 'Sa_tbot', tbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_rain', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_rain' , rain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm, compatm), 'Faxa_rainl', rainl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm, compatm), 'Faxa_rainl', rainc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(rain(nmax))
       rain = rainl + rainc
    endif
       
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrain', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrain', hrain, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hrain(nmax))
    endif
    
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap' , evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_evap' , evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hevap', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hevap', hevap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hevap(nmax))
    endif
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hcond', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hcond', hcond, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hcond(nmax))
    endif

    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_snow', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Faxa_snow' , snow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc' , snowc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl' , snowl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(snow(nmax))
       snow = snowc + snowl
    endif

    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', hsnow, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hsnow(nmax))
    endif
    
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , rofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl', rofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', hrofl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hrofl(nmax))
    endif

    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , rofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi', rofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    if(FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', rc)) then
       call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', hrofi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(hrofi(nmax))
    endif

    call fldbun_getdata1d(is_local%wrap%FBImp(compocn,compocn), 'So_omask', ofrac, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n=1,nmax
       ! for F cases (docn) tocn is non-zero over land and so ofrac must be included
       ! so that only ocean points are included in calculation
       ! Need max to ensure that will not have an enthalpy contribution if the water is below 0C
       hrain(n)  = max((tbot(n) - tkfrz), 0._r8) * rain(n)  * cpfw  * ofrac(n)
       hsnow(n)  = min((tbot(n) - tkfrz), 0._r8) * snow(n)  * cpice * ofrac(n)
       hevap(n)  = (tocn(n) - tkfrz) * min(evap(n), 0._r8)  * cpwv  * ofrac(n)
       hcond(n)  = (tocn(n) - tkfrz) * max(evap(n), 0._r8)  * cpwv  * ofrac(n)
       hrofl(n)  = max((tocn(n) - tkfrz), 0._r8) * rofl(n)  * cpfw  * ofrac(n)
       hrofi(n)  = min((tocn(n) - tkfrz), 0._r8) * rofi(n)  * cpice * ofrac(n)
    end do
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_rain', rc)) deallocate(rain)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Faxa_snow', rc)) deallocate(snow)

    
    ! Determine enthalpy correction factor that will be added to the sensible heat flux sent to the atm
    ! Areas here in radians**2 - this is an instantaneous snapshot that will be sent to the atm - only
    ! need to calculate this if data is sent back to the atm
    
    if (FB_fldchk(is_local%wrap%FBExp(compatm), 'Faxx_sen', rc=rc)) then
       allocate(hcorr(nmax))
       areas => is_local%wrap%mesh_info(compocn)%areas
       do n = 1,nmax
          hcorr(n) = (hrain(n) + hsnow(n) + hcond(n) + hevap(n) + hrofl(n) + hrofi(n)) * &
               areas(n) * glob_area_inv
       end do

       ! Determine sum of enthalpy correction for each hcorr index locally
       local_htot_corr(1) = 0._r8
       do n = 1,size(hcorr)
          local_htot_corr(1) = local_htot_corr(1) + hcorr(n)
       end do
       call ESMF_VMAllreduce(is_local%wrap%vm, senddata=local_htot_corr, recvdata=global_htot_corr, count=1, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (maintask) write(logunit, '(a,a,f21.13)') trim(subname),' global enthalpy correction: ',global_htot_corr(1)
       deallocate(hcorr)
    endif
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hsnow', rc)) deallocate(hsnow)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofl', rc)) deallocate(hrofl)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrofi', rc)) deallocate(hrofi)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hrain', rc)) deallocate(hrain)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hevap', rc)) deallocate(hevap)
    if(.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_hcond', rc)) deallocate(hcond)
    
    call t_stopf(subname)

  end subroutine med_compute_enthalpy

end module med_enthalpy_mod
