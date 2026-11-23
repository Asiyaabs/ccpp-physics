!> \file water_isotopes.F90
!! Water isotopes (HDO, H2-18O): fractionation factors and isotope helper functions.
!!
!! Plain Fortran module (no CCPP metadata), for any scheme that moves water: PBL, convection,
!! microphysics, photochemistry, and later the land (Noah-MP). A scheme that uses it lists
!! ../../tools/water_isotopes.F90 in "dependencies" of its .meta file.
!!
!! Ported from iCESM1.2 (public domain, CCSM copyright notice):
!!   models/csm_share/shr/water_isotopes.F90   (D. Noone, J. Nusbaumer)
!!   models/csm_share/shr/water_types.F90      (wtype_get_alpha, C. Bardeen)
!!   models/atm/cam/src/physics/cam/water_tracers.F90 (wtrc_efac, wtrc_dqequil,
!!     wtrc_liqvap_equil, wtrc_vap_distil, wtrc_equil_time, stewart_isoevap; D. Noone, J. Nusbaumer)
!! Formulas, coefficients and tuning constants are the iCESM1.2 defaults. Changes from iCESM:
!!  - species ids follow water_tracers_config (1 H2O, 2 HDO, 3 H218O; no H216O);
!!  - functions return values (no CAM state, no qsat calls, no write/endrun); inputs that
!!    iCESM computed inside (saturation vapour pressure, relative humidity) are arguments;
!!  - no qtiny in ratios: ratio = qiso/qtot for qtot > 0, else Rstd;
!!  - no NOFRAC/NOKIN compile flags: the caller decides whether to apply alpha.
!!
!! Conventions:
!!  - alpha = R_condensate / R_vapour at equilibrium (> 1 for HDO and H218O).
!!    In wtr_transfer (moved tracer = alpha * src ratio * dq) use alpha for vapour -> liquid/ice
!!    and 1/alpha for liquid -> vapour; wiso_alpha_phase returns the right one per direction.
!!  - Rstd = 1 for all species (iCESM default, "best numerics"): tracers are carried as
!!    R/R_VSMOW, so delta = 1000 (R - 1). The natural ratios are in wiso_get_rnat.
module water_isotopes

  use machine,  only: kind_phys
  use physcons, only: con_g, con_rd, con_rv, con_ttp, karman, rhowater

  implicit none
  private

  ! species and phases
  public :: wiso_nspec, wiso_h2o, wiso_hdo, wiso_h218o, wiso_get_ispec
  ! equilibrium fractionation
  public :: wiso_alpl, wiso_alpl_majoube, wiso_alpi, wiso_alpha_phase
  ! kinetic fractionation
  public :: wiso_akel, wiso_akci, wiso_ssatf, wiso_heff, wiso_kmol, wiso_kmolv10
  ! surface (ocean) evaporation
  public :: wiso_flxoce, wiso_cg_ratio, wiso_get_roce
  ! rain evaporation and equilibration
  public :: wiso_efac, wiso_dqequil, wiso_liqvap_equil, wiso_vap_distil, &
            wiso_equil_time, wiso_drop_radius, wiso_stewart_ratio
  ! ratios, delta, constants
  public :: wiso_ratio, wiso_delta, wiso_get_rstd, wiso_get_rnat, wiso_get_fisub, wiso_get_difrm

  !> Species ids (same as water_tracers_config: wtr_species_h2o/hdo/h218o)
  integer, parameter :: wiso_nspec = 3
  integer, parameter :: wiso_h2o   = 1
  integer, parameter :: wiso_hdo   = 2
  integer, parameter :: wiso_h218o = 3
  character(len=5), parameter :: spnam(wiso_nspec) = (/ 'H2O  ', 'HDO  ', 'H218O' /)

  !> Phase ids for wiso_alpha_phase (same order as water_tracers_config wtr_phase_names)
  integer, parameter :: ph_qv = 1, ph_ql = 2, ph_qi = 3, ph_qr = 4, ph_qs = 5, ph_qg = 6

  ! ---- Tunable parameters (iCESM1.2 values) ----
  real(kind_phys), parameter :: dkfac  = 0.58_kind_phys     !< diffusive evaporation kinetic power law
  real(kind_phys), parameter :: tkini  = 253.15_kind_phys   !< kinetic effects for ice below this T (Jouzel & Merlivat 1984)
  real(kind_phys), parameter :: recrit = 1.0_kind_phys      !< critical Reynolds number for kmol
  real(kind_phys), parameter :: fsata  = 1.000_kind_phys    !< supersaturation s = fsata + fsatb*(T - T0) (Hoffmann)
  real(kind_phys), parameter :: fsatb  = -0.002_kind_phys   !< tuned to Antarctic d-excess (iCESM; Hoffmann -0.003)
  real(kind_phys), parameter :: ssatmx = 2.00_kind_phys     !< maximum supersaturation
  real(kind_phys), parameter :: fkhum  = 0.25_kind_phys     !< effective humidity factor
  real(kind_phys), parameter :: tzero  = con_ttp            !< triple point (supersaturation reference)
  real(kind_phys), parameter :: phi_stw = 0.9_kind_phys     !< RH enhancement near drops (Bony et al. 2008), Stewart

  ! ---- Physical constants per species (H2O, HDO, H218O) ----
  real(kind_phys), parameter :: fisub(wiso_nspec) = (/ 1._kind_phys, 2._kind_phys, 1._kind_phys /)    !< isotope atoms per molecule
  real(kind_phys), parameter :: difrm(wiso_nspec) = &  !< diffusivity ratio Di/D (Merlivat 1978)
       (/ 1._kind_phys, 0.9757_kind_phys, 0.9727_kind_phys /)
  real(kind_phys), parameter :: rnat(wiso_nspec)  = &  !< natural (VSMOW) ratios
       (/ 1._kind_phys, 155.76e-6_kind_phys, 2005.20e-6_kind_phys /)
  real(kind_phys), parameter :: rstd(wiso_nspec)  = (/ 1._kind_phys, 1._kind_phys, 1._kind_phys /)   !< model standard ratio
  real(kind_phys), parameter :: boce(wiso_nspec)  = (/ 1._kind_phys, 1._kind_phys, 1._kind_phys /)   !< ocean surface enrichment

  ! Ocean surface kinetic fractionation, Merlivat & Jouzel (1979) 10 m wind method
  real(kind_phys), parameter :: aksmc(wiso_nspec) = (/ 0._kind_phys, 0.00528_kind_phys,   0.006_kind_phys    /)
  real(kind_phys), parameter :: akrfa(wiso_nspec) = (/ 0._kind_phys, 0.2508e-3_kind_phys, 0.285e-3_kind_phys /)
  real(kind_phys), parameter :: akrfb(wiso_nspec) = (/ 0._kind_phys, 0.7216e-3_kind_phys, 0.82e-3_kind_phys  /)

  ! Liquid/vapour, Horita & Wesolowski (1994) (iCESM default)
  real(kind_phys), parameter :: alpal(wiso_nspec) = (/ 0._kind_phys, 1158.8e-12_kind_phys,  0.35041e+6_kind_phys /)
  real(kind_phys), parameter :: alpbl(wiso_nspec) = (/ 0._kind_phys, -1620.1e-9_kind_phys, -1.6664e+3_kind_phys  /)
  real(kind_phys), parameter :: alpcl(wiso_nspec) = (/ 0._kind_phys, 794.84e-6_kind_phys,   6.7123_kind_phys     /)
  real(kind_phys), parameter :: alpdl(wiso_nspec) = (/ 0._kind_phys, -161.04e-3_kind_phys, -7.685e-3_kind_phys   /)
  real(kind_phys), parameter :: alpel(wiso_nspec) = (/ 0._kind_phys, 2.9992e+6_kind_phys,   0._kind_phys         /)

  ! Liquid/vapour, Majoube (1971a) (iCESM alternative)
  real(kind_phys), parameter :: majal(wiso_nspec) = (/ 0._kind_phys, 24.844e+3_kind_phys,  1.137e+3_kind_phys   /)
  real(kind_phys), parameter :: majbl(wiso_nspec) = (/ 0._kind_phys, -76.248_kind_phys,   -0.4156_kind_phys     /)
  real(kind_phys), parameter :: majcl(wiso_nspec) = (/ 0._kind_phys, 52.612e-3_kind_phys, -2.0667e-3_kind_phys  /)

  ! Ice/vapour, Merlivat & Nief (1967) for HDO, Majoube (1971b) for H218O
  real(kind_phys), parameter :: alpai(wiso_nspec) = (/ 0._kind_phys, 16289._kind_phys,     0._kind_phys         /)
  real(kind_phys), parameter :: alpbi(wiso_nspec) = (/ 0._kind_phys, 0._kind_phys,         11.839_kind_phys     /)
  real(kind_phys), parameter :: alpci(wiso_nspec) = (/ 0._kind_phys, -9.45e-2_kind_phys,  -28.224e-3_kind_phys  /)

contains

  ! =====================================================================
  ! Species
  ! =====================================================================

  !> Species id from its name ('H2O', 'HDO', 'H218O'); 0 if unknown.
  pure integer function wiso_get_ispec(name)
    character(len=*), intent(in) :: name
    integer :: m
    wiso_get_ispec = 0
    do m = 1, wiso_nspec
       if (trim(adjustl(name)) == trim(spnam(m))) wiso_get_ispec = m
    end do
  end function wiso_get_ispec

  ! =====================================================================
  ! Equilibrium fractionation
  ! =====================================================================

  !> Liquid/vapour equilibrium fractionation, Horita & Wesolowski (1994). 1 for H2O.
  elemental real(kind_phys) function wiso_alpl(isp, tk)
    integer,         intent(in) :: isp  !< species id
    real(kind_phys), intent(in) :: tk   !< temperature (K)
    if (isp == wiso_h2o) then
       wiso_alpl = 1._kind_phys
    else if (isp == wiso_hdo) then
       wiso_alpl = exp(alpal(isp)*tk**3 + alpbl(isp)*tk**2 + alpcl(isp)*tk + alpdl(isp) + alpel(isp)/tk**3)
    else
       wiso_alpl = exp(alpal(isp)/tk**3 + alpbl(isp)/tk**2 + alpcl(isp)/tk + alpdl(isp))
    end if
  end function wiso_alpl

  !> Liquid/vapour equilibrium fractionation, Majoube (1971a). 1 for H2O.
  elemental real(kind_phys) function wiso_alpl_majoube(isp, tk)
    integer,         intent(in) :: isp  !< species id
    real(kind_phys), intent(in) :: tk   !< temperature (K)
    if (isp == wiso_h2o) then
       wiso_alpl_majoube = 1._kind_phys
    else
       wiso_alpl_majoube = exp(majal(isp)/tk**2 + majbl(isp)/tk + majcl(isp))
    end if
  end function wiso_alpl_majoube

  !> Ice/vapour equilibrium fractionation, Merlivat & Nief (1967) HDO, Majoube (1971b) H218O.
  elemental real(kind_phys) function wiso_alpi(isp, tk)
    integer,         intent(in) :: isp  !< species id
    real(kind_phys), intent(in) :: tk   !< temperature (K)
    if (isp == wiso_h2o) then
       wiso_alpi = 1._kind_phys
    else
       wiso_alpi = exp(alpai(isp)/tk**2 + alpbi(isp)/tk + alpci(isp))
    end if
  end function wiso_alpi

  !> Fractionation factor for water moving from phase src to phase dst (phase ids as
  !! water_tracers_config: 1 qv, 2 ql, 3 qi, 4 qr, 5 qs, 6 qg), as iCESM wtype_get_alpha:
  !!   vapour -> liquid: alpl (with kinetic: akel(rh)); vapour -> ice: alpi (with kinetic: akci)
  !!   liquid -> vapour: 1/alpl (with kinetic: 1/akel(rh)); ice -> vapour: 1 (no fractionation)
  !!   between condensed phases (freezing, melting, riming, autoconversion, ...): 1
  !! The result is the alpha of wtr_transfer (moved tracer = alpha * src ratio * dq).
  elemental real(kind_phys) function wiso_alpha_phase(isp, src, dst, tk, rh, kinetic)
    integer,         intent(in) :: isp      !< species id
    integer,         intent(in) :: src      !< source phase id
    integer,         intent(in) :: dst      !< destination phase id
    real(kind_phys), intent(in) :: tk       !< temperature (K)
    real(kind_phys), intent(in) :: rh       !< relative humidity over liquid (fraction), for akel
    logical,         intent(in) :: kinetic  !< include kinetic fractionation
    logical :: liq_dst, liq_src

    wiso_alpha_phase = 1._kind_phys
    if (src == dst) return
    liq_dst = (dst == ph_ql .or. dst == ph_qr)
    liq_src = (src == ph_ql .or. src == ph_qr)
    if (src == ph_qv) then
       if (liq_dst) then
          wiso_alpha_phase = wiso_alpl(isp, tk)
          if (kinetic) wiso_alpha_phase = wiso_akel(isp, rh, wiso_alpha_phase)
       else
          wiso_alpha_phase = wiso_alpi(isp, tk)
          if (kinetic) wiso_alpha_phase = wiso_akci(isp, tk, wiso_alpha_phase)
       end if
    else if (dst == ph_qv) then
       if (liq_src) then
          wiso_alpha_phase = wiso_alpl(isp, tk)
          if (kinetic) wiso_alpha_phase = wiso_akel(isp, rh, wiso_alpha_phase)
          wiso_alpha_phase = 1._kind_phys / wiso_alpha_phase
       end if
    end if
  end function wiso_alpha_phase

  ! =====================================================================
  ! Kinetic fractionation
  ! =====================================================================

  !> Effective humidity near an evaporating drop (Jouzel type).
  elemental real(kind_phys) function wiso_heff(h0)
    real(kind_phys), intent(in) :: h0  !< relative humidity (fraction)
    wiso_heff = min(1.0_kind_phys, fkhum*h0 + 1.0_kind_phys - fkhum)
  end function wiso_heff

  !> Supersaturation over ice as a function of temperature (Hoffmann), 1 <= s <= ssatmx.
  elemental real(kind_phys) function wiso_ssatf(tk)
    real(kind_phys), intent(in) :: tk  !< temperature (K)
    wiso_ssatf = fsata + fsatb*(tk - tzero)
    wiso_ssatf = max(wiso_ssatf, 1.0_kind_phys)
    wiso_ssatf = min(wiso_ssatf, ssatmx)
  end function wiso_ssatf

  !> Effective fractionation for liquid evaporation into unsaturated air (kinetic effect).
  !! (iCESM wiso_akel; its temperature argument was not used and is dropped.)
  elemental real(kind_phys) function wiso_akel(isp, hum0, alpeq)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: hum0   !< relative humidity (fraction)
    real(kind_phys), intent(in) :: alpeq  !< equilibrium fractionation factor
    real(kind_phys) :: heff, dondi
    heff  = wiso_heff(min(1.0_kind_phys, hum0))
    dondi = (1._kind_phys/difrm(isp))**dkfac
    wiso_akel = alpeq*heff / (alpeq*dondi*(heff - 1._kind_phys) + 1._kind_phys)
  end function wiso_akel

  !> Effective fractionation for vapour deposition to ice (kinetic effect from supersaturation),
  !! below tkini; alpeq above.
  elemental real(kind_phys) function wiso_akci(isp, tk, alpeq)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: tk     !< temperature (K)
    real(kind_phys), intent(in) :: alpeq  !< equilibrium fractionation factor
    real(kind_phys) :: sat1, dondi
    if (tk < tkini) then
       sat1  = max(1._kind_phys, wiso_ssatf(tk))
       dondi = 1._kind_phys/difrm(isp)
       wiso_akci = alpeq*sat1 / (alpeq*dondi*(sat1 - 1._kind_phys) + 1._kind_phys)
    else
       wiso_akci = alpeq
    end if
  end function wiso_akci

  !> Kinetic factor (1 - k_mol) for ocean evaporation, Merlivat & Jouzel (1979), solving the
  !! Brutsaert equations for smooth/rough flow. Roughness from Charnock.
  elemental real(kind_phys) function wiso_kmol(isp, rbot, zbot, ustar)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: rbot   !< air density of the lowest layer (kg m-3)
    real(kind_phys), intent(in) :: zbot   !< height of the lowest level (m)
    real(kind_phys), intent(in) :: ustar  !< friction velocity (m s-1)
    real(kind_phys), parameter  :: difair = 2.36e-5_kind_phys  !< molecular diffusivity of air
    real(kind_phys), parameter  :: muair  = 1.7e-5_kind_phys   !< dynamic viscosity of air
    real(kind_phys) :: z0, vmu, sc, reno, tmr, enn, difn, kmol

    z0   = ustar**2 / (81.1_kind_phys*con_g)          ! Charnock
    vmu  = muair / rbot                               ! kinematic viscosity
    sc   = vmu / difair                               ! Schmidt number
    reno = ustar*z0 / vmu                             ! roughness Reynolds number
    if (reno < recrit) then                           ! smooth
       enn = 2._kind_phys/3._kind_phys
       tmr = ((1._kind_phys/karman)*log(ustar*zbot/(30._kind_phys*vmu))) / (13.6_kind_phys*sc**(2._kind_phys/3._kind_phys))
    else                                              ! rough
       enn = 0.5_kind_phys
       tmr = ((1._kind_phys/karman)*log(zbot/z0) - 5._kind_phys) / (7.3_kind_phys*reno**0.25_kind_phys*sc**0.5_kind_phys)
    end if
    difn = (1._kind_phys/difrm(isp))**enn
    kmol = (difn - 1._kind_phys) / (difn + tmr)
    wiso_kmol = 1._kind_phys - kmol
  end function wiso_kmol

  !> Kinetic factor (1 - k_mol) for ocean evaporation from the 10 m wind,
  !! Merlivat & Jouzel (1979) empirical smooth (< 7 m/s) / rough regimes.
  elemental real(kind_phys) function wiso_kmolv10(isp, ustar)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: ustar  !< friction velocity (m s-1)
    real(kind_phys) :: z0, v10, kmol
    z0  = ustar**2 / (81.1_kind_phys*con_g)
    v10 = ustar*log(10._kind_phys/z0)/karman
    if (v10 < 7.0_kind_phys) then
       kmol = aksmc(isp)
    else
       kmol = akrfa(isp)*v10 + akrfb(isp)
    end if
    wiso_kmolv10 = 1._kind_phys - kmol
  end function wiso_kmolv10

  ! =====================================================================
  ! Surface (ocean) evaporation
  ! =====================================================================

  !> Ocean surface isotope ratio (default, no ocean isotope data): boce * Rstd.
  elemental real(kind_phys) function wiso_get_roce(isp)
    integer, intent(in) :: isp
    wiso_get_roce = boce(isp)*rstd(isp)
  end function wiso_get_roce

  !> Isotope flux from the ocean, bulk formula (Merlivat & Jouzel 1979, D. Noone version):
  !!   flux = rho u* (1 - k_mol) re (q_iso,bot - qs(Ts) Roce/alpha_l(Ts))
  !! CESM sign convention: positive downward, negative = evaporation (kg m-2 s-1).
  !! re is the moisture exchange coefficient of the bulk flux (E = rho u* re (q_bot - qs)).
  elemental real(kind_phys) function wiso_flxoce(isp, rbot, zbot, wtbot, ts, rocn, ustar, re, ssq)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: rbot   !< air density of the lowest layer (kg m-3)
    real(kind_phys), intent(in) :: zbot   !< height of the lowest level (m)
    real(kind_phys), intent(in) :: wtbot  !< isotope vapour of the lowest layer (kg kg-1)
    real(kind_phys), intent(in) :: ts     !< sea surface temperature (K)
    real(kind_phys), intent(in) :: rocn   !< ocean isotope ratio / Rstd (0: use the default)
    real(kind_phys), intent(in) :: ustar  !< friction velocity (m s-1)
    real(kind_phys), intent(in) :: re     !< moisture exchange coefficient
    real(kind_phys), intent(in) :: ssq    !< saturation specific humidity at ts (kg kg-1)
    real(kind_phys) :: alpha, alpkn, roce, delq, qstar, tau

    alpha = wiso_alpl(isp, ts)
    alpkn = wiso_kmol(isp, rbot, zbot, ustar)
    if (rocn == 0._kind_phys) then
       roce = wiso_get_roce(isp)
    else
       roce = rstd(isp)*rocn
    end if
    delq  = wtbot - ssq*roce/alpha
    qstar = re*delq
    tau   = rbot*ustar*ustar
    wiso_flxoce = tau*alpkn*qstar/ustar
  end function wiso_flxoce

  !> Isotope ratio of the evaporation flux, Craig-Gordon in ratio form (Schmidt 1999; in
  !! iCESM wiso_flxoce as the commented alternative). Multiply by the bulk evaporation:
  !!   R_E = (1 - k_mol) (Roce/alpha_l(Ts) - h Rv) / (1 - h),  h = q_bot/qs(Ts)
  !! Returns 0 for h >= 1 (no evaporation; dew is a separate vapour -> surface transfer).
  elemental real(kind_phys) function wiso_cg_ratio(isp, ts, rocn, h, rv, alpkn)
    integer,         intent(in) :: isp    !< species id
    real(kind_phys), intent(in) :: ts     !< surface temperature (K)
    real(kind_phys), intent(in) :: rocn   !< surface water isotope ratio / Rstd (0: ocean default)
    real(kind_phys), intent(in) :: h      !< relative humidity of the lowest layer w.r.t. qs(ts)
    real(kind_phys), intent(in) :: rv     !< isotope ratio of the lowest-layer vapour
    real(kind_phys), intent(in) :: alpkn  !< kinetic factor (1 - k_mol), wiso_kmol or wiso_kmolv10
    real(kind_phys) :: roce
    if (rocn == 0._kind_phys) then
       roce = wiso_get_roce(isp)
    else
       roce = rstd(isp)*rocn
    end if
    if (h < 1._kind_phys) then
       wiso_cg_ratio = alpkn*(roce/wiso_alpl(isp, ts) - h*rv) / (1._kind_phys - h)
    else
       wiso_cg_ratio = 0._kind_phys
    end if
  end function wiso_cg_ratio

  ! =====================================================================
  ! Rain evaporation and liquid/vapour equilibration
  ! =====================================================================

  !> Equilibration factor: fraction of the total isotope in the vapour when liquid and vapour
  !! are in equilibrium (R_l = alpha R_v): efac = 1/(alpha l/v + 1), 0 <= efac <= 1.
  elemental real(kind_phys) function wiso_efac(alpha, vapnew, liqnew)
    real(kind_phys), intent(in) :: alpha   !< fractionation factor
    real(kind_phys), intent(in) :: vapnew  !< total vapour
    real(kind_phys), intent(in) :: liqnew  !< total liquid
    real(kind_phys) :: vfrac, alov
    if (vapnew <= 0._kind_phys) then
       wiso_efac = 0._kind_phys
    else if (vapnew + liqnew <= 0._kind_phys) then
       wiso_efac = 1._kind_phys
    else
       vfrac = vapnew/(vapnew + liqnew)
       alov  = alpha*(1.0_kind_phys/vfrac - 1.0_kind_phys)    ! alpha * l/v
       wiso_efac = min(1._kind_phys, max(0._kind_phys, 1.0_kind_phys/(alov + 1.0_kind_phys)))
    end if
  end function wiso_efac

  !> Change of isotope vapour when vapour and liquid (partially, fraction feq0) equilibrate,
  !! with dv + dl = 0 and R_l = alpha R_v. Limited to the isotope available.
  elemental real(kind_phys) function wiso_dqequil(alpha, feq0, vtotnew, ltotnew, visoold, lisoold)
    real(kind_phys), intent(in) :: alpha    !< fractionation factor
    real(kind_phys), intent(in) :: feq0     !< fraction equilibrated (0-1)
    real(kind_phys), intent(in) :: vtotnew  !< new total vapour
    real(kind_phys), intent(in) :: ltotnew  !< new total liquid
    real(kind_phys), intent(in) :: visoold  !< old isotope vapour
    real(kind_phys), intent(in) :: lisoold  !< old isotope liquid
    real(kind_phys) :: qiso, vieql, vinof, dviso
    qiso  = visoold + lisoold
    vieql = qiso*wiso_efac(alpha, vtotnew, ltotnew)          ! fractionating
    vinof = qiso*wiso_efac(1.0_kind_phys, vtotnew, ltotnew)  ! non-fractionating
    dviso = feq0*vieql + (1._kind_phys - feq0)*vinof - visoold
    if (dviso < 0._kind_phys) then
       dviso = max(dviso, -visoold)
    else
       dviso = min(dviso,  lisoold)
    end if
    wiso_dqequil = dviso
  end function wiso_dqequil

  !> Equilibrate isotope vapour and liquid (partially, fraction feq0), keeping the totals.
  !! Inputs must be >= 0.
  elemental subroutine wiso_liqvap_equil(alpha, feq0, vaptot, liqtot, vapiso, liqiso, dliqiso)
    real(kind_phys), intent(in)    :: alpha    !< fractionation factor
    real(kind_phys), intent(in)    :: feq0     !< fraction equilibrated (0-1)
    real(kind_phys), intent(in)    :: vaptot   !< total vapour
    real(kind_phys), intent(in)    :: liqtot   !< total liquid
    real(kind_phys), intent(inout) :: vapiso   !< isotope vapour
    real(kind_phys), intent(inout) :: liqiso   !< isotope liquid
    real(kind_phys), intent(out)   :: dliqiso  !< change in isotope liquid
    real(kind_phys), parameter :: qtiny = 1.e-36_kind_phys

    dliqiso = 0._kind_phys
    if (vaptot + liqtot < qtiny) return              ! no water
    if (vapiso + liqiso < qtiny) return              ! no isotope
    if (liqtot < qtiny) then                         ! no liquid: all isotope to vapour
       dliqiso = -liqiso
    else if (vaptot < qtiny) then                    ! no vapour: all isotope to liquid
       dliqiso = vapiso
    else
       dliqiso = -wiso_dqequil(alpha, feq0, vaptot, liqtot, vapiso, liqiso)
    end if
    liqiso = liqiso + dliqiso
    vapiso = vapiso - dliqiso
  end subroutine wiso_liqvap_equil

  !> Rayleigh distillation of isolated vapour from vtotold to vtotnew (vtotnew <= vtotold):
  !!   visonew = visoold (vtotnew/vtotold)**alpha
  elemental subroutine wiso_vap_distil(alpha, vtotold, vtotnew, visoold, visonew, dvapiso)
    real(kind_phys), intent(in)  :: alpha    !< fractionation factor
    real(kind_phys), intent(in)  :: vtotold  !< initial total vapour
    real(kind_phys), intent(in)  :: vtotnew  !< final total vapour
    real(kind_phys), intent(in)  :: visoold  !< initial isotope vapour
    real(kind_phys), intent(out) :: visonew  !< final isotope vapour
    real(kind_phys), intent(out) :: dvapiso  !< change in isotope vapour
    real(kind_phys), parameter :: qtiny = 1.e-22_kind_phys
    if (vtotold < qtiny) then
       visonew = visoold
    else
       visonew = visoold*(vtotnew/vtotold)**alpha
    end if
    dvapiso = visonew - visoold
  end subroutine wiso_vap_distil

  !> Fraction of a falling drop that equilibrates isotopically while crossing a layer
  !! (Lee & Fung 2008 adjustment time, Pruppacher & Klett ventilation, Foote & du Toit fall speed):
  !!   fequil = 1 - exp(-t_exposure / t_adjust)
  !! (iCESM wtrc_equil_time; esat is an argument instead of a qsat call.)
  elemental real(kind_phys) function wiso_equil_time(tk, pres, esat, rdrop, zdel, alpha, difrm_i)
    real(kind_phys), intent(in) :: tk       !< temperature (K)
    real(kind_phys), intent(in) :: pres     !< pressure (Pa)
    real(kind_phys), intent(in) :: esat     !< saturation vapour pressure over liquid (Pa)
    real(kind_phys), intent(in) :: rdrop    !< drop radius (m)
    real(kind_phys), intent(in) :: zdel     !< layer thickness (m)
    real(kind_phys), intent(in) :: alpha    !< equilibrium fractionation factor
    real(kind_phys), intent(in) :: difrm_i  !< diffusivity ratio of the species (wiso_get_difrm)
    real(kind_phys), parameter  :: cddrop = 0.6_kind_phys   !< drop drag coefficient (Straka 2009)
    real(kind_phys) :: rhoa, difa, mu, vterm, re, sc, xx, fvent, tadjust, texposure

    rhoa  = pres/(con_rd*tk)
    difa  = 2.11e-5_kind_phys*difrm_i*(tk/273.15_kind_phys)**1.94_kind_phys*(101325._kind_phys/pres)
    mu    = 1.72e-5_kind_phys*((tk/273.0_kind_phys)**1.5_kind_phys)*393.0_kind_phys/(tk + 120.0_kind_phys)
    vterm = sqrt((4._kind_phys/3._kind_phys)*con_g*(2._kind_phys*rdrop)*rhowater/(cddrop*rhoa))
    re    = 2._kind_phys*rdrop*rhoa*vterm/mu
    sc    = mu/(rhoa*difa)
    xx    = re**0.5_kind_phys * sc**(1._kind_phys/3._kind_phys)
    if (xx >= 1.4_kind_phys) then
       fvent = 0.78_kind_phys + 0.308_kind_phys*xx
    else
       fvent = 1.0_kind_phys + 0.108_kind_phys*xx*xx
    end if
    tadjust   = alpha*(rdrop*rdrop)*rhowater*con_rv*tk/(3._kind_phys*fvent*difa*esat)
    texposure = zdel/vterm
    wiso_equil_time = min(1._kind_phys, max(0._kind_phys, 1._kind_phys - exp(-texposure/tadjust)))
  end function wiso_equil_time

  !> Mass-weighted mean raindrop radius (m) from the rain rate, Marshall-Palmer
  !! (mean diameter in mm = 4/(4.1 R**-0.21), R in mm/h), as in iCESM stewart_isoevap.
  elemental real(kind_phys) function wiso_drop_radius(rain_mmh)
    real(kind_phys), intent(in) :: rain_mmh  !< rain rate (mm h-1), > 0
    wiso_drop_radius = 2._kind_phys/(4.1_kind_phys*rain_mmh**(-0.21_kind_phys)) / 1000._kind_phys
  end function wiso_drop_radius

  !> Isotope ratio of rain after part of it evaporates while falling through a layer,
  !! Stewart (1975) with effective humidity heff = phi + (1 - phi) rh (Bony et al. 2008):
  !!   R_r = fequil R_Stewart + (1 - fequil) R_r0,  R_Stewart = (R_r0 - gam R_v0) fr**bet + gam R_v0
  !! Returns R_r0 when heff >= 1 or nothing evaporates (fr >= 1); for heff >= 1 iCESM then
  !! equilibrates rain and vapour with wiso_liqvap_equil (fraction fequil).
  elemental real(kind_phys) function wiso_stewart_ratio(isp, tk, rh, rr0, rv0, fr, fequil)
    integer,         intent(in) :: isp     !< species id
    real(kind_phys), intent(in) :: tk      !< air temperature (K)
    real(kind_phys), intent(in) :: rh      !< relative humidity over liquid (fraction), before evaporation
    real(kind_phys), intent(in) :: rr0     !< rain isotope ratio before evaporation
    real(kind_phys), intent(in) :: rv0     !< vapour isotope ratio before evaporation
    real(kind_phys), intent(in) :: fr      !< fraction of the rain mass remaining (0-1)
    real(kind_phys), intent(in) :: fequil  !< fraction that equilibrates (wiso_equil_time)
    real(kind_phys) :: heff, ae, ak, gam, bet, frc, rstw

    heff = phi_stw + (1._kind_phys - phi_stw)*rh
    frc  = min(1._kind_phys, max(0._kind_phys, fr))
    wiso_stewart_ratio = rr0
    if (heff < 1._kind_phys .and. frc < 1._kind_phys) then
       ae  = wiso_alpl(isp, tk)
       ak  = (1._kind_phys/difrm(isp))**dkfac
       gam = (ae*heff) / (1._kind_phys - ae*ak*(1._kind_phys - heff))
       bet = (1._kind_phys - ae*ak*(1._kind_phys - heff)) / (ae*ak*(1._kind_phys - heff))
       if (.not. (frc == 0._kind_phys .and. bet < 0._kind_phys)) then
          rstw = (rr0 - gam*rv0)*frc**bet + gam*rv0
          wiso_stewart_ratio = fequil*rstw + (1._kind_phys - fequil)*rr0
       end if
    end if
  end function wiso_stewart_ratio

  ! =====================================================================
  ! Ratios, delta and constants
  ! =====================================================================

  !> Isotope ratio qiso/qtot; Rstd when qtot <= 0.
  elemental real(kind_phys) function wiso_ratio(isp, qiso, qtot)
    integer,         intent(in) :: isp
    real(kind_phys), intent(in) :: qiso, qtot
    if (qtot > 0._kind_phys) then
       wiso_ratio = qiso/qtot
    else
       wiso_ratio = rstd(isp)
    end if
  end function wiso_ratio

  !> delta value (permil) relative to the standard: 1000 (R/Rstd - 1).
  elemental real(kind_phys) function wiso_delta(isp, qiso, qtot)
    integer,         intent(in) :: isp
    real(kind_phys), intent(in) :: qiso, qtot
    wiso_delta = 1000._kind_phys*(wiso_ratio(isp, qiso, qtot)/rstd(isp) - 1._kind_phys)
  end function wiso_delta

  !> Model standard ratio (1 for all species).
  elemental real(kind_phys) function wiso_get_rstd(isp)
    integer, intent(in) :: isp
    wiso_get_rstd = rstd(isp)
  end function wiso_get_rstd

  !> Natural (VSMOW) ratio.
  elemental real(kind_phys) function wiso_get_rnat(isp)
    integer, intent(in) :: isp
    wiso_get_rnat = rnat(isp)
  end function wiso_get_rnat

  !> Number of isotope atoms per molecule.
  elemental real(kind_phys) function wiso_get_fisub(isp)
    integer, intent(in) :: isp
    wiso_get_fisub = fisub(isp)
  end function wiso_get_fisub

  !> Diffusivity ratio Di/D (Merlivat 1978).
  elemental real(kind_phys) function wiso_get_difrm(isp)
    integer, intent(in) :: isp
    wiso_get_difrm = difrm(isp)
  end function wiso_get_difrm

end module water_isotopes
