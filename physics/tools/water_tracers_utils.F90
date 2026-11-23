!> \file water_tracers_utils.F90
!! Water tracers (tagging, isotopes): the one helper used inside the physics schemes.
!!
!! Plain Fortran module (no CCPP metadata). A scheme that moves water calls wtr_transfer once
!! per water transfer, for every water tracer set, with the same amount of model water dq:
!!
!!   moved tracer = alpha * (src_trc / src_q) * dq          (alpha = 1: exactly this)
!!
!! i.e. the tracer leaves the source store at the source ratio (tracer/model water), times
!! the fractionation factor alpha of that process (1 for H2O and tags; wtr_alpha below).
!! For alpha /= 1 (isotopes) the moved water is in equilibrium with the water left in the
!! source (closed-system step), which is the formula above for small dq/src_q and never
!! moves more tracer than the source holds.
!!
!! Usage, e.g. condensation of vapour to cloud water (dq > 0) for set s:
!!
!!   call wtr_transfer(qv, q(:,k,wtr_idx(s,1)), dq, q(:,k,wtr_idx(s,2)), &
!!                     wtr_alpha(wtr_species(s), wtr_proc_cond_liq, t))
!!
!! Rules:
!!  - src_q is the model water of the source store BEFORE this transfer.
!!  - dq >= 0 is the model water leaving the source. For the opposite direction, swap the
!!    stores (e.g. evaporation of cloud water: source ql, destination qv).
!!  - dst_trc absent: the tracer leaves the atmosphere (e.g. precipitation reaching the surface).
!!  - Transfers that leave the same source store without fractionation can be done in one
!!    call with the summed dq; every fractionating process gets its own call.
module water_tracers_utils

  use machine,        only: kind_phys
  use water_isotopes, only: wiso_alpl, wiso_alpi

  implicit none
  private

  public :: wtr_transfer, wtr_alpha, wtr_exchange, wtr_net, wtr_sediment, wtr_ratio_limit, wtr_cnv_substep, &
            wtr_col_ratio
  public :: wtr_proc_none, wtr_proc_cond_liq, wtr_proc_dep_ice, wtr_proc_evap_liq, &
            wtr_proc_subl_ice, wtr_proc_freeze, wtr_proc_melt, wtr_proc_sfc_evap,   &
            wtr_proc_seaice_growth

  !> Process ids for wtr_alpha (where fractionation can happen)
  integer, parameter :: wtr_proc_none          = 0  !< no fractionation (mixing, fall-out, runoff, ...)
  integer, parameter :: wtr_proc_cond_liq      = 1  !< vapour -> liquid (equilibrium)
  integer, parameter :: wtr_proc_dep_ice       = 2  !< vapour -> ice (equilibrium + kinetic)
  integer, parameter :: wtr_proc_evap_liq      = 3  !< liquid -> vapour (rain evaporation, partial)
  integer, parameter :: wtr_proc_subl_ice      = 4  !< ice -> vapour
  integer, parameter :: wtr_proc_freeze        = 5  !< liquid -> ice
  integer, parameter :: wtr_proc_melt          = 6  !< ice -> liquid
  integer, parameter :: wtr_proc_sfc_evap      = 7  !< surface -> vapour (Craig-Gordon)
  integer, parameter :: wtr_proc_seaice_growth = 8  !< seawater -> sea ice

contains

  !> Move the tracer share of a model-water transfer dq from a source store to a
  !! destination store (or out of the atmosphere if dst_trc is absent).
  elemental subroutine wtr_transfer(src_q, src_trc, dq, dst_trc, alpha, rdef)
    real(kind_phys), intent(in)              :: src_q    !< model water in the source store, before the transfer
    real(kind_phys), intent(inout)           :: src_trc  !< tracer in the source store
    real(kind_phys), intent(in)              :: dq       !< model water moved out of the source (>= 0)
    real(kind_phys), intent(inout), optional :: dst_trc  !< tracer in the destination store
    real(kind_phys), intent(in),    optional :: alpha    !< fractionation factor (default 1)
    real(kind_phys), intent(in),    optional :: rdef     !< tracer/model ratio of water the model takes
                                                         !! beyond what the source holds (wtr_col_ratio)

    real(kind_phys) :: dtrc, f

    if (dq <= 0.0_kind_phys) return

    if (present(rdef) .and. dq > max(src_q, 0.0_kind_phys)) then
      ! The model takes more than the source holds (e.g. a cloud top where the parcel water
      ! or the saturation value is < 0), i.e. it creates the part beyond max(src_q, 0):
      ! the source tracer gives all it holds (an emptied store: no fractionation), the
      ! created part gets the ratio rdef.
      dtrc = (dq - max(src_q, 0.0_kind_phys)) * rdef
      if (src_q > 0.0_kind_phys .and. src_trc > 0.0_kind_phys) dtrc = dtrc + src_trc
    else
      if (src_q <= 0.0_kind_phys .or. src_trc <= 0.0_kind_phys) return
      ! Default (alpha absent or 1: H2O, tags), and wherever the model takes all of a store
      ! (f = 1) or more than it holds (f > 1): the tracer follows at the source ratio, not
      ! capped (a bulk set mirrors the model where it drives a store negative; tags are
      ! clipped by the host; an emptied store has nothing left to fractionate against).
      f = dq / src_q
      dtrc = f * src_trc
      if (present(alpha)) then
        if (alpha /= 1.0_kind_phys .and. f < 1.0_kind_phys) then
          ! Fractionation (isotopes): the water moved is in equilibrium with the water left
          ! in the source, (dtrc/dq) / ((src_trc-dtrc)/(src_q-dq)) = alpha. Equal to
          ! alpha*f*src_trc for small f, and never more than the source holds (f -> 1: all).
          dtrc = src_trc * alpha * f / (1.0_kind_phys + (alpha - 1.0_kind_phys) * f)
        end if
      end if
    end if

    src_trc = src_trc - dtrc
    if (present(dst_trc)) dst_trc = dst_trc + dtrc
  end subroutine wtr_transfer

  !> Two-way exchange between stores a and b of model water dq (> 0: a -> b, < 0: b -> a),
  !! with the fractionation factor of each direction (one fractionating process).
  elemental subroutine wtr_exchange(qa, ta, qb, tb, dq, alpha_ab, alpha_ba)
    real(kind_phys), intent(in)    :: qa, qb             !< model water in a and b, before the exchange
    real(kind_phys), intent(inout) :: ta, tb             !< tracer in a and b
    real(kind_phys), intent(in)    :: dq                 !< model water moved from a to b (signed)
    real(kind_phys), intent(in)    :: alpha_ab, alpha_ba !< fractionation factors a -> b and b -> a

    if (dq > 0.0_kind_phys) then
      call wtr_transfer(qa, ta, dq, tb, alpha_ab)
    else if (dq < 0.0_kind_phys) then
      call wtr_transfer(qb, tb, -dq, ta, alpha_ba)
    end if
  end subroutine wtr_exchange

  !> Net change of n stores at one point after a list of non-fractionating steps.
  !! Stores that lose water give tracer at their own ratio into a common pool; the gains,
  !! up to the water lost, share the pool. Gains beyond the water lost (water created by
  !! the model, e.g. a floor) get the ratio rdef; pool water not regained is removed.
  pure subroutine wtr_net(n, q0, dq, trc, rdef)
    integer,         intent(in)    :: n
    real(kind_phys), intent(in)    :: q0(n)   !< model water before the steps
    real(kind_phys), intent(in)    :: dq(n)   !< model net change over the steps
    real(kind_phys), intent(inout) :: trc(n)  !< tracer
    real(kind_phys), intent(in)    :: rdef    !< tracer/model ratio for created water
    real(kind_phys) :: lost, gain, pool, fpool
    integer :: p

    lost = 0.0_kind_phys
    gain = 0.0_kind_phys
    pool = 0.0_kind_phys
    do p = 1, n
      if (dq(p) < 0.0_kind_phys) then
        lost = lost - dq(p)
        call wtr_transfer(q0(p), trc(p), -dq(p), pool)
      else
        gain = gain + dq(p)
      end if
    end do
    if (gain <= 0.0_kind_phys) return
    fpool = min(lost, gain) / gain        ! part of each gain that comes from the pool
    do p = 1, n
      if (dq(p) > 0.0_kind_phys) then
        trc(p) = trc(p) + dq(p) * (1.0_kind_phys - fpool) * rdef
        if (lost > 0.0_kind_phys) trc(p) = trc(p) + dq(p) * fpool * pool / lost
      end if
    end do
  end subroutine wtr_net

  !> Net sedimentation of one species in a column (index 1 = bottom, nk = top), from the
  !! model profiles before and after (all sub-steps together, no fractionation). The flux
  !! out of the bottom of each level is rebuilt from the model mass budget, top down, and
  !! carries the tracer ratio of the level's water plus the inflow from above.
  pure subroutine wtr_sediment(nk, w, q0, q1, trc, sfc)
    integer,         intent(in)    :: nk
    real(kind_phys), intent(in)    :: w(nk)   !< mass weight of each level (e.g. rho*dz)
    real(kind_phys), intent(in)    :: q0(nk)  !< model mixing ratio before sedimentation
    real(kind_phys), intent(in)    :: q1(nk)  !< model mixing ratio after sedimentation
    real(kind_phys), intent(inout) :: trc(nk) !< tracer mixing ratio
    real(kind_phys), intent(out)   :: sfc     !< tracer mass reaching the surface (weight units)
    real(kind_phys) :: fin, tin, mix, tmix, fout
    integer :: k

    fin = 0.0_kind_phys
    tin = 0.0_kind_phys
    do k = nk, 1, -1
      mix  = q0(k)*w(k) + fin
      tmix = trc(k)*w(k) + tin
      fout = mix - q1(k)*w(k)
      tin  = 0.0_kind_phys
      call wtr_transfer(mix, tmix, fout, tin)
      trc(k) = tmix / w(k)
      fin = max(fout, 0.0_kind_phys)
    end do
    sfc = tin
  end subroutine wtr_sediment

  !> Keep the tracer/model ratio of a column after a scheme within the range it had before
  !! the scheme (levels with model water > 0). For explicit transports that can move more
  !! than a layer's mass per step: exact for a bulk set (ratio 1), bounded for tags.
  pure subroutine wtr_ratio_limit(n, q0, t0, q1, t1)
    integer,         intent(in)    :: n
    real(kind_phys), intent(in)    :: q0(n), t0(n)  !< model water and tracer before the scheme
    real(kind_phys), intent(in)    :: q1(n)         !< model water after the scheme
    real(kind_phys), intent(inout) :: t1(n)         !< tracer after the scheme
    real(kind_phys) :: rmin, rmax
    integer :: k

    rmin =  huge(1.0_kind_phys)
    rmax = -huge(1.0_kind_phys)
    do k = 1, n
      if (q0(k) > 0.0_kind_phys) then
        rmin = min(rmin, t0(k)/q0(k))
        rmax = max(rmax, t0(k)/q0(k))
      end if
    end do
    if (rmax < rmin) return
    do k = 1, n
      if (q1(k) > 0.0_kind_phys) t1(k) = q1(k) * min(rmax, max(rmin, t1(k)/q1(k)))
    end do
  end subroutine wtr_ratio_limit

  !> Convective (SAMF) transport of one vapour field in a column, sub-stepped: the
  !! scheme's own source terms s (dellae before the environment term: updraft and
  !! downdraft fluxes, entrainment, detrainment) and the environment subsidence/uplift
  !! by the net cloud mass flux mf, first-order upwind, in ns sub-steps of at most one
  !! layer each (the scheme takes one explicit step of up to several layers, which is
  !! unstable for ratios: isotope compositions then alternate level by level). Same
  !! level ranges and kb/jmin double-counting corrections as the scheme. Used for
  !! isotope sets: ratio = (isotope transported) / (model vapour transported).
  pure subroutine wtr_cnv_substep(nk, ktop, kb, jmin, mf, ekb, djm, dp, grav, xdt, s, x)
    integer,         intent(in)    :: nk, ktop, kb, jmin  !< levels, cloud top, cloud base, downdraft start (0: none)
    real(kind_phys), intent(in)    :: mf(nk)    !< net cloud mass flux at the top of each layer, per unit xmb (> 0: environment subsides)
    real(kind_phys), intent(in)    :: ekb, djm  !< eta(kb) and edto*etad(jmin) of the scheme's double-counting corrections
    real(kind_phys), intent(in)    :: dp(nk)    !< layer pressure thickness (Pa)
    real(kind_phys), intent(in)    :: grav, xdt !< gravity; cloud-base mass flux times time step
    real(kind_phys), intent(in)    :: s(nk)     !< source tendency per unit xmb (scheme's dellae without the environment term)
    real(kind_phys), intent(inout) :: x(nk)     !< field before (in) and after (out) the transport
    real(kind_phys) :: flx(0:nk), t(nk), c, xi
    integer :: ns, it, k

    c = 0.0_kind_phys
    do k = 1, ktop-1
      c = max(c, abs(mf(k))*xdt*grav/min(dp(k), dp(k+1)))
    end do
    ns = max(1, ceiling(c))
    do it = 1, ns
      flx = 0.0_kind_phys
      t   = s
      do k = 1, ktop-1
        if (mf(k) > 0.0_kind_phys) then
          xi = x(k+1)
        else if (mf(k) < 0.0_kind_phys) then
          xi = x(k)
        else
          xi = 0.5_kind_phys*(x(k) + x(k+1))
        end if
        flx(k) = mf(k)*xi
        if (k == jmin) t(k+1) = t(k+1) - djm*xi*grav/dp(k+1)
        if (k == kb)   t(k)   = t(k)   - ekb*xi*grav/dp(k)
      end do
      do k = 1, ktop
        x(k) = x(k) + (t(k) + (flx(k) - flx(k-1))*grav/dp(k))*xdt/ns
      end do
    end do
  end subroutine wtr_cnv_substep

  !> Tracer/model ratio of a column (levels with model water > 0), for water the model
  !! creates where there is no ratio to follow (rdef of wtr_transfer). rdef here: the value
  !! when no level has water (1 for a bulk set, 0 for a tag).
  pure real(kind_phys) function wtr_col_ratio(q, t, rdef)
    real(kind_phys), intent(in) :: q(:), t(:)  !< model water and tracer of the column
    real(kind_phys), intent(in) :: rdef

    if (any(q > 0.0_kind_phys)) then
      wtr_col_ratio = sum(t, mask = q > 0.0_kind_phys) / sum(q, mask = q > 0.0_kind_phys)
    else
      wtr_col_ratio = rdef
    end if
  end function wtr_col_ratio

  !> Fractionation factor of a process for a species, at temperature temp (K): the alpha
  !! of wtr_transfer for water leaving the source store of that process.
  !! Equilibrium fractionation only (isotopes step 2a), as iCESM (water_isotopes):
  !!   vapour -> liquid (cond_liq):   alpha_l(T)    Horita & Wesolowski (1994)
  !!   vapour -> ice    (dep_ice):    alpha_i(T)    Merlivat & Nief (1967), Majoube (1971b)
  !!   liquid -> vapour (evap_liq):   1/alpha_l(T)
  !!   ice -> vapour    (subl_ice):   1             (ice sublimates layer by layer)
  !!   freezing, melting, surface evaporation, sea-ice growth, none: 1 (for now)
  !! Exactly 1 for H2O (bulk and tags). Kinetic effects (ice supersaturation, evaporation
  !! into unsaturated air) and the surface values (Craig-Gordon) come in later steps.
  elemental real(kind_phys) function wtr_alpha(species, process, temp)
    integer,         intent(in) :: species  !< species id (water_tracers_config: 1 H2O, 2 HDO, 3 H218O)
    integer,         intent(in) :: process  !< process id (wtr_proc_*)
    real(kind_phys), intent(in) :: temp     !< temperature (K)

    wtr_alpha = 1.0_kind_phys
    if (species <= 1) return
    select case (process)
    case (wtr_proc_cond_liq)
      wtr_alpha = wiso_alpl(species, temp)
    case (wtr_proc_dep_ice)
      wtr_alpha = wiso_alpi(species, temp)
    case (wtr_proc_evap_liq)
      wtr_alpha = 1.0_kind_phys / wiso_alpl(species, temp)
    end select
  end function wtr_alpha

end module water_tracers_utils
