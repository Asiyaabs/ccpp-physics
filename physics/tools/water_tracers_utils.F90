!> \file water_tracers_utils.F90
!! Water tracers (tagging, later isotopes): the one helper used inside the physics schemes.
!!
!! Plain Fortran module (no CCPP metadata). A scheme that moves water calls wtr_transfer once
!! per water transfer, for every water tracer set, with the same amount of model water dq:
!!
!!   moved tracer = alpha * (src_trc / src_q) * dq
!!
!! i.e. the tracer leaves the source store at the source ratio (tracer/model water), times
!! the fractionation factor alpha of that process (1 for H2O and tags; wtr_alpha below).
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

  use machine, only: kind_phys

  implicit none
  private

  public :: wtr_transfer, wtr_alpha
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
  elemental subroutine wtr_transfer(src_q, src_trc, dq, dst_trc, alpha)
    real(kind_phys), intent(in)              :: src_q    !< model water in the source store, before the transfer
    real(kind_phys), intent(inout)           :: src_trc  !< tracer in the source store
    real(kind_phys), intent(in)              :: dq       !< model water moved out of the source (>= 0)
    real(kind_phys), intent(inout), optional :: dst_trc  !< tracer in the destination store
    real(kind_phys), intent(in),    optional :: alpha    !< fractionation factor (default 1)

    real(kind_phys) :: dtrc

    if (dq <= 0.0_kind_phys .or. src_q <= 0.0_kind_phys .or. src_trc <= 0.0_kind_phys) return

    dtrc = min(dq / src_q, 1.0_kind_phys) * src_trc
    if (present(alpha)) dtrc = min(alpha * dtrc, src_trc)

    src_trc = src_trc - dtrc
    if (present(dst_trc)) dst_trc = dst_trc + dtrc
  end subroutine wtr_transfer

  !> Fractionation factor of a process for a species, at temperature temp (K).
  !! 1 for H2O (bulk and tags) and, for now, for all species and processes.
  !! The isotope values (HDO, H218O) are filled in Phase 4.
  elemental real(kind_phys) function wtr_alpha(species, process, temp)
    integer,         intent(in) :: species  !< species id (water_tracers_config: 1 H2O, 2 HDO, 3 H218O)
    integer,         intent(in) :: process  !< process id (wtr_proc_*)
    real(kind_phys), intent(in) :: temp     !< temperature (K)

    wtr_alpha = 1.0_kind_phys
  end function wtr_alpha

end module water_tracers_utils
