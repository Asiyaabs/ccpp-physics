!> \file water_tracers_config.F90
!! Water tracers (tagging, later isotopes): namelist, tracer names and index tables.
!!
!! Plain Fortran module (no CCPP metadata), used by the host (to add the tracers to the
!! tracer list) and by GFS_typedefs (to build the index tables passed to the schemes).
!!
!! Namelist group, in the physics namelist file (optional; absent = tagging off):
!!
!!   &water_tracer_nml
!!     tagging        = .false.   ! master switch: .true. adds the water tracers
!!     wtr_species    = 'H2O'     ! species to track ('H2O'; 'HDO', 'H218O' later)
!!     wtr_tag_names  = 'COL'     ! one tag set per name (SCM: COL = the column)
!!     wtr_tag_source = 'surface' ! where tagged water enters (per tag)
!!     wtr_bulk_check = .true.    ! H2O bulk set, must equal the model water exactly
!!     wtr_init_tags  = 'zero'    ! initial tags: 'zero' | 'bulk' | 'fraction'
!!     wtr_init_frac  = 0.5       ! initial tag fraction for wtr_init_tags = 'fraction'
!!     wtr_tag_moisten = .false.  ! test mode: tags also get the forcing moistening (with
!!                                ! wtr_init_tags = 'bulk', every tag must then equal the bulk set)
!!   /
!!
!! Each set holds one tracer per water phase (qv, ql, qi, qr, qs, qg). Sets, in order:
!! for each species: its bulk set (for H2O only if wtr_bulk_check), then one tag set per tag.
!! Tracer names: <phase>_<species>[_<tag>], e.g. qv_H2O, qv_H2O_COL.
module water_tracers_config

  use machine, only: kind_phys

  implicit none
  private

  public :: wtr_config_type, wtr_read_nml, wtr_nset, wtr_set_info, wtr_tracer_names, wtr_index_table
  public :: wtr_nphase, wtr_phase_names, wtr_max_species, wtr_max_tags
  public :: wtr_species_h2o, wtr_species_hdo, wtr_species_h218o, wtr_species_names
  public :: wtr_kind_bulk, wtr_kind_tag

  integer,          parameter :: wtr_nphase = 6                          !< phases per set
  character(len=2), parameter :: wtr_phase_names(wtr_nphase) = (/ 'qv', 'ql', 'qi', 'qr', 'qs', 'qg' /)

  integer,          parameter :: wtr_max_species   = 3
  integer,          parameter :: wtr_species_h2o   = 1
  integer,          parameter :: wtr_species_hdo   = 2
  integer,          parameter :: wtr_species_h218o = 3
  character(len=5), parameter :: wtr_species_names(wtr_max_species) = (/ 'H2O  ', 'HDO  ', 'H218O' /)

  integer,          parameter :: wtr_max_tags  = 10
  integer,          parameter :: wtr_kind_bulk = 1
  integer,          parameter :: wtr_kind_tag  = 2

  !> User choices from the namelist, checked
  type wtr_config_type
     logical            :: tagging    = .false.
     integer            :: nspecies   = 0
     integer            :: species(wtr_max_species) = 0          !< species ids, in namelist order
     integer            :: ntags      = 0
     character(len=16)  :: tag_names(wtr_max_tags)  = ' '
     character(len=16)  :: tag_source(wtr_max_tags) = ' '
     logical            :: bulk_check = .true.
     character(len=16)  :: init_tags  = 'zero'
     real(kind_phys)    :: init_frac  = 0.5_kind_phys
     logical            :: tag_moisten = .false.
  end type wtr_config_type

contains

  !> Read &water_tracer_nml from a namelist file (unit opened here) or, if lines is
  !! present, from an internal namelist file. A missing group means tagging off.
  subroutine wtr_read_nml(cfg, errmsg, errflg, filename, lines)
    type(wtr_config_type),                      intent(out) :: cfg
    character(len=*),                           intent(out) :: errmsg
    integer,                                    intent(out) :: errflg
    character(len=*),                 optional, intent(in)  :: filename
    character(len=*), dimension(:),   optional, intent(in)  :: lines

    logical            :: tagging
    character(len=8)   :: wtr_species(wtr_max_species)
    character(len=16)  :: wtr_tag_names(wtr_max_tags), wtr_tag_source(wtr_max_tags)
    logical            :: wtr_bulk_check
    character(len=16)  :: wtr_init_tags
    real(kind_phys)    :: wtr_init_frac
    logical            :: wtr_tag_moisten
    integer            :: unit, ios

    namelist /water_tracer_nml/ tagging, wtr_species, wtr_tag_names, wtr_tag_source, &
                                wtr_bulk_check, wtr_init_tags, wtr_init_frac, wtr_tag_moisten

    errmsg = ''
    errflg = 0

    tagging           = .false.
    wtr_species       = ' '
    wtr_species(1)    = 'H2O'
    wtr_tag_names     = ' '
    wtr_tag_names(1)  = 'COL'
    wtr_tag_source    = 'surface'
    wtr_bulk_check    = .true.
    wtr_init_tags     = 'zero'
    wtr_init_frac     = 0.5_kind_phys
    wtr_tag_moisten   = .false.

    if (present(lines)) then
       read(lines, nml=water_tracer_nml, iostat=ios)
    else if (present(filename)) then
       open(newunit=unit, file=trim(filename), action='read', status='old', iostat=ios)
       if (ios /= 0) then
          errmsg = 'water_tracers_config: cannot open namelist file ' // trim(filename)
          errflg = 1
          return
       end if
       read(unit, nml=water_tracer_nml, iostat=ios)
       close(unit)
    else
       errmsg = 'water_tracers_config: wtr_read_nml needs filename or lines'
       errflg = 1
       return
    end if

    if (ios < 0) then
       ! group not present: tagging off, nothing else to check
       cfg%tagging = .false.
       return
    else if (ios > 0) then
       write(errmsg,'(a,i0)') 'water_tracers_config: error reading &water_tracer_nml, iostat = ', ios
       errflg = 1
       return
    end if

    cfg%tagging = tagging
    if (.not. cfg%tagging) return

    call set_species(cfg, wtr_species, errmsg, errflg)
    if (errflg /= 0) return
    call set_tags(cfg, wtr_tag_names, wtr_tag_source, errmsg, errflg)
    if (errflg /= 0) return

    cfg%bulk_check = wtr_bulk_check
    cfg%init_tags  = adjustl(wtr_init_tags)
    cfg%init_frac  = wtr_init_frac
    cfg%tag_moisten = wtr_tag_moisten

    select case (trim(cfg%init_tags))
    case ('zero', 'bulk')
    case ('fraction')
       if (cfg%init_frac < 0.0_kind_phys .or. cfg%init_frac > 1.0_kind_phys) then
          errmsg = 'water_tracers_config: wtr_init_frac must be between 0 and 1'
          errflg = 1
          return
       end if
    case default
       errmsg = 'water_tracers_config: wtr_init_tags must be zero, bulk or fraction, got ' // trim(cfg%init_tags)
       errflg = 1
       return
    end select

    if (wtr_nset(cfg) == 0) then
       errmsg = 'water_tracers_config: tagging = .true. but no tracer set (no tags and wtr_bulk_check = .false.)'
       errflg = 1
       return
    end if

  end subroutine wtr_read_nml

  subroutine set_species(cfg, names, errmsg, errflg)
    type(wtr_config_type), intent(inout) :: cfg
    character(len=*),      intent(in)    :: names(:)
    character(len=*),      intent(out)   :: errmsg
    integer,               intent(out)   :: errflg
    integer :: n, m, id

    errmsg = ''
    errflg = 0
    cfg%nspecies = 0
    do n = 1, size(names)
       if (len_trim(names(n)) == 0) cycle
       id = 0
       do m = 1, wtr_max_species
          if (trim(adjustl(names(n))) == trim(wtr_species_names(m))) id = m
       end do
       if (id == 0) then
          errmsg = 'water_tracers_config: unknown species in wtr_species: ' // trim(names(n))
          errflg = 1
          return
       end if
       if (id /= wtr_species_h2o) then
          errmsg = 'water_tracers_config: species ' // trim(names(n)) // ' not implemented yet (only H2O)'
          errflg = 1
          return
       end if
       if (any(cfg%species(1:cfg%nspecies) == id)) then
          errmsg = 'water_tracers_config: species listed twice in wtr_species: ' // trim(names(n))
          errflg = 1
          return
       end if
       cfg%nspecies = cfg%nspecies + 1
       cfg%species(cfg%nspecies) = id
    end do
    if (cfg%nspecies == 0) then
       errmsg = 'water_tracers_config: tagging = .true. but wtr_species is empty'
       errflg = 1
    end if
  end subroutine set_species

  subroutine set_tags(cfg, names, sources, errmsg, errflg)
    type(wtr_config_type), intent(inout) :: cfg
    character(len=*),      intent(in)    :: names(:), sources(:)
    character(len=*),      intent(out)   :: errmsg
    integer,               intent(out)   :: errflg
    integer :: n

    errmsg = ''
    errflg = 0
    cfg%ntags = 0
    do n = 1, size(names)
       if (len_trim(names(n)) == 0) cycle
       if (any(cfg%tag_names(1:cfg%ntags) == adjustl(names(n)))) then
          errmsg = 'water_tracers_config: tag listed twice in wtr_tag_names: ' // trim(names(n))
          errflg = 1
          return
       end if
       if (trim(adjustl(sources(n))) /= 'surface') then
          errmsg = 'water_tracers_config: wtr_tag_source ' // trim(sources(n)) // &
                   ' not implemented yet (only surface), tag ' // trim(names(n))
          errflg = 1
          return
       end if
       cfg%ntags = cfg%ntags + 1
       cfg%tag_names(cfg%ntags)  = adjustl(names(n))
       cfg%tag_source(cfg%ntags) = adjustl(sources(n))
    end do
  end subroutine set_tags

  !> Number of tracer sets
  pure integer function wtr_nset(cfg)
    type(wtr_config_type), intent(in) :: cfg
    integer :: s
    wtr_nset = 0
    if (.not. cfg%tagging) return
    do s = 1, cfg%nspecies
       if (has_bulk_set(cfg, cfg%species(s))) wtr_nset = wtr_nset + 1
       wtr_nset = wtr_nset + cfg%ntags
    end do
  end function wtr_nset

  pure logical function has_bulk_set(cfg, species)
    type(wtr_config_type), intent(in) :: cfg
    integer,               intent(in) :: species
    ! isotopes always need their bulk set; for H2O it is the optional check set
    has_bulk_set = (species /= wtr_species_h2o) .or. cfg%bulk_check
  end function has_bulk_set

  !> Per set: name (e.g. H2O, H2O_COL), species id, kind (bulk/tag), and the set index of
  !! the bulk set of the same species (0: the model water itself, for H2O tags without check set)
  subroutine wtr_set_info(cfg, set_name, set_species, set_kind, set_bulk)
    type(wtr_config_type), intent(in)  :: cfg
    character(len=*),      intent(out) :: set_name(:)
    integer,               intent(out) :: set_species(:), set_kind(:), set_bulk(:)
    integer :: s, t, iset, ibulk

    iset = 0
    do s = 1, cfg%nspecies
       ibulk = 0
       if (has_bulk_set(cfg, cfg%species(s))) then
          iset = iset + 1
          ibulk = iset
          set_name(iset)    = trim(wtr_species_names(cfg%species(s)))
          set_species(iset) = cfg%species(s)
          set_kind(iset)    = wtr_kind_bulk
          set_bulk(iset)    = iset
       end if
       do t = 1, cfg%ntags
          iset = iset + 1
          set_name(iset)    = trim(wtr_species_names(cfg%species(s))) // '_' // trim(cfg%tag_names(t))
          set_species(iset) = cfg%species(s)
          set_kind(iset)    = wtr_kind_tag
          set_bulk(iset)    = ibulk
       end do
    end do
  end subroutine wtr_set_info

  !> Tracer names to add to the tracer list, set by set, phases in wtr_phase_names order
  subroutine wtr_tracer_names(cfg, names)
    type(wtr_config_type), intent(in)  :: cfg
    character(len=*),      intent(out) :: names(:)   !< size wtr_nset(cfg) * wtr_nphase
    character(len=32) :: set_name(wtr_nset(cfg))
    integer           :: set_species(wtr_nset(cfg)), set_kind(wtr_nset(cfg)), set_bulk(wtr_nset(cfg))
    integer           :: iset, ip

    call wtr_set_info(cfg, set_name, set_species, set_kind, set_bulk)
    do iset = 1, wtr_nset(cfg)
       do ip = 1, wtr_nphase
          names((iset-1)*wtr_nphase + ip) = wtr_phase_names(ip) // '_' // trim(set_name(iset))
       end do
    end do
  end subroutine wtr_tracer_names

  !> Index table: idx(iset, iphase) = position of that tracer in tracer_names
  subroutine wtr_index_table(cfg, tracer_names, idx, errmsg, errflg)
    type(wtr_config_type), intent(in)  :: cfg
    character(len=*),      intent(in)  :: tracer_names(:)
    integer,               intent(out) :: idx(:,:)
    character(len=*),      intent(out) :: errmsg
    integer,               intent(out) :: errflg
    character(len=32) :: names(wtr_nset(cfg)*wtr_nphase)
    integer           :: iset, ip, n

    errmsg = ''
    errflg = 0
    idx = 0
    call wtr_tracer_names(cfg, names)
    do iset = 1, wtr_nset(cfg)
       do ip = 1, wtr_nphase
          do n = 1, size(tracer_names)
             if (trim(tracer_names(n)) == trim(names((iset-1)*wtr_nphase + ip))) then
                idx(iset, ip) = n
                exit
             end if
          end do
          if (idx(iset, ip) == 0) then
             errmsg = 'water_tracers_config: water tracer not in tracer list: ' // trim(names((iset-1)*wtr_nphase + ip))
             errflg = 1
             return
          end if
       end do
    end do
  end subroutine wtr_index_table

end module water_tracers_config
