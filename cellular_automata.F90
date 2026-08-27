!L.Bengtsson, 2017-06
!P.Pegion, 2021-09
! swtich to new random number generator and improve computational efficiency
! and remove unsued code. Also add restart capability ca_global
!K.Barton, 2026-06 : modify to run on gaussian grid

!This program evolves a cellular automaton uniform over the globe

module cellular_automata_mod

use constants_mod, only : radius
use kinddef
use netcdf

implicit none
private

public :: ggca_grid_t
public :: ggca_init
public :: ggca_write
public :: cellular_automata

! Gaussian grid for cellular automata
type :: ggca_grid_t
  integer :: ntrunc = 0 ! Spectral truncation wavenumber
  integer :: nlat = 0   ! Number gaussian latitudes
  integer :: nlon = 0   ! Number longitudes
  real(kind=kind_io8), allocatable :: lats(:) ! latitude coords
  real(kind=kind_io8), allocatable :: lons(:) ! longitude coords
  real(kind=kind_io8) :: wlon = 0.0_kind_io8
  real(kind=kind_io8) :: rnlat = 0.0_kind_io8
end type ggca_grid_t

! Module variables for CA state persistence
integer, allocatable, save :: board_g(:,:,:)  ! CA board state (alive=1,dead=0)
integer, allocatable, save :: lives_g(:,:,:)  ! Cell lifetime counter
type(ggca_grid_t), save    :: ggrid           ! Gaussian grid information 
real(kind=kind_dbl_prec), allocatable, save :: ca_field(:,:,:) ! Normalized CA field for output

contains

subroutine cellular_automata(     &
        kstep,           & ! Current model time step 
        restart,         & ! Whether to initialize from restart (to do...)
        first_time_step, & ! Whether first time step  
        nca,             & ! Number of independent CA fields
        l_scale,         & ! Scale factor for determining GG resolution
        nlives,          & ! Maximum lifetime of CA cells
        nfracseed,       & ! Probabily threshold for seeding
        nseed,           & ! Frequency (in timesteps) for new seeding
        iseed_ca,        & ! random number seed
        ca_smooth,       & ! Whether to apply spatial smoothing
        nspinup,         & ! Number of spinup iterations on first step 
        nsmooth,         & ! Number of smoothing passes if ca_smooth=.true.
        ca_amplitude,    & ! Amplitude scaling factor for normalized output
        l_grid,          & ! Parent grid size (GG resolution based on l_grid / l_scale)
        ca_field_out,    & ! Normalized CA field for output
        grid_out)          ! Gaussian grid information (for remapping ca_field_out back to parent grid)

use kinddef,       only: kind_dbl_prec, kind_phys
use random_numbers,  only: random_01_CB

implicit none

integer, intent(in) :: kstep,nca,nlives,nseed,nspinup,nsmooth,l_scale
integer(kind=kind_dbl_prec), intent(in) :: iseed_ca
real(kind=kind_phys), intent(in) :: nfracseed,ca_amplitude,l_grid
logical, intent(in) :: ca_smooth,first_time_step, restart
real(kind=kind_dbl_prec), allocatable, intent(out) :: ca_field_out(:,:,:)
type(ggca_grid_t), intent(out) :: grid_out

integer :: i,j,k,nf, count4, ct, ntrunc
integer(8) :: csum
integer, allocatable :: iini_g(:,:,:), ilives_g(:,:,:)
real(kind=kind_dbl_prec) :: CAmean, CAstdv, psum, sq_diff
real(kind=kind_dbl_prec), allocatable :: noise(:,:,:)

if (nca .LT. 1) return

if (first_time_step) then
  print *, "Initializing Gaussian grid for CA..."
  call ggca_init(l_grid, l_scale, ggrid)

  ! Allocate module-level CA state fields
  if (allocated(ca_field)) deallocate(ca_field)
  allocate(ca_field(ggrid%nlon, ggrid%nlat, nca))
  ca_field(:,:,:) = 0.0

  if (allocated(board_g)) deallocate(board_g)
  allocate(board_g(ggrid%nlon,ggrid%nlat,nca))
  board_g(:,:,:) = 0

  if (allocated(lives_g)) deallocate(lives_g)
  allocate(lives_g(ggrid%nlon,ggrid%nlat,nca))
  lives_g(:,:,:) = 0
endif

! Allocated only for a single time step
allocate(iini_g  (ggrid%nlon, ggrid%nlat, nca)) !
allocate(ilives_g(ggrid%nlon, ggrid%nlat, nca)) ! Max lives in CA
allocate(noise   (ggrid%nlon, ggrid%nlat, nca)) ! Random noise
iini_g  (:,:,:) = 0
ilives_g(:,:,:) = 0
noise   (:,:,:) = 0.0

! Generate random noise for seeding/perturbation
do j=1,ggrid%nlat
  do i=1,ggrid%nlon
    call count_generator(iseed_ca, i+ggrid%nlon*(j-1), count4)
    do nf=1,nca
      noise(i,j,nf)=real(random_01_CB(nf*kstep,count4),kind=8)
    enddo
  enddo
enddo

! Initiate cellular automaton with random numbers larger than nfracseed
do nf=1,nca
  do j=1,ggrid%nlat
    do i=1,ggrid%nlon
      if (noise(i,j,nf) > nfracseed) then
        iini_g(i,j,nf)=1
      else
        iini_g(i,j,nf)=0
      endif
      ilives_g(i,j,nf)=int(real(nlives)*1.5*noise(i,j,nf))

      ! Initialize board and lives only on first time step
      if (first_time_step) then
        board_g(i,j,nf) = iini_g(i,j,nf)
        lives_g(i,j,nf) = ilives_g(i,j,nf)*iini_g(i,j,nf)
      endif
    enddo
  enddo
enddo

csum = int(ggrid%nlon, 8) * int(ggrid%nlat, 8) ! total # gaussian grid cells

do nf=1,nca ! Run update for each CA
  call update_cells(kstep, first_time_step, iseed_ca, nseed, nspinup, ggrid%nlon, ggrid%nlat, nca, nf, ilives_g)

!  ! Normalize output
!  psum    = SUM(ca_field(:,:,nf))
!  CAmean  = psum / real(csum, kind=kind_dbl_prec)
!  sq_diff = SUM((ca_field(:,:,nf) - CAmean)**2.0_kind_dbl_prec)
!  CAstdv  = sqrt(sq_diff / real(csum, kind=kind_dbl_prec))
!  ! Transform to mean of 1 and ca_amplitude standard deviation
!  ca_field(:,:,nf) = 1.0_kind_dbl_prec + (ca_field(:,:,nf) - CAmean) * (real(ca_amplitude, kind=kind_dbl_prec) / CAstdv)
!  ca_field(:,:,nf) = min(max(ca_field(:,:,nf), 0.0_kind_dbl_prec), 2.0_kind_dbl_prec)
enddo

! Populate output field and grid
ca_field_out = ca_field
grid_out = ggrid

end subroutine cellular_automata

subroutine update_cells( &  
        kstep,           & ! Current model time step 
        first_time_step, & ! Whether first time step  
        iseed_ca,        & ! random number seed
        nseed,           & ! Frequency (in timesteps) for new seeding
        nspinup,         & ! Number of spinup iterations on first step 
        nlon,            & ! Number longitudes on Gaussian grid
        nlat,            & ! Number latitudes on Gaussian grid
        nca,             & ! Total number of independent CAs
        nf,              & ! Current CA number
        ilives_g)          ! Max lifetime for cells

use random_numbers,  only: random_01_CB
use kinddef,         only: kind_dbl_prec

integer, intent(in)           :: kstep, nlon, nlat, nca, nf, nseed, nspinup
integer, intent(in)           :: ilives_g(nlon,nlat,nca)
integer(8), intent(in)        :: iseed_ca
logical, intent(in)           :: first_time_step
real(kind=kind_dbl_prec), dimension(nlon,nlat)    :: noise_b
integer, dimension(nlon,nlat) :: newseed, neighbors, birth, newcell
integer, allocatable          :: board_halo(:,:)
integer                       :: i,j,k,it,spinup,count4

allocate(board_halo(nlon+2,nlat+2)) ! For periodicity

! Seed new active cells
newseed=0
if (mod(kstep,nseed) == 0) then
  do j=1,nlat
    do i=1,nlon
      call count_generator(iseed_ca, i+nlon*(j-1), count4)
      noise_b(i,j)=real(random_01_CB(kstep,count4),kind=8)

      if(board_g(i,j,nf) == 0 .and. noise_b(i,j)>0.75) then
        newseed(i,j)=1
      endif

      board_g(i,j,nf) = board_g(i,j,nf) + newseed(i,j)
    enddo
  enddo
endif

if (first_time_step) then
  spinup = nspinup
else
  spinup = 1
endif

! Main CA iteration loop
do it=1,spinup
  neighbors=0
  birth=0
  newcell=0
  board_halo=0

  ! Setup periodicity
  board_halo(2:nlon+1,2:nlat+1) = board_g(:,:,nf)
  ! North pole : crossing pole shifts 180 degrees
  board_halo(1:nlon/2+1, 1) = board_halo(nlon/2+1:nlon+1, 2) 
  board_halo(nlon/2+2:nlon+2, 1) = board_halo(2:nlon/2+1, 2)
  ! South pole : crossing pole shifts 180 degrees
  board_halo(1:nlon/2+1, nlat+2) = board_halo(nlon/2+1:nlon+1, nlat+1) 
  board_halo(nlon/2+2:nlon+2, nlat+2) = board_halo(2:nlon/2+1, nlat+1)
  ! East-West periodicity
  board_halo(1, :) = board_halo(nlon+1, :)
  board_halo(nlon+2, :) = board_halo(2, :)
  ! Gaussian grid nlon will always be even

  ! Get neighbor count
  do j=2,nlat+1
    do i=2,nlon+1
      neighbors(i-1,j-1) = board_halo(i-1,j-1) + board_halo(i-1,j) + board_halo(i-1,j+1) + &
                           board_halo(i,  j-1)           +           board_halo(i,  j+1) + &
                           board_halo(i+1,j-1) + board_halo(i+1,j) + board_halo(i+1,j+1)
    enddo
  enddo

  ! Apply game of life logic
  do j=1,nlat
    do i=1,nlon

      if (neighbors(i,j)==2 .or. neighbors(i,j)==3) then
        birth(i,j) = 1
      else if (neighbors(i,j)<2 .or. neighbors(i,j)>3) then
        lives_g(i,j,nf)=lives_g(i,j,nf) - 1
      endif

      if (lives_g(i,j,nf)<0) then
        lives_g(i,j,nf) = 0
      endif

      if (birth(i,j)==1 .and. lives_g(i,j,nf)==0) then
        newcell(i,j) = 1
      endif

      lives_g(i,j,nf) = lives_g(i,j,nf) + newcell(i,j)*ilives_g(i,j,nf)

      if (neighbors(i,j)==3 .or. (board_g(i,j,nf)==1 .and. neighbors(i,j)==2)) then
        board_g(i,j,nf)=1
      else
        board_g(i,j,nf)=0
      endif
    enddo
  enddo
enddo ! spinup

ca_field(:,:,nf)=real(lives_g(:,:,nf), kind=kind_dbl_prec)

end subroutine update_cells

subroutine ggca_init(l_grid, l_scale, grid)

real(kind=kind_phys), intent(in) :: l_grid  ! Size of parent grid cell
integer, intent(in) :: l_scale ! Scale GG resolution to l_grid / l_scale
type(ggca_grid_t), intent(inout) :: grid

! Circumference of Earth
real(kind=kind_io8), parameter :: circ = 2.0_kind_io8 * acos(-1.0_kind_io8) * real(radius, kind_io8)
real(kind=kind_phys) :: l_min

l_min = l_grid / real(l_scale, kind=kind_phys)

! Generate truncation value from minimum resolved length scale 
print *, "Computing spectral truncation from l_min =", l_min, "meters"
grid%ntrunc = int(circ / l_min)
grid%ntrunc = ((grid%ntrunc + 1)/4)*4 + 2

grid%nlat = 2*grid%ntrunc + 2
grid%nlon = 2*grid%nlat

if (allocated(grid%lats)) deallocate(grid%lats)
allocate(grid%lats(grid%nlat))

if (allocated(grid%lons)) deallocate(grid%lons)
allocate(grid%lons(grid%nlon))

call build_gaussian_latitudes(grid%nlat, grid%lats)
call build_regular_longitudes(grid%nlon, grid%lons)

grid%wlon = grid%lons(1) - (grid%lons(2) - grid%lons(1))
grid%rnlat = grid%lats(1)*2.0_kind_io8 - grid%lats(2)

end subroutine ggca_init

subroutine ggca_write(outfile, grid, field, kstep, mype, root_pe)
! In/Out
character(len=*), intent(in) :: outfile
type(ggca_grid_t), intent(in) :: grid
real(kind=kind_io8), intent(in) :: field(grid%nlon, grid%nlat)
integer, intent(in) :: kstep
integer, intent(in), optional :: mype, root_pe

! Local
integer :: ncid, dlon, dlat, dtime, vlon, vlat, vtime, vfield, ntime, ierr
integer :: start_time(1), start_field(3), count_field(3)
real(kind=kind_io8), dimension(grid%nlon, grid%nlat, 1) :: field_out
logical :: file_exists

if (present(mype) .and. present(root_pe)) then
  if (mype /= root_pe) return
end if

if (.not. allocated(grid%lats)) then
  write(0,*) 'ggca_write: grid%lats not allocated'
  stop 1
end if
if (.not. allocated(grid%lons)) then
  write(0,*) 'ggca_write: grid%lons not allocated'
  stop 1
end if

inquire(file=trim(outfile), exist=file_exists)

if (.not. file_exists) then
  ierr = nf90_create(trim(outfile), NF90_CLOBBER, ncid)
  call check_nf90(ierr, 'nf90_create('//trim(outfile)//')')

  ierr = nf90_def_dim(ncid, "lon", grid%nlon, dlon)
  call check_nf90(ierr, 'nf90_def_dim(lon)')
  ierr = nf90_def_dim(ncid, "lat", grid%nlat, dlat)
  call check_nf90(ierr, 'nf90_def_dim(lat)')
  ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dtime)
  call check_nf90(ierr, 'nf90_def_dim(time)')

  ierr = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/dlon/), vlon)
  call check_nf90(ierr, 'nf90_def_var(lon)')
  ierr = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/dlat/), vlat)
  call check_nf90(ierr, 'nf90_def_var(lat)')
  ierr = nf90_def_var(ncid, "time", NF90_INT, (/dtime/), vtime)
  call check_nf90(ierr, 'nf90_def_var(time)')
  ierr = nf90_def_var(ncid, "field", NF90_DOUBLE, (/dlon, dlat, dtime/), vfield)
  call check_nf90(ierr, 'nf90_def_var(field)')

  ierr = nf90_enddef(ncid)
  call check_nf90(ierr, 'nf90_enddef')

  ierr = nf90_put_var(ncid, vlon, grid%lons)
  call check_nf90(ierr, 'nf90_put_var(lon)')
  ierr = nf90_put_var(ncid, vlat, grid%lats)
  call check_nf90(ierr, 'nf90_put_var(lat)')
else
  ierr = nf90_open(trim(outfile), NF90_WRITE, ncid)
  call check_nf90(ierr, 'nf90_open('//trim(outfile)//')')

  ierr = nf90_inq_dimid(ncid, "lon", dlon)
  call check_nf90(ierr, 'nf90_inq_dimid(lon)')
  ierr = nf90_inq_dimid(ncid, "lat", dlat)
  call check_nf90(ierr, 'nf90_inq_dimid(lat)')
  ierr = nf90_inq_dimid(ncid, "time", dtime)
  call check_nf90(ierr, 'nf90_inq_dimid(time)')

  ierr = nf90_inq_varid(ncid, "time", vtime)
  call check_nf90(ierr, 'nf90_inq_varid(time)')
  ierr = nf90_inq_varid(ncid, "field", vfield)
  call check_nf90(ierr, 'nf90_inq_varid(field)')
end if

ierr = nf90_inquire_dimension(ncid, dtime, len=ntime)
call check_nf90(ierr, 'nf90_inquire_dimension(time)')

start_time = (/ntime + 1/)
ierr = nf90_put_var(ncid, vtime, kstep, start=start_time)
call check_nf90(ierr, 'nf90_put_var(time)')

field_out(:,:,1) = field
start_field = (/1, 1, ntime + 1/)
count_field = (/grid%nlon, grid%nlat, 1/)
ierr = nf90_put_var(ncid, vfield, field_out, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_put_var(field)')

ierr = nf90_close(ncid)
call check_nf90(ierr, 'nf90_close')
end subroutine ggca_write

subroutine build_gaussian_latitudes(nlat, lats)
! Based on glats_stochy in spectral_transforms.F90
integer, intent(in) :: nlat
real(kind=kind_io8), intent(out) :: lats(nlat)

integer :: k, lgghaf
real(kind=kind_io8), allocatable :: colrad(:), wgt(:), rcs2(:)
real(kind=kind_io8), parameter   :: rad2deg = 180.0_kind_io8 / acos(-1.0_kind_io8)

integer :: k1, l2
real(kind=kind_qdt_prec) :: drad, dradz, p1, p2, pi, rad, rc
real(kind=kind_qdt_prec) :: rl2, scale, si, w, x

real(kind=kind_io8), parameter :: cons1 = 1.0_kind_io8
real(kind=kind_io8), parameter :: cons2 = 2.0_kind_io8
real(kind=kind_io8), parameter :: cons4 = 4.0_kind_io8
real(kind=kind_io8), parameter :: cons0p25 = 0.25_kind_io8

!#ifdef NO_QUAD_PRECISION
!real(kind=kind_qdt_prec), parameter :: eps = 1.0e-12_kind_qdt_prec
!#else
real(kind=kind_qdt_prec), parameter :: eps = 1.0e-20_kind_qdt_prec
!#endif

print *, "building gaussain lats"

lgghaf = nlat / 2
allocate(colrad(lgghaf), wgt(lgghaf), rcs2(lgghaf))

si = 1.0_kind_qdt_prec
l2 = 2*lgghaf
rl2 = real(l2, kind_qdt_prec)
scale = 2.0_kind_qdt_prec / (rl2 * rl2)
k1 = l2 - 1
pi = atan(si) * real(cons4, kind_qdt_prec)

dradz = pi / real(lgghaf, kind_qdt_prec) / 200.0_kind_qdt_prec
rad = 0.0_kind_qdt_prec

do k = 1, lgghaf
  drad = dradz
1 call poly(l2, rad, p2)
2 p1 = p2
  rad = rad + drad
  call poly(l2, rad, p2)
  if (sign(si, p1) == sign(si, p2)) go to 2
  if (drad < eps) go to 3
  rad = rad - drad
  drad = drad * real(cons0p25, kind_qdt_prec)
  go to 1
3 continue
  colrad(k) = real(rad, kind_io8)
  call poly(k1, rad, p1)
  x = cos(rad)
  w = scale * (1.0_kind_qdt_prec - x*x) / (p1*p1)
  wgt(k) = real(w, kind_io8)
  rc = 1.0_kind_qdt_prec / (sin(rad) * sin(rad))
  rcs2(k) = real(rc, kind_io8)
end do
! end glats_stochy

! Convert to degrees
do k = 1, lgghaf
  !lats(k) = -1.0_kind_io8 * colrad(lgghaf - k + 1) * rad2deg
  !lats(nlat - k + 1) = -1.0_kind_io8 * lats(k)
  lats(k) = (colrad(k) * rad2deg) - 90.0_kind_io8
  lats(nlat - k + 1) = -1.0_kind_io8 * lats(k)
end do

deallocate(colrad, wgt, rcs2)
end subroutine build_gaussian_latitudes

subroutine poly(n, rad, p)
! Based on poly routine from spectral_transforms
integer, intent(in) :: n
real(kind=kind_qdt_prec), intent(in) :: rad
real(kind=kind_qdt_prec), intent(out) :: p

integer :: i
real(kind=kind_qdt_prec) :: g, x, y1, y2, y3

x = cos(rad)
y1 = 1.0_kind_qdt_prec
y2 = x
y3 = y2

do i = 2, n
  g = x*y2
  y3 = g - y1 + g - (g-y1) / real(i, kind_qdt_prec)
  y1 = y2
  y2 = y3
end do

p = y3
end subroutine poly

! based on init_stochastic_physics in stochastic_physics.F90
subroutine build_regular_longitudes(nlon, lons)
integer, intent(in) :: nlon
real(kind=kind_io8), intent(out) :: lons(nlon)

integer :: i
real(kind=kind_io8) :: dx

print *, "building regular lons"

dx = 360.0_kind_io8 / real(nlon, kind_io8)
do i = 1, nlon
  lons(i) = dx * real(i - 1, kind_io8)
end do
end subroutine build_regular_longitudes

subroutine count_generator(iseed_ca, loc_val, count4)
integer(kind=kind_dbl_prec), intent(in) :: iseed_ca
integer, intent(in)  :: loc_val
integer, intent(out) ::count4 
integer(8) :: count, count_rate, count_max, count_trunc, iscale=10000000000_8

if (iseed_ca <= 0) then ! Generate random seed from system cloc
  call system_clock(count, count_rate, count_max)
  ! iseed is elapsed time since unix epoch began (secs)
  ! truncate to 4 byte integer
  count_trunc = iscale*(count/iscale)
  count4 = count - count_trunc*loc_val
else ! Use seed for reproducibility
  count4 = int(mod(int(iseed_ca*loc_val, 8) + 2147483648_8, 4294967296_8) - 2147483648_8)
endif
end subroutine count_generator

subroutine check_nf90(status, where)
use netcdf
integer, intent(in) :: status
character(len=*), intent(in) :: where
if (status /= nf90_noerr) then
  write(0,*) 'NetCDF error in ', trim(where), ': ', trim(nf90_strerror(status))
  stop 2
end if
end subroutine check_nf90

end module cellular_automata_mod
