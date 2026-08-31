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

public :: ca_grid_t                    ! Cellular Automata Object Type
public :: initialize_cellular_automata ! Initialize brand new CA
public :: restart_cellular_automata    ! Initialize CA from restart files 
public :: update_cellular_automata     ! Run CA update logic
public :: output_cellular_automata     ! Output CA data to netcdf

! Cellular Automata object
type :: ca_grid_t 
  integer :: nlat = 0   ! Number gaussian latitudes
  integer :: nlon = 0   ! Number longitudes
  integer :: nca = 0    ! Number of independent CAs to run
  integer, allocatable :: board(:,:,:)  ! CA board state (alive=1,dead=0)
  integer, allocatable :: lives(:,:,:)  ! Cell lifetime counter
  real(kind=kind_io8), allocatable :: lats(:)   ! latitude coords
  real(kind=kind_io8), allocatable :: lons(:)   ! longitude coords
  real(kind=kind_dbl_prec), allocatable :: field(:,:,:) ! Normalized CA field for output
  real(kind=kind_io8) :: wlon = 0.0_kind_io8
  real(kind=kind_io8) :: rnlat = 0.0_kind_io8
  logical :: first_time_step = .True.
end type ca_grid_t

contains

! -------------------------------------- !
! Initialize brand new cellular automata !
! -------------------------------------- !
subroutine initialize_cellular_automata(     &
        l_grid,          & ! Parent grid size (GG resolution based on l_grid / l_scale)
        l_scale,         & ! Scale factor for determining GG resolution
        nca,             & ! Number of independent CA fields
        grid)              ! CA grid for output

use kinddef,       only: kind_dbl_prec, kind_phys
use random_numbers,  only: random_01_CB

implicit none
! Input/Output Variables
real(kind=kind_phys), intent(in) :: l_grid
integer, intent(in) :: l_scale
integer, intent(in) :: nca
type(ca_grid_t), intent(out) :: grid

! Private Variables
real(kind=kind_io8), parameter :: circ = 2.0_kind_io8 * acos(-1.0_kind_io8) * real(radius, kind_io8)
real(kind=kind_phys) :: l_min, wlon, rnlat
real(kind=kind_io8), allocatable :: lats(:)   ! latitude coords
real(kind=kind_io8), allocatable :: lons(:)   ! longitude coords
integer :: ntrunc, nlon, nlat

if (nca .LT. 1) return

! Build Gaussian Grid
l_min = l_grid / real(l_scale, kind=kind_phys)

print *, "Computing spectral truncation from l_min =", l_min, "meters"
ntrunc = int(circ / l_min)
ntrunc = ((ntrunc + 1)/4)*4 + 2

nlat = 2*ntrunc + 2
nlon = 2*nlat

allocate(lats(nlat))
lats(:) = 0.0
allocate(lons(nlon))
lons(:) = 0.0

call build_gaussian_latitudes(nlat, lats)
call build_regular_longitudes(nlon, lons)

wlon  = lons(1) - (lons(2) - lons(1))
rnlat = lats(1)*2.0_kind_io8 - lats(2)

! Populate CA grid object
grid%nca = nca
grid%nlat = nlat
grid%nlon = nlon
grid%wlon = wlon
grid%rnlat = rnlat
grid%first_time_step = .True.

if (allocated(grid%lats)) deallocate(grid%lats)
allocate(grid%lats(nlat))
grid%lats(:) = lats(:)

if (allocated(grid%lons)) deallocate(grid%lons)
allocate(grid%lons(nlon))
grid%lons(:) = lons(:)

if (allocated(grid%board)) deallocate(grid%board)
allocate(grid%board(nlon, nlat, nca))
grid%board(:,:,:) = 0

if (allocated(grid%lives)) deallocate(grid%lives)
allocate(grid%lives(nlon, nlat, nca))
grid%lives(:,:,:) = 0

if (allocated(grid%field)) deallocate(grid%field)
allocate(grid%field(nlon, nlat, nca))
grid%field(:,:,:) = 0.0

end subroutine initialize_cellular_automata

! ----------------------------------------------- !
! Initialize cellular automata from restart files !
! ----------------------------------------------- !
subroutine restart_cellular_automata(     &
        restart_file,    &
        time_index,      &
        grid)

implicit none
! Input/Output Variables
character(len=*), intent(in) :: restart_file 
integer, intent(in) :: time_index
type(ca_grid_t), intent(out) :: grid

! Local variables
integer :: ncid, dlon, dlat, dnca, dtime, ierr
integer :: vwlon, vrnlat, vlon, vlat, vtime, vboard, vlives, vfield
integer :: start_field(4), count_field(4)
logical :: file_exists

! Check file exists, initialize from scratch if it doesn't
inquire(file=trim(restart_file), exist=file_exists)
if (.not. file_exists) then
  write(0,*) 'restart_cellular_automata: file not found - ', trim(restart_file)
  stop 1
endif

ierr = nf90_open(trim(restart_file), NF90_NOWRITE, ncid)
call check_nf90(ierr, 'nf90_open('//trim(restart_file)//')')

! Gather dimensions
ierr = nf90_inq_dimid(ncid, "lon", dlon)
call check_nf90(ierr, 'nf90_inq_dimid(lon)')
ierr = nf90_inquire_dimension(ncid, dlon, len=grid%nlon)

ierr = nf90_inq_dimid(ncid, "lat", dlat)
call check_nf90(ierr, 'nf90_inq_dimid(lat)')
ierr = nf90_inquire_dimension(ncid, dlat, len=grid%nlat)

ierr = nf90_inq_dimid(ncid, "nca", dnca)
call check_nf90(ierr, 'nf90_inq_dimid(nca)')
ierr = nf90_inquire_dimension(ncid, dnca, len=grid%nca)

! Gather variables
ierr = nf90_inq_varid(ncid, "wlon", vwlon)
call check_nf90(ierr, 'nf90_inq_varid(wlon)')
ierr = nf90_inq_varid(ncid, "rnlat", vrnlat)
call check_nf90(ierr, 'nf90_inq_varid(rnlat)')
ierr = nf90_inq_varid(ncid, "lon", vlon)
call check_nf90(ierr, 'nf90_inq_varid(lon)')
ierr = nf90_inq_varid(ncid, "lat", vlat)
call check_nf90(ierr, 'nf90_inq_varid(lat)')
ierr = nf90_inq_varid(ncid, "board", vboard)
call check_nf90(ierr, 'nf90_inq_varid(board)')
ierr = nf90_inq_varid(ncid, "lives", vlives)
call check_nf90(ierr, 'nf90_inq_varid(lives)')
ierr = nf90_inq_varid(ncid, "field", vfield)
call check_nf90(ierr, 'nf90_inq_varid(field)')

! Populate CA object
if (allocated(grid%lons))  deallocate(grid%lons)
allocate(grid%lons(grid%nlon))
if (allocated(grid%lats))  deallocate(grid%lats)
allocate(grid%lats(grid%nlat))
if (allocated(grid%board)) deallocate(grid%board)
allocate(grid%board(grid%nlon, grid%nlat, grid%nca))
if (allocated(grid%lives)) deallocate(grid%lives)
allocate(grid%lives(grid%nlon, grid%nlat, grid%nca))
if (allocated(grid%field)) deallocate(grid%field)
allocate(grid%field(grid%nlon, grid%nlat, grid%nca))

ierr = nf90_get_var(ncid, vwlon, grid%wlon)
call check_nf90(ierr, 'nf90_get_var(wlon)')
ierr = nf90_get_var(ncid, vrnlat, grid%rnlat)
call check_nf90(ierr, 'nf90_get_var(rnlat)')
ierr = nf90_get_var(ncid, vlon, grid%lons)
call check_nf90(ierr, 'nf90_get_var(lon)')
ierr = nf90_get_var(ncid, vlat, grid%lats)
call check_nf90(ierr, 'nf90_get_var(lat)')

start_field = (/ 1, 1, 1, time_index /)
count_field = (/ grid%nlon, grid%nlat, grid%nca, 1 /)

ierr = nf90_get_var(ncid, vboard, grid%board, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_get_var(board)')
ierr = nf90_get_var(ncid, vlives, grid%lives, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_get_var(lives)')
ierr = nf90_get_var(ncid, vfield, grid%field, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_get_var(field)')

ierr = nf90_close(ncid)
call check_nf90(ierr, 'nf90_close()')

! For restart it is not first_time_step
grid%first_time_step = .false.

end subroutine restart_cellular_automata

! ------------------------------------------ !
! Advance cellular automata by one time step !
! ------------------------------------------ !
subroutine update_cellular_automata( &  
        kstep,           & ! Current model time step 
        iseed_ca,        & ! random number seed
        nfracseed,       & ! 
        nseed,           & ! Frequency (in timesteps) for new seeding
        nspinup,         & ! Number of spinup iterations on first step
        nlives,          &
        ca_grid)

use random_numbers,  only: random_01_CB
use kinddef,         only: kind_dbl_prec

! In/out variables
integer, intent(in)             :: kstep, nseed, nspinup, nlives
integer(8), intent(in)          :: iseed_ca
real(kind=kind_phys), intent(in):: nfracseed
type(ca_grid_t), intent(inout)  :: ca_grid

! Private variables
real(kind=kind_dbl_prec), allocatable :: noise(:,:,:), noise_b(:,:)
integer, allocatable          :: board_halo(:,:), newseed(:,:), neighbors(:,:), birth(:,:), newcell(:,:)
integer, allocatable          :: ilives(:,:,:), iini(:,:,:)
integer                       :: i,j,k,it,spinup,count4,nf

! Allocated only for a single time step
allocate(iini     (ca_grid%nlon, ca_grid%nlat, ca_grid%nca)) !
allocate(ilives   (ca_grid%nlon, ca_grid%nlat, ca_grid%nca)) ! Max lives in CA
allocate(noise    (ca_grid%nlon, ca_grid%nlat, ca_grid%nca)) ! Random noise
allocate(noise_b  (ca_grid%nlon, ca_grid%nlat)) ! Random noise
allocate(newseed  (ca_grid%nlon, ca_grid%nlat)) 
allocate(neighbors(ca_grid%nlon, ca_grid%nlat)) 
allocate(birth    (ca_grid%nlon, ca_grid%nlat)) 
allocate(newcell  (ca_grid%nlon, ca_grid%nlat)) 
allocate(board_halo(ca_grid%nlon+2,ca_grid%nlat+2)) ! For periodicity
iini     (:,:,:) = 0
ilives   (:,:,:) = 0
noise    (:,:,:) = 0.0
noise_b  (:,:) = 0.0

! Generate random noise for seeding/perturbation
do j=1,ca_grid%nlat
  do i=1,ca_grid%nlon
    call count_generator(iseed_ca, i+ca_grid%nlon*(j-1), count4)
    do nf=1,ca_grid%nca
      noise(i,j,nf)=real(random_01_CB(nf*kstep,count4),kind=8)
    enddo
  enddo
enddo

! Initiate cellular automata with random numbers larger than nfracseed
do nf=1,ca_grid%nca
  newseed  (:,:) = 0
  do j=1,ca_grid%nlat
    do i=1,ca_grid%nlon

      if (noise(i,j,nf) > nfracseed) then
        iini  (i,j,nf)=1
      else
        iini  (i,j,nf)=0
      endif

      ilives  (i,j,nf)=int(real(nlives)*1.5*noise(i,j,nf))

      if (ca_grid%first_time_step) then
        ! Initialize board and lives only on first time step
        ca_grid%board(i,j,nf) = iini(i,j,nf)
        ca_grid%lives(i,j,nf) = ilives(i,j,nf)*iini(i,j,nf)
      endif

      ! Seed new active cells
      if (mod(kstep, nseed) == 0) then
        call count_generator(iseed_ca, i+ca_grid%nlon*(j-1), count4)
  
        noise_b(i,j)=real(random_01_CB(kstep,count4),kind=8)
  
        if(ca_grid%board(i,j,nf) == 0 .and. noise_b(i,j)>0.75) then
          newseed(i,j)=1
        endif
  
        ca_grid%board(i,j,nf) = ca_grid%board(i,j,nf) + newseed(i,j)
      endif

    enddo
  enddo
enddo

if (ca_grid%first_time_step) then
  spinup = nspinup
else
  spinup = 1
endif

! Main CA iteration loop
do nf=1,ca_grid%nca
  do it=1,spinup
    neighbors(:,:) = 0
    birth    (:,:) = 0
    newcell  (:,:) = 0
    board_halo(:,:) = 0
  
    ! Setup periodicity
    board_halo(2:ca_grid%nlon+1,2:ca_grid%nlat+1) = ca_grid%board(:,:,nf)
    ! North pole : crossing pole shifts 180 degrees
    board_halo(1:ca_grid%nlon/2+1, 1) = board_halo(ca_grid%nlon/2+1:ca_grid%nlon+1, 2) 
    board_halo(ca_grid%nlon/2+2:ca_grid%nlon+2, 1) = board_halo(2:ca_grid%nlon/2+1, 2)
    ! South pole : crossing pole shifts 180 degrees
    board_halo(1:ca_grid%nlon/2+1, ca_grid%nlat+2) = board_halo(ca_grid%nlon/2+1:ca_grid%nlon+1, ca_grid%nlat+1) 
    board_halo(ca_grid%nlon/2+2:ca_grid%nlon+2, ca_grid%nlat+2) = board_halo(2:ca_grid%nlon/2+1, ca_grid%nlat+1)
    ! East-West periodicity
    board_halo(1, :) = board_halo(ca_grid%nlon+1, :)
    board_halo(ca_grid%nlon+2, :) = board_halo(2, :)
    ! Gaussian grid nlon will always be even
  
    ! Get neighbor count
    do j=2,ca_grid%nlat+1
      do i=2,ca_grid%nlon+1
        neighbors(i-1,j-1) = board_halo(i-1,j-1) + board_halo(i-1,j) + board_halo(i-1,j+1) + &
                             board_halo(i,  j-1)           +           board_halo(i,  j+1) + &
                             board_halo(i+1,j-1) + board_halo(i+1,j) + board_halo(i+1,j+1)
      enddo
    enddo
  
    ! Apply game of life logic
    do j=1,ca_grid%nlat
      do i=1,ca_grid%nlon
  
        if (neighbors(i,j)==2 .or. neighbors(i,j)==3) then
          birth(i,j) = 1
        else if (neighbors(i,j)<2 .or. neighbors(i,j)>3) then
          ca_grid%lives(i,j,nf)=ca_grid%lives(i,j,nf) - 1
        endif
  
        if (ca_grid%lives(i,j,nf)<0) then
          ca_grid%lives(i,j,nf) = 0
        endif
  
        if (birth(i,j)==1 .and. ca_grid%lives(i,j,nf)==0) then
          newcell(i,j) = 1
        endif
  
        ca_grid%lives(i,j,nf) = ca_grid%lives(i,j,nf) + newcell(i,j)*ilives(i,j,nf)
  
        if (neighbors(i,j)==3 .or. (ca_grid%board(i,j,nf)==1 .and. neighbors(i,j)==2)) then
          ca_grid%board(i,j,nf)=1
        else
          ca_grid%board(i,j,nf)=0
        endif
      enddo
    enddo
  enddo 
enddo! spinup

ca_grid%field(:,:,:) = real(ca_grid%lives(:,:,:), kind=kind_dbl_prec)

!! Normalize output
!csum = int(ca_grid%nlon, 8) * int(ca_grid%nlat, 8) ! total # gaussian grid cells
!do nf=1,ca_grid%nca
!  psum    = SUM(ca_grid%field(:,:,nf))
!  CAmean  = psum / real(csum, kind=kind_dbl_prec)
!  sq_diff = SUM((ca_grid%field(:,:,nf) - CAmean)**2.0_kind_dbl_prec)
!  CAstdv  = sqrt(sq_diff / real(csum, kind=kind_dbl_prec))
!  ! Transform to mean of 1 and ca_amplitude standard deviation
!  field(:,:,nf) = 1.0_kind_dbl_prec + (ca_grid%field(:,:,nf) - CAmean) * (real(ca_amplitude, kind=kind_dbl_prec) / CAstdv)
!  field(:,:,nf) = min(max(ca_grid%field(:,:,nf), 0.0_kind_dbl_prec), 2.0_kind_dbl_prec)
!enddo

deallocate(iini)
deallocate(ilives)
deallocate(noise)
deallocate(noise_b)
deallocate(newseed)  
deallocate(neighbors)
deallocate(birth) 
deallocate(newcell)
deallocate(board_halo)

if (ca_grid%first_time_step) ca_grid%first_time_step = .False.

end subroutine update_cellular_automata 

! ------------------------------------------- !
! Build gg lats based on glats_stochy routine !
! ------------------------------------------- !
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

! Helper function for building gg lats
! Based on poly routine from spectral_transforms
subroutine poly(n, rad, p)
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

! ------------------------------------------ !
! Build gg lons based on spectral_transforms !
! ------------------------------------------ !
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

! ----------------------------- !
! Generates random count4 value !
! ----------------------------- !
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

! ----------------------------------- !
! Output cellular automata to netcdf  !
! ----------------------------------- !
subroutine output_cellular_automata(outfile, grid, time)
! In/Out
character(len=*), intent(in) :: outfile
type(ca_grid_t), intent(in) :: grid
integer, intent(in) :: time 

! Local
integer :: ncid, dlon, dlat, dnca, dtime, vwlon, vrnlat, vlon, vlat, vtime, vboard, vlives, vfield, ntime, ierr
integer :: start_time(1), start_field(4), count_field(3)
real(kind=kind_io8), dimension(grid%nlon, grid%nlat, 1) :: field_out
logical :: file_exists

! CA type vars
  integer :: nlat = 0   ! Number gaussian latitudes
  integer :: nlon = 0   ! Number longitudes
  integer :: nca = 0    ! Number of independent CAs to run
  integer, allocatable :: board(:,:,:)  ! CA board state (alive=1,dead=0)
  integer, allocatable :: lives(:,:,:)  ! Cell lifetime counter
  real(kind=kind_io8), allocatable :: lats(:)   ! latitude coords
  real(kind=kind_io8), allocatable :: lons(:)   ! longitude coords
  real(kind=kind_dbl_prec), allocatable :: field(:,:,:) ! Normalized CA field for output
  real(kind=kind_io8) :: wlon = 0.0_kind_io8
  real(kind=kind_io8) :: rnlat = 0.0_kind_io8
  logical :: first_time_step = .True.

if (.not. allocated(grid%lats)) then
  write(0,*) 'output_cellular_automata: grid%lats not allocated'
  stop 1
end if
if (.not. allocated(grid%lons)) then
  write(0,*) 'output_cellular_automata: grid%lons not allocated'
  stop 1
end if

inquire(file=trim(outfile), exist=file_exists)

if (.not. file_exists) then
  ierr = nf90_create(trim(outfile), NF90_CLOBBER, ncid)
  call check_nf90(ierr, 'nf90_create('//trim(outfile)//')')

  ! Define dimensions
  ierr = nf90_def_dim(ncid, "lon", grid%nlon, dlon)
  call check_nf90(ierr, 'nf90_def_dim(lon)')
  ierr = nf90_def_dim(ncid, "lat", grid%nlat, dlat)
  call check_nf90(ierr, 'nf90_def_dim(lat)')
  ierr = nf90_def_dim(ncid, "nca", grid%nca, dnca)
  call check_nf90(ierr, 'nf90_def_dim(nca)')
  ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dtime)
  call check_nf90(ierr, 'nf90_def_dim(time)')

  ! Define variables
  ierr = nf90_def_var(ncid, "wlon", NF90_DOUBLE, vwlon)
  call check_nf90(ierr, 'nf90_dev_var(wlon)')
  ierr = nf90_def_var(ncid, "rnlat", NF90_DOUBLE, vrnlat)
  call check_nf90(ierr, 'nf90_dev_var(rnlat)')
  ierr = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/dlon/), vlon)
  call check_nf90(ierr, 'nf90_def_var(lon)')
  ierr = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/dlat/), vlat)
  call check_nf90(ierr, 'nf90_def_var(lat)')
  ierr = nf90_def_var(ncid, "time", NF90_INT, (/dtime/), vtime)
  call check_nf90(ierr, 'nf90_def_var(time)')
  ierr = nf90_def_var(ncid, "board", NF90_INT, (/dlon, dlat, dnca, dtime/), vboard)
  call check_nf90(ierr, 'nf90_def_var(board)')
  ierr = nf90_def_var(ncid, "lives", NF90_INT, (/dlon, dlat, dnca, dtime/), vlives)
  call check_nf90(ierr, 'nf90_def_var(lives)')
  ierr = nf90_def_var(ncid, "field", NF90_DOUBLE, (/dlon, dlat, dnca, dtime/), vfield)
  call check_nf90(ierr, 'nf90_def_var(field)')

  ierr = nf90_enddef(ncid)
  call check_nf90(ierr, 'nf90_enddef')

  ! Place vars that aren't time-dependent
  ierr = nf90_put_var(ncid, vlon, grid%lons)
  call check_nf90(ierr, 'nf90_put_var(lon)')
  ierr = nf90_put_var(ncid, vlat, grid%lats)
  call check_nf90(ierr, 'nf90_put_var(lat)')
  ierr = nf90_put_var(ncid, vwlon, grid%wlon)
  call check_nf90(ierr, 'nf90_put_var(wlon)')
  ierr = nf90_put_var(ncid, vrnlat, grid%rnlat)
  call check_nf90(ierr, 'nf90_put_var(rnlat)')

else
  ierr = nf90_open(trim(outfile), NF90_WRITE, ncid)
  call check_nf90(ierr, 'nf90_open('//trim(outfile)//')')

  ierr = nf90_inq_dimid(ncid, "lon", dlon)
  call check_nf90(ierr, 'nf90_inq_dimid(lon)')
  ierr = nf90_inq_dimid(ncid, "lat", dlat)
  call check_nf90(ierr, 'nf90_inq_dimid(lat)')
  ierr = nf90_inq_dimid(ncid, "nca", dnca)
  call check_nf90(ierr, 'nf90_inq_dimid(nca)')
  ierr = nf90_inq_dimid(ncid, "time", dtime)
  call check_nf90(ierr, 'nf90_inq_dimid(time)')

  ierr = nf90_inq_varid(ncid, "time", vtime)
  call check_nf90(ierr, 'nf90_inq_varid(time)')
  ierr = nf90_inq_varid(ncid, "board", vboard)
  call check_nf90(ierr, 'nf90_inq_varid(board)')
  ierr = nf90_inq_varid(ncid, "lives", vlives)
  call check_nf90(ierr, 'nf90_inq_varid(lives)')
  ierr = nf90_inq_varid(ncid, "field", vfield)
  call check_nf90(ierr, 'nf90_inq_varid(field)')
end if

ierr = nf90_inquire_dimension(ncid, dtime, len=ntime)
call check_nf90(ierr, 'nf90_inquire_dimension(time)')

start_time = (/ntime + 1/)
ierr = nf90_put_var(ncid, vtime, time, start=start_time)
call check_nf90(ierr, 'nf90_put_var(time)')

start_field = (/1, 1, 1, ntime + 1/)
count_field = (/grid%nlon, grid%nlat, grid%nca/)
ierr = nf90_put_var(ncid, vfield, grid%field, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_put_var(field)')
ierr = nf90_put_var(ncid, vboard, grid%board, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_put_var(board)')
ierr = nf90_put_var(ncid, vlives, grid%lives, start=start_field, count=count_field)
call check_nf90(ierr, 'nf90_put_var(lives)')

ierr = nf90_close(ncid)
call check_nf90(ierr, 'nf90_close')
end subroutine output_cellular_automata

! NetCDF output helper function
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
