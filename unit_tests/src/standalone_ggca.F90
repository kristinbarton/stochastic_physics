program standalone_ggca

use netcdf
use mpi_f08
use cellular_automata_mod
use atmosphere_stub_mod, only: Atm, atmosphere_init_stub
use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe
use fms_mod, only: fms_init
use xgrid_mod, only: grid_box_type
use spectral_transforms, only: stochy_la2ga
use kinddef, only: kind_phys, kind_dbl_prec

implicit none

integer :: ntasks, nblks, blksz, ierr, my_id, root_pe, i, j, nx, ny
integer :: isc, iec, jsc, jec, nb
integer :: istart, iend, dump_time, ct, nsteps
real(kind=4) :: ts

type(grid_box_type) :: grid_box
type(ca_grid_t) :: ca_object

real(kind=kind_dbl_prec), allocatable :: fv3_ca(:,:,:)
real(kind=kind_dbl_prec), allocatable :: fv3_lons(:,:), fv3_lats(:,:)

character(len=128) :: output_file, filename
character(len=4) :: strid

! Cellular automata controls
integer :: nca, nlives, ncells
real(kind=kind_phys) :: nfracseed, rcell
integer :: nseed, nseed_g
logical :: do_ca, ca_sgs, ca_global, ca_smooth
integer*8 :: iseed_ca
integer :: nspinup, nsmooth
real :: ca_amplitude
real(kind=kind_phys) :: l_grid
logical :: ca_closure, ca_entr, ca_trigger, warm_start, write_restart
integer :: time_idx

namelist /gfs_physics_nml/ do_ca, ca_sgs, ca_global, nca, ncells, nlives, nseed,       &
                          nfracseed, rcell, ca_trigger, ca_entr, ca_closure, nca,    &
                          ncells, nlives, nseed, ca_smooth, nspinup, iseed_ca,   &
                          nsmooth, ca_amplitude, warm_start

! Defaults
warm_start     = .false.
write_restart  = .true.
nca            = 3
ncells         = 1
nlives         = 1
nfracseed      = 0.5
nseed          = 100000
iseed_ca       = 0
nspinup        = 1
do_ca          = .true.
ca_sgs         = .false.
ca_smooth      = .false.
ca_amplitude   = 500.
rcell          = 0.0
nsteps         = 100

! open namelist file
open (unit=565, file='input.nml', status='OLD', iostat=ierr)
read(565,gfs_physics_nml)
close(565)

print *, 'ca_sgs,ca_global',ca_sgs,ca_global
if (.not. do_ca) nca = 0

! initialize fms
call mpp_init()
call fms_init
my_id=mpp_pe()
root_pe = mpp_root_pe()
ntasks = mpp_npes()

call atmosphere_init_stub(grid_box)
!define domain
isc = Atm(1)%bd%isc
iec = Atm(1)%bd%iec
jsc = Atm(1)%bd%jsc
jec = Atm(1)%bd%jec

print *, 'ATM npx,npy=', Atm(1)%npx, Atm(1)%npy
nx = iec - isc + 1
ny = jec - jsc + 1
print *, 'after init', my_id, Atm(1)%tile_of_mosaic, isc, jec

! for this simple test, nblocks = ny, blksz=ny
blksz=nx
nblks=ny

! set gaussian grid scale l_grid based on fv3 grid size
l_grid = 0.5_kind_phys * (sum(grid_box%dx(isc:iec,jsc:jec)) + sum(grid_box%dy(isc:iec,jsc:jec))) / real(nx*ny, kind=kind_phys) * 5
print *, 'average FV3 cell size (m) = ', l_grid

if (do_ca) then
  allocate(fv3_ca(nblks,blksz,max(1,nca)))
  allocate(fv3_lons(nblks,blksz))
  allocate(fv3_lats(nblks,blksz))

  do j=jsc,jec
    do i=isc,iec
      nb = j-jsc+1
      fv3_lons(nb,i-isc+1) = Atm(1)%gridstruct%agrid_64(i,j,1) * 180.0_kind_dbl_prec / acos(-1.0_kind_dbl_prec)
      fv3_lats(nb,i-isc+1) = Atm(1)%gridstruct%agrid_64(i,j,2) * 180.0_kind_dbl_prec / acos(-1.0_kind_dbl_prec)
    enddo
  enddo

  ! Setup
  if (warm_start) then
    time_idx=nsteps/5
    call restart_cellular_automata("RESTART/ca_restart.nc", time_idx, ca_object)
    istart = nsteps 
    iend = istart + nsteps 
  else
    call initialize_cellular_automata(l_grid, ncells, nca, ca_object)
    istart = 0
    iend = istart + nsteps
  endif
  ct = 1
  strid = get_strid(my_id, ntasks)
  
  do i = istart, iend 
    print *, "kstep = ", i
  
    call update_cellular_automata(i, iseed_ca, nfracseed, nseed, nspinup, nlives, ca_object)
    call remap_ca_to_fv3(ca_object, fv3_lons, fv3_lats, fv3_ca)
  
    if (mod(i,5) == 0) then
      write(filename, '("RESTART/ca_out.tile", A, ".nc")') trim(strid)
      call fv3ca_write(trim(filename), fv3_ca, i, ierr)
      if (write_restart .and. my_id == 0) then
        print *, "writing restart: my_id",my_id,"time",i
        write(filename, '("RESTART/ca_restart.nc")')
        call output_cellular_automata(trim(filename), ca_object, i)
      endif
    endif
  enddo
endif
contains

subroutine remap_ca_to_fv3(gaussian_grid, target_lons, target_lats, fv3_fields)
  type(ca_grid_t), intent(in) :: gaussian_grid
  real(kind=kind_dbl_prec), intent(in) :: target_lons(:,:), target_lats(:,:)
  real(kind=kind_dbl_prec), intent(out) :: fv3_fields(:,:,:)

  integer :: block, nf, nblocks, blocksize

  nblocks = size(target_lons, 1)
  blocksize = size(target_lons, 2)
  do nf=1,size(gaussian_grid%field, 3)
    do block=1,nblocks
      call stochy_la2ga(                                  &
             gaussian_grid%field(:,:,nf),                 &
             gaussian_grid%nlon, gaussian_grid%nlat,      &
             gaussian_grid%lons, gaussian_grid%lats,      &
             gaussian_grid%wlon, gaussian_grid%rnlat,     &
             fv3_fields(block,1:blocksize, nf),           &
             blocksize,                                   &
             target_lats(block,:), target_lons(block,:)   &
           )
    enddo
  enddo
end subroutine remap_ca_to_fv3

subroutine fv3ca_write(filename, ca, time, ierr)
  use netcdf
  use kinddef, only : kind_phys
  ! In/Out variables
  character(len=*), intent(in)         :: filename
  real(kind=kind_dbl_prec), intent(in) :: ca(:,:,:)
  integer, intent(in)             :: time
  integer, intent(out)                 :: ierr

  ! Local Variables
  integer :: ncid, xdim_id, ydim_id, time_dim_id
  integer :: xt_var_id, yt_var_id, time_var_id
  integer :: nx_out, ny_out, nc_out, k
  integer, allocatable :: ca_ids(:)
  real(kind=4), allocatable :: grid_xt(:), grid_yt(:), workg(:,:)
  real(kind=4) :: undef
  logical :: file_exists
  character(len=10) :: var_name

  nc_out = size(ca, 3)
  nx_out = size(ca, 2)
  ny_out = size(ca, 1)

  allocate(ca_ids(nc_out))
  allocate(grid_xt(nx_out), grid_yt(ny_out), workg(nx_out, ny_out))

  undef = 9.99e+20
  grid_xt = [(real(i, kind=kind_phys), i=1,nx_out)]
  grid_yt = [(real(j, kind=kind_phys), j=1,ny_out)]

  inquire(file=trim(filename), exist=file_exists)

if (.not. file_exists) then
  ierr = nf90_create(filename, NF90_CLOBBER, ncid)

  ! Define Dimensions
  ierr = nf90_def_dim(ncid, 'grid_xt', nx_out, xdim_id)
  ierr = nf90_def_dim(ncid, 'grid_yt', ny_out, ydim_id)
  ierr = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim_id)

  ! Define Coordinates
  ierr = nf90_def_var(ncid, 'grid_xt', NF90_FLOAT, (/xdim_id/), xt_var_id)
  ierr = nf90_put_att(ncid, xt_var_id, 'long_name', 'T-cell longitude')
  ierr = nf90_put_att(ncid, xt_var_id, 'cartesian_axis', 'X')
  ierr = nf90_put_att(ncid, xt_var_id, 'units', 'degrees_E')

  ierr = nf90_def_var(ncid, 'grid_yt', NF90_FLOAT, (/ydim_id/), yt_var_id)
  ierr = nf90_put_att(ncid, yt_var_id, 'long_name', 'T-cell latitude')
  ierr = nf90_put_att(ncid, yt_var_id, 'cartesian_axis', 'Y')
  ierr = nf90_put_att(ncid, yt_var_id, 'units', 'degrees_N')

  ierr = nf90_def_var(ncid, 'time', NF90_INT, (/time_dim_id/), time_var_id)
  ierr = nf90_put_att(ncid, time_var_id, 'long_name', 'time')
  ierr = nf90_put_att(ncid, time_var_id, 'units', 'hours since 2014-08-01 00:00:00')
  ierr = nf90_put_att(ncid, time_var_id, 'cartesian_axis', 'T')
  ierr = nf90_put_att(ncid, time_var_id, 'calendar_type', 'JULIAN')
  ierr = nf90_put_att(ncid, time_var_id, 'calendar', 'JULIAN')

  ! Define all CAs
  do k = 1, nc_out
    write(var_name, '("ca", I0)') k
    ierr = nf90_def_var(ncid, trim(var_name), NF90_FLOAT, (/xdim_id, ydim_id, time_dim_id/), ca_ids(k))
    ierr = nf90_put_att(ncid, ca_ids(k), 'long_name', 'random pattern' // trim(var_name))
    ierr = nf90_put_att(ncid, ca_ids(k), 'units', 'None')
    ierr = nf90_put_att(ncid, ca_ids(k), 'missing_value', undef)
    ierr = nf90_put_att(ncid, ca_ids(k), '_FillValue', undef)
    ierr = nf90_put_att(ncid, ca_ids(k), 'cell_methods', 'time: point')
  enddo

  ierr = nf90_enddef(ncid)

  ! Write coordinates
  ierr = nf90_put_var(ncid, xt_var_id, grid_xt)
  ierr = nf90_put_var(ncid, yt_var_id, grid_yt)

else
  ierr = nf90_open(filename, NF90_WRITE, ncid)
  ierr = nf90_inq_varid(ncid, 'time', time_var_id)
  do k=1,nc_out
    write(var_name, '("ca",I0)') k
    ierr = nf90_inq_varid(ncid, trim(var_name), ca_ids(k))
  enddo
endif

! Add new timestep
ierr = nf90_put_var(ncid, time_var_id, time, start=(/time/))

! Add all CAs
do k = 1,nc_out
  workg = transpose(real(ca(:,:,k), kind=4))
  ierr = nf90_put_var(ncid, ca_ids(k), workg, start=(/1,1,time/))
enddo

ierr = nf90_close(ncid)

deallocate(grid_xt, grid_yt, workg, ca_ids)

end subroutine fv3ca_write 

function get_strid(id, ntasks) result(str)

integer, intent(in) :: id, ntasks 
character(len=4) :: str

if (ntasks > 1000) then
  write(str, '(I4.4)') id+1
else if (ntasks> 100) then
  write(str, '(I3.3)') id+1
else if (ntasks> 10) then
  write(str, '(I2.2)') id+1
else
  write(str, '(I1.1)') id+1
endif

end function get_strid


end program standalone_ggca
