program standalone_ggca

use netcdf
use mpi_f08
use ESMF
use cellular_automata_mod
use atmosphere_stub_mod, only: Atm, atmosphere_init_stub
use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error, FATAL
use fms_mod, only: fms_init
use xgrid_mod, only: grid_box_type
use spectral_transforms, only: stochy_la2ga
use fv3_to_gaussian_remap_mod, only: remap_fv3_to_gaussian
use kinddef, only: kind_phys, kind_dbl_prec

implicit none

integer :: ntasks, nblks, blksz, ierr, my_id, root_pe, i, j, nx, ny
integer :: isc, iec, jsc, jec, nb, n
integer :: istart, iend, dump_time, nsteps
integer :: rc, local_de_count, tile
integer :: decomptile(2,6)
type(MPI_Datatype) :: grid_stats_type
type(ESMF_Grid) :: fv3_grid
type(ESMF_DistGrid) :: fv3_distgrid
type(ESMF_CoordSys_Flag) :: fv3_coord_sys
type(ESMF_TypeKind_Flag) :: fv3_coord_kind
real(ESMF_KIND_R8), pointer :: fv3_lon_ptr(:,:), fv3_lat_ptr(:,:)
real(kind_dbl_prec) :: local_grid_stats(2), global_grid_stats(2)
real(kind=4) :: ts

type(grid_box_type) :: grid_box
type(ca_grid_t) :: ca_object

real(kind=kind_dbl_prec), allocatable :: fv3_ca(:,:,:)
real(kind=kind_dbl_prec), allocatable :: gg_ca_esmf(:,:,:)  ! (nlon,nlat,nca) FV3 -> Gaussian
real(ESMF_KIND_R8), allocatable :: remapped(:,:)
real(kind=kind_dbl_prec), allocatable :: fv3_lons(:,:), fv3_lats(:,:)

character(len=128) :: output_file, filename, mosaic_file
character(len=4) :: strid

! Cellular automata controls
integer :: nca, nlives, ncells
real(kind=kind_phys) :: nfracseed, rcell
integer :: nseed, nseed_g
logical :: do_ca, ca_sgs, ca_global, ca_smooth
integer*8 :: iseed_ca, clock_count, clock_rate
integer :: nspinup, nsmooth
real(kind=kind_phys) :: l_grid, ca_amplitude
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

! All PETs must initialize the same Gaussian grid. A local mean can differ
! among tiles, so calculate one global cell-size mean before initialization.
local_grid_stats(1) = sum(real(grid_box%dx(isc:iec,jsc:jec),kind_dbl_prec)) + &
                      sum(real(grid_box%dy(isc:iec,jsc:jec),kind_dbl_prec))
local_grid_stats(2) = real(nx,kind_dbl_prec) * real(ny,kind_dbl_prec)
call MPI_Type_match_size(MPI_TYPECLASS_REAL, storage_size(0.0_kind_dbl_prec)/8, grid_stats_type, ierr)
if (ierr /= MPI_SUCCESS) call mpp_error(FATAL, 'MPI_Type_match_size failed')
call MPI_Allreduce(local_grid_stats, global_grid_stats, 2, grid_stats_type, MPI_SUM, MPI_COMM_WORLD, ierr)
if (ierr /= MPI_SUCCESS) call mpp_error(FATAL, 'MPI_Allreduce for Gaussian grid scale failed')
l_grid = real(2.5_kind_dbl_prec * global_grid_stats(1) / global_grid_stats(2), kind_phys)
print *, 'average FV3 cell size (m) = ', l_grid

if (do_ca) then
  ! Choose a clock-based seed once on the root and share it: calling the
  ! clock independently on each PET would produce different CA fields.
  if (iseed_ca <= 0) then
    if (my_id == root_pe) then
      call system_clock(count=clock_count, count_rate=clock_rate)
      if (clock_rate <= 0 .or. clock_count < 0) call mpp_error(FATAL, 'system clock unavailable for CA seed')
      iseed_ca = modulo(clock_count, int(huge(0), kind=8)) + 1_8
      print *, 'Clock-generated CA seed: ', iseed_ca
    endif
    call MPI_Bcast(iseed_ca, 1, MPI_INTEGER8, root_pe, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call mpp_error(FATAL, 'MPI_Bcast for CA seed failed')
  endif
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

  ! Read the same FV3 mosaic used by the stub. This test assumes one whole
  ! tile per PET; check the actual ESMF tile rather than relying on rank order.
  call ESMF_Initialize(rc=rc)
  call check_esmf(rc, 'ESMF_Initialize')
  write(mosaic_file, '(A,I0,A)') 'INPUT/C', Atm(1)%npx-1, '_mosaic.nc'
  decomptile = 1
  fv3_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file), tileFilePath='INPUT/', &
                                   regDecompPTile=decomptile, &
                                   staggerlocList=(/ESMF_STAGGERLOC_CENTER/), rc=rc)
  call check_esmf(rc, 'ESMF_GridCreateMosaic')
  call ESMF_GridGet(fv3_grid, localDECount=local_de_count, coordSys=fv3_coord_sys, &
                    coordTypeKind=fv3_coord_kind, rc=rc)
  call check_esmf(rc, 'ESMF_GridGet')
  if (local_de_count /= 1) call mpp_error(FATAL, 'expected one FV3 DE per PET')
  if (fv3_coord_kind /= ESMF_TYPEKIND_R8) call mpp_error(FATAL, 'expected R8 FV3 grid coordinates')
  call ESMF_GridGet(fv3_grid, staggerloc=ESMF_STAGGERLOC_CENTER, distgrid=fv3_distgrid, rc=rc)
  call check_esmf(rc, 'ESMF_GridGet center DistGrid')
  call ESMF_DistGridGet(fv3_distgrid, localDe=0, tile=tile, rc=rc)
  call check_esmf(rc, 'ESMF_DistGridGet tile')
  if (tile /= Atm(1)%tile_of_mosaic) call mpp_error(FATAL, 'FMS and ESMF FV3 tiles do not match')

  ! Use the stub's FV3 cell centers, where stochy_la2ga produced fv3_ca.
  ! The mosaic supplies the cubed-sphere topology and decomposition; its
  ! precomputed centers need not be exactly the centers calculated by FV3.
  call ESMF_GridGetCoord(fv3_grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
                         localDe=0, coordDim=1, farrayPtr=fv3_lon_ptr, rc=rc)
  call check_esmf(rc, 'ESMF_GridGetCoord FV3 longitude')
  call ESMF_GridGetCoord(fv3_grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
                         localDe=0, coordDim=2, farrayPtr=fv3_lat_ptr, rc=rc)
  call check_esmf(rc, 'ESMF_GridGetCoord FV3 latitude')
  if (any(shape(fv3_lon_ptr) /= (/nx,ny/)) .or. any(shape(fv3_lat_ptr) /= (/nx,ny/))) &
    call mpp_error(FATAL, 'FMS and ESMF FV3 grid shapes do not match')
  if (fv3_coord_sys == ESMF_COORDSYS_SPH_DEG) then
    fv3_lon_ptr = transpose(fv3_lons)
    fv3_lat_ptr = transpose(fv3_lats)
  else if (fv3_coord_sys == ESMF_COORDSYS_SPH_RAD) then
    fv3_lon_ptr = transpose(fv3_lons) * acos(-1.0_kind_dbl_prec) / 180.0_kind_dbl_prec
    fv3_lat_ptr = transpose(fv3_lats) * acos(-1.0_kind_dbl_prec) / 180.0_kind_dbl_prec
  else
    call mpp_error(FATAL, 'FV3 ESMF grid must use spherical coordinates')
  endif
  allocate(gg_ca_esmf(ca_object%nlon, ca_object%nlat, max(1,ca_object%nca)))

  strid = get_strid(my_id, ntasks)
  
  do i = istart, iend 
    print *, "kstep = ", i
  
    call update_cellular_automata(i, iseed_ca, nfracseed, nseed, nspinup, nlives, ca_amplitude, ca_object)
    call remap_ca_to_fv3(ca_object, fv3_lons, fv3_lats, fv3_ca)

    if (mod(i,5) == 0) then
      ! Round trip: the existing Gaussian -> FV3 interpolation above followed
      ! by the new FV3 -> Gaussian ESMF remap. All PETs must participate.
      do n = 1, ca_object%nca
        call remap_fv3_to_gaussian(fv3_grid, real(transpose(fv3_ca(:,:,n)), kind_phys), &
                                   ca_object, remapped, MPI_COMM_WORLD)
        gg_ca_esmf(:,:,n) = remapped
      enddo
      write(filename, '("RESTART/ca_out.tile", A, ".nc")') trim(strid)
      call fv3ca_write(trim(filename), fv3_ca, i, ierr)
      if (my_id == root_pe) then
        call fv3ca_write("RESTART/gg_out.nc", transpose3(ca_object%field), i, ierr, ca_object%lons, ca_object%lats)
        call fv3ca_write("RESTART/gg_esmf_out.nc", transpose3(gg_ca_esmf), i, ierr, ca_object%lons, ca_object%lats)
      endif
      if (write_restart .and. my_id == root_pe) then
        print *, "writing restart: my_id",my_id,"time",i
        write(filename, '("RESTART/ca_restart.nc")')
        call write_cellular_automata(trim(filename), ca_object, i)
      endif
    endif
  enddo

  call ESMF_GridDestroy(fv3_grid, rc=rc)
  call check_esmf(rc, 'ESMF_GridDestroy')
  call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=rc)
  call check_esmf(rc, 'ESMF_Finalize')
endif
contains

! fv3ca_write expects (ny, nx, nca) and transposes internally;
! flip (nx, ny, nca)-ordered arrays to match
function transpose3(a) result(b)
  real(kind=kind_dbl_prec), intent(in) :: a(:,:,:)
  real(kind=kind_dbl_prec) :: b(size(a,2), size(a,1), size(a,3))
  integer :: k

  do k = 1, size(a,3)
    b(:,:,k) = transpose(a(:,:,k))
  enddo
end function transpose3

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

subroutine fv3ca_write(filename, ca, time, ierr, lons, lats)
  use netcdf
  use kinddef, only : kind_phys
  ! In/Out variables
  character(len=*), intent(in)         :: filename
  real(kind=kind_dbl_prec), intent(in) :: ca(:,:,:)
  integer, intent(in)             :: time
  integer, intent(out)                 :: ierr
  real(kind=kind_dbl_prec), optional, intent(in) :: lons(:), lats(:)

  ! Local Variables
  integer :: ncid, xdim_id, ydim_id, time_dim_id
  integer :: xt_var_id, yt_var_id, time_var_id
  integer :: nx_out, ny_out, nc_out, k, ix, iy, record, existing_nx, existing_ny
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
  grid_xt = [(real(ix, kind=4), ix=1,nx_out)]
  grid_yt = [(real(iy, kind=4), iy=1,ny_out)]
  if (present(lons) .neqv. present(lats)) call mpp_error(FATAL, 'both output coordinates are required')
  if (present(lons)) then
    if (size(lons) /= nx_out .or. size(lats) /= ny_out) call mpp_error(FATAL, 'output coordinate size mismatch')
    grid_xt = real(lons, kind=4)
    grid_yt = real(lats, kind=4)
  endif

  inquire(file=trim(filename), exist=file_exists)

if (.not. file_exists) then
  call check_nc(nf90_create(filename, NF90_CLOBBER, ncid), 'creating '//trim(filename))
  record = 1

  ! Define Dimensions
  call check_nc(nf90_def_dim(ncid, 'grid_xt', nx_out, xdim_id), 'defining grid_xt')
  call check_nc(nf90_def_dim(ncid, 'grid_yt', ny_out, ydim_id), 'defining grid_yt')
  call check_nc(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim_id), 'defining time')

  ! Define Coordinates
  call check_nc(nf90_def_var(ncid, 'grid_xt', NF90_FLOAT, (/xdim_id/), xt_var_id), 'defining x coordinate')
  call check_nc(nf90_def_var(ncid, 'grid_yt', NF90_FLOAT, (/ydim_id/), yt_var_id), 'defining y coordinate')
  if (present(lons)) then
    call check_nc(nf90_put_att(ncid, xt_var_id, 'units', 'degrees_E'), 'setting longitude units')
    call check_nc(nf90_put_att(ncid, yt_var_id, 'units', 'degrees_N'), 'setting latitude units')
  else
    call check_nc(nf90_put_att(ncid, xt_var_id, 'long_name', 'FV3 local i index'), 'setting x name')
    call check_nc(nf90_put_att(ncid, yt_var_id, 'long_name', 'FV3 local j index'), 'setting y name')
  endif
  call check_nc(nf90_put_att(ncid, xt_var_id, 'cartesian_axis', 'X'), 'setting x axis')
  call check_nc(nf90_put_att(ncid, yt_var_id, 'cartesian_axis', 'Y'), 'setting y axis')
  call check_nc(nf90_def_var(ncid, 'time', NF90_INT, (/time_dim_id/), time_var_id), 'defining time variable')
  call check_nc(nf90_put_att(ncid, time_var_id, 'long_name', 'model step'), 'setting time name')
  call check_nc(nf90_put_att(ncid, time_var_id, 'units', '1'), 'setting time units')

  ! Define all CAs
  do k = 1, nc_out
    write(var_name, '("ca", I0)') k
    call check_nc(nf90_def_var(ncid, trim(var_name), NF90_FLOAT, (/xdim_id, ydim_id, time_dim_id/), &
                               ca_ids(k)), 'defining '//trim(var_name))
    call check_nc(nf90_put_att(ncid, ca_ids(k), 'long_name', 'random pattern' // trim(var_name)), 'setting name')
    call check_nc(nf90_put_att(ncid, ca_ids(k), 'units', 'None'), 'setting units')
    call check_nc(nf90_put_att(ncid, ca_ids(k), 'missing_value', undef), 'setting missing value')
    call check_nc(nf90_put_att(ncid, ca_ids(k), '_FillValue', undef), 'setting fill value')
    call check_nc(nf90_put_att(ncid, ca_ids(k), 'cell_methods', 'time: point'), 'setting cell methods')
  enddo

  call check_nc(nf90_enddef(ncid), 'ending NetCDF definition')

  ! Write coordinates
  call check_nc(nf90_put_var(ncid, xt_var_id, grid_xt), 'writing x coordinate')
  call check_nc(nf90_put_var(ncid, yt_var_id, grid_yt), 'writing y coordinate')

else
  call check_nc(nf90_open(filename, NF90_WRITE, ncid), 'opening '//trim(filename))
  call check_nc(nf90_inq_dimid(ncid, 'grid_xt', xdim_id), 'finding x dimension')
  call check_nc(nf90_inq_dimid(ncid, 'grid_yt', ydim_id), 'finding y dimension')
  call check_nc(nf90_inquire_dimension(ncid, xdim_id, len=existing_nx), 'reading x dimension')
  call check_nc(nf90_inquire_dimension(ncid, ydim_id, len=existing_ny), 'reading y dimension')
  if (existing_nx /= nx_out .or. existing_ny /= ny_out) &
    call mpp_error(FATAL, 'existing output dimensions differ: '//trim(filename))
  call check_nc(nf90_inq_varid(ncid, 'grid_xt', xt_var_id), 'finding x coordinate')
  call check_nc(nf90_inq_varid(ncid, 'grid_yt', yt_var_id), 'finding y coordinate')
  block
    real(kind=4) :: old_x(nx_out), old_y(ny_out)
    call check_nc(nf90_get_var(ncid, xt_var_id, old_x), 'reading x coordinate')
    call check_nc(nf90_get_var(ncid, yt_var_id, old_y), 'reading y coordinate')
    if (any(old_x /= grid_xt) .or. any(old_y /= grid_yt)) &
      call mpp_error(FATAL, 'existing output coordinates differ: '//trim(filename))
  end block
  call check_nc(nf90_inq_dimid(ncid, 'time', time_dim_id), 'finding time dimension')
  call check_nc(nf90_inquire_dimension(ncid, time_dim_id, len=record), 'reading time dimension')
  record = record + 1
  call check_nc(nf90_inq_varid(ncid, 'time', time_var_id), 'finding time variable')
  do k=1,nc_out
    write(var_name, '("ca",I0)') k
    call check_nc(nf90_inq_varid(ncid, trim(var_name), ca_ids(k)), 'finding '//trim(var_name))
  enddo
endif

! NetCDF records start at one; model step zero is a value, not an index.
call check_nc(nf90_put_var(ncid, time_var_id, time, start=(/record/)), 'writing model step')

! Add all CAs
do k = 1,nc_out
  workg = transpose(real(ca(:,:,k), kind=4))
  write(var_name, '("ca",I0)') k
  call check_nc(nf90_put_var(ncid, ca_ids(k), workg, start=(/1,1,record/)), 'writing '//trim(var_name))
enddo

call check_nc(nf90_close(ncid), 'closing '//trim(filename))
ierr = nf90_noerr

deallocate(grid_xt, grid_yt, workg, ca_ids)

end subroutine fv3ca_write 

subroutine check_nc(status, operation)
  integer, intent(in) :: status
  character(len=*), intent(in) :: operation
  if (status /= nf90_noerr) then
    write(0,*) trim(operation), ': ', trim(nf90_strerror(status))
    call mpp_error(FATAL, 'NetCDF output failed')
  endif
end subroutine check_nc

subroutine check_esmf(status, operation)
  integer, intent(in) :: status
  character(len=*), intent(in) :: operation
  if (status /= ESMF_SUCCESS) then
    write(0,*) trim(operation), ': ESMF rc = ', status
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif
end subroutine check_esmf

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

