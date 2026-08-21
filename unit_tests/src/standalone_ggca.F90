program standalone_ggca

use netcdf
use mpi_f08
use cellular_automata_gg_mod, only : cellular_automata_gg, ggca_grid_t, ggca_write
use update_ca, only: read_ca_restart
use atmosphere_stub_mod, only: Atm, atmosphere_init_stub
use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe
use fms_mod, only: fms_init
use xgrid_mod, only: grid_box_type
use spectral_transforms, only: stochy_la2ga
use kinddef, only: kind_phys, kind_dbl_prec

implicit none
integer :: ntasks, nblks, blksz, ierr, my_id, root_pe, i, j, nx, ny
integer :: isc, iec, jsc, jec, nb, npts, i1, i2
integer :: istart, dump_time, ct
real(kind=4) :: ts
logical :: first_time_step
type(grid_box_type) :: grid_box
type(ggca_grid_t) :: gaussian_grid
real(kind=kind_dbl_prec), allocatable :: gaussian_ca(:,:,:)
real(kind=kind_dbl_prec), allocatable :: fv3_ca(:,:,:)
real(kind=kind_dbl_prec), allocatable :: fv3_lons(:,:), fv3_lats(:,:)
character(len=128) :: output_file
character(len=4) :: strid

! Cellular automata controls
integer :: nca, nlives, ncells
integer :: nca_g, nlives_g, ncells_g
real(kind=kind_phys) :: nfracseed, rcell
integer :: nseed, nseed_g
logical :: do_ca, ca_sgs, ca_global, ca_smooth
integer*8 :: iseed_ca
integer :: nspinup, nsmooth
real :: ca_amplitude
real(kind=kind_phys) :: l_grid
logical :: ca_closure, ca_entr, ca_trigger, warm_start
real(kind=kind_phys), allocatable :: ca1(:,:), ca2(:,:), ca3(:,:)

! Some SGS fields
real(kind=kind_phys), allocatable :: cond_in(:,:), condition(:,:)
real(kind=kind_phys), allocatable :: sst(:,:), lmsk(:,:), lake(:,:)
real(kind=kind_phys), allocatable :: ca_deep(:,:), ca_turb(:,:), ca_shal(:,:)

namelist /gfs_physics_nml/ do_ca, ca_sgs, ca_global, nca, ncells, nlives, nseed,       &
                          nfracseed, rcell, ca_trigger, ca_entr, ca_closure, nca_g,    &
                          ncells_g, nlives_g, nseed_g, ca_smooth, nspinup, iseed_ca,   &
                          nsmooth, ca_amplitude, warm_start

! Defaults
first_time_step = .true.
warm_start = .false.
nca            = 0
nca_g          = 0
ncells_g       = 1
nlives_g       = 1
nfracseed      = 0.5
nseed          = 100000
iseed_ca       = 0
nspinup        = 1
do_ca          = .false.
ca_sgs         = .false.
ca_global      = .false.
ca_smooth      = .false.
ca_amplitude   = 500.
rcell          = 0.0

! open namelist file
open (unit=565, file='input.nml', status='OLD', iostat=ierr)
read(565,gfs_physics_nml)
close(565)
! define stuff
print *, 'ca_sgs,ca_global',ca_sgs,ca_global
if (.not. ca_sgs) then
   nca=0
endif
if (.not. ca_global) then
   nca_g=0
endif

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

if (ca_global) then
  allocate(ca1 (nblks,blksz))
  allocate(ca2 (nblks,blksz))
  allocate(ca3 (nblks,blksz))
  allocate(fv3_ca(nblks,blksz,max(1,nca_g)))
  allocate(fv3_lons(nblks,blksz))
  allocate(fv3_lats(nblks,blksz))

  do j=jsc,jec
    do i=isc,iec
      nb = j-jsc+1
      fv3_lons(nb,i-isc+1) = Atm(1)%gridstruct%agrid_64(i,j,1) * 180.0_kind_dbl_prec / acos(-1.0_kind_dbl_prec)
      fv3_lats(nb,i-isc+1) = Atm(1)%gridstruct%agrid_64(i,j,2) * 180.0_kind_dbl_prec / acos(-1.0_kind_dbl_prec)
    enddo
  enddo
endif

if (ca_sgs) then
  allocate(ca_deep  (nblks,blksz))
  allocate(ca_turb  (nblks,blksz))
  allocate(ca_shal  (nblks,blksz))
  allocate(sst      (nblks,blksz))
  allocate(lmsk     (nblks,blksz))
  allocate(lake     (nblks,blksz))
  allocate(condition(nblks,blksz))
!  allocate(dx       (nblks,blksz))
!  allocate(uwind    (nblks,blksz, levs))
!  allocate(vwind    (nblks,blksz, levs))
!  allocate(height   (nblks,blksz, levs))
  allocate(cond_in(isc:iec,jsc:jec))
  
  sst (:,:) = 303.
  lmsk(:,:) = 0.
  lake(:,:) = 0
 i1=isc
 j=jsc
 do nb=1,nblks
   i2=i1+blksz-1
   if (i2 .le. iec) then  
     condition(nb,1:blksz) = cond_in(i1:i2,j)
     i1=i1+blksz
   else
     npts=iec-i1+1
     condition(nb,1:npts) = cond_in(i1:iec,j)
     if (j.LT. jec) then
       condition(nb,npts+1:blksz) = cond_in(isc:isc+(blksz-npts+1),j+1) 
     endif
     i1=npts+1
     j=j+1
   endif
   if (i2.EQ.iec) then
     i1=isc
     j=j+1
   endif
 end do
endif

dump_time=50
if (warm_start) then
  istart=dump_time+1
  call read_ca_restart(Atm(1)%domain,0,ncells,nca,ncells_g,nca_g)
else
  istart=1
endif
ct = 1

do i = istart, 101
  print *, "kstep = ", i
  if (ca_global) then
    call cellular_automata_gg( &
           i,                  &
           warm_start,         &
           first_time_step,    &
           nca_g,              &
           ncells_g,           & ! Passes into l_scale
           nlives_g,           &
           nfracseed,          &
           nseed_g,            &
           iseed_ca,           &
           ca_smooth,          &
           nspinup,            &
           nsmooth,            &
           ca_amplitude,       &
           l_grid,             &
           gaussian_ca,        &
           gaussian_grid       &
         )

    call remap_ca_to_fv3(gaussian_ca, gaussian_grid, fv3_lons, fv3_lats, fv3_ca)

    ca1(:,:) = 0.0_kind_phys
    ca2(:,:) = 0.0_kind_phys
    ca3(:,:) = 0.0_kind_phys
    if (nca_g >= 1) ca1(:,:) = real(fv3_ca(:,:,1), kind=kind_phys)
    if (nca_g >= 2) ca2(:,:) = real(fv3_ca(:,:,2), kind=kind_phys)
    if (nca_g >= 3) ca3(:,:) = real(fv3_ca(:,:,3), kind=kind_phys)

    do j=1,nca_g
      write(output_file, '("ca",I1,"_gaussian.nc")') j
      call ggca_write(trim(output_file), gaussian_grid, gaussian_ca(:,:,j), i, my_id, root_pe)
    enddo

    if (ntasks > 1000) then
      write(strid,'(I4.4)') my_id+1
    else if (ntasks > 100) then
      write(strid,'(I3.3)') my_id+1
    else if (ntasks > 10) then
      write(strid,'(I2.2)') my_id+1
    else
      write(strid,'(I1.1)') my_id+1
    endif
    ts = i / 4.0
    if (mod(i-1,5) == 0) then
      call fv3ca_write('ca_out.tile'//trim(strid)//'.nc', ca1, ca2, ca3, ts, ct, ierr)
    endif

  endif

  first_time_step=.false.

enddo

contains

subroutine remap_ca_to_fv3(gaussian_fields, gaussian_grid, target_lons, target_lats, fv3_fields)
  real(kind=kind_dbl_prec), intent(in) :: gaussian_fields(:,:,:)
  type(ggca_grid_t), intent(in) :: gaussian_grid
  real(kind=kind_dbl_prec), intent(in) :: target_lons(:,:), target_lats(:,:)
  real(kind=kind_dbl_prec), intent(out) :: fv3_fields(:,:,:)

  integer :: block, field_number, nblocks, blocksize

  nblocks = size(target_lons, 1)
  blocksize = size(target_lons, 2)
  do field_number=1,size(gaussian_fields, 3)
    do block=1,nblocks
      call stochy_la2ga(                                  &
             gaussian_fields(:,:,field_number),           &
             gaussian_grid%nlon, gaussian_grid%nlat,      &
             gaussian_grid%lons, gaussian_grid%lats,      &
             gaussian_grid%wlon, gaussian_grid%rnlat,     &
             fv3_fields(block,1:blocksize, field_number), &
             blocksize,                                   &
             target_lats(block,:), target_lons(block,:)   &
           )
    enddo
  enddo
end subroutine remap_ca_to_fv3

subroutine fv3ca_write(filename, ca1, ca2, ca3, time, record, ierr)
  use netcdf
  use kinddef, only : kind_phys
  character(len=*), intent(in) :: filename
  real(kind=kind_dbl_prec), intent(in) :: ca1(:,:), ca2(:,:), ca3(:,:)
  real(kind=4), intent(in) :: time
  integer, intent(inout) :: record
  integer, intent(out) :: ierr

  integer :: ncid, xdim_id, ydim_id, time_dim_id
  integer :: xt_var_id, yt_var_id, time_var_id
  integer :: ca1_id, ca2_id, ca3_id
  integer :: nx_out, ny_out
  real(kind=4), allocatable :: grid_xt(:), grid_yt(:), workg(:,:)
  real(kind=4) :: undef
  logical :: file_exists

  nx_out = size(ca1, 2)
  ny_out = size(ca1, 1)
  undef = 9.99e+20
  inquire(file=trim(filename), exist=file_exists)
  allocate(grid_xt(nx_out), grid_yt(ny_out), workg(nx_out, ny_out))
  grid_xt = [(real(i, kind=kind_phys), i=1,nx_out)]
  grid_yt = [(real(j, kind=kind_phys), j=1,ny_out)]

  if (record == 1 .or. .not. file_exists) then
    ierr = nf90_create(filename, NF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, 'grid_xt', nx_out, xdim_id)
    ierr = nf90_def_dim(ncid, 'grid_yt', ny_out, ydim_id)
    ierr = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim_id)
    ierr = nf90_def_var(ncid, 'grid_xt', NF90_FLOAT, (/xdim_id/), xt_var_id)
    ierr = nf90_put_att(ncid, xt_var_id, 'long_name', 'T-cell longitude')
    ierr = nf90_put_att(ncid, xt_var_id, 'cartesian_axis', 'X')
    ierr = nf90_put_att(ncid, xt_var_id, 'units', 'degrees_E')
    ierr = nf90_def_var(ncid, 'grid_yt', NF90_FLOAT, (/ydim_id/), yt_var_id)
    ierr = nf90_put_att(ncid, yt_var_id, 'long_name', 'T-cell latitude')
    ierr = nf90_put_att(ncid, yt_var_id, 'cartesian_axis', 'Y')
    ierr = nf90_put_att(ncid, yt_var_id, 'units', 'degrees_N')
    ierr = nf90_def_var(ncid, 'time', NF90_FLOAT, (/time_dim_id/), time_var_id)
    ierr = nf90_put_att(ncid, time_var_id, 'long_name', 'time')
    ierr = nf90_put_att(ncid, time_var_id, 'units', 'hours since 2014-08-01 00:00:00')
    ierr = nf90_put_att(ncid, time_var_id, 'cartesian_axis', 'T')
    ierr = nf90_put_att(ncid, time_var_id, 'calendar_type', 'JULIAN')
    ierr = nf90_put_att(ncid, time_var_id, 'calendar', 'JULIAN')
    ierr = nf90_def_var(ncid, 'ca1', NF90_FLOAT, (/xdim_id, ydim_id, time_dim_id/), ca1_id)
    ierr = nf90_put_att(ncid, ca1_id, 'long_name', 'random pattern')
    ierr = nf90_put_att(ncid, ca1_id, 'units', 'None')
    ierr = nf90_put_att(ncid, ca1_id, 'missing_value', undef)
    ierr = nf90_put_att(ncid, ca1_id, '_FillValue', undef)
    ierr = nf90_put_att(ncid, ca1_id, 'cell_methods', 'time: point')
    ierr = nf90_def_var(ncid, 'ca2', NF90_FLOAT, (/xdim_id, ydim_id, time_dim_id/), ca2_id)
    ierr = nf90_put_att(ncid, ca2_id, 'long_name', 'random pattern')
    ierr = nf90_put_att(ncid, ca2_id, 'units', 'None')
    ierr = nf90_put_att(ncid, ca2_id, 'missing_value', undef)
    ierr = nf90_put_att(ncid, ca2_id, '_FillValue', undef)
    ierr = nf90_put_att(ncid, ca2_id, 'cell_methods', 'time: point')
    ierr = nf90_def_var(ncid, 'ca3', NF90_FLOAT, (/xdim_id, ydim_id, time_dim_id/), ca3_id)
    ierr = nf90_put_att(ncid, ca3_id, 'long_name', 'random pattern')
    ierr = nf90_put_att(ncid, ca3_id, 'units', 'None')
    ierr = nf90_put_att(ncid, ca3_id, 'missing_value', undef)
    ierr = nf90_put_att(ncid, ca3_id, '_FillValue', undef)
    ierr = nf90_put_att(ncid, ca3_id, 'cell_methods', 'time: point')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, xt_var_id, grid_xt)
    ierr = nf90_put_var(ncid, yt_var_id, grid_yt)
  else
    ierr = nf90_open(filename, NF90_WRITE, ncid)
    ierr = nf90_inq_varid(ncid, 'time', time_var_id)
    ierr = nf90_inq_varid(ncid, 'ca1', ca1_id)
    ierr = nf90_inq_varid(ncid, 'ca2', ca2_id)
    ierr = nf90_inq_varid(ncid, 'ca3', ca3_id)
  endif

  ierr = nf90_put_var(ncid, time_var_id, time, start=(/record/))
  workg = transpose(ca1)
  ierr = nf90_put_var(ncid, ca1_id, workg, start=(/1,1,record/))
  workg = transpose(ca2)
  ierr = nf90_put_var(ncid, ca2_id, workg, start=(/1,1,record/))
  workg = transpose(ca3)
  ierr = nf90_put_var(ncid, ca3_id, workg, start=(/1,1,record/))
  ierr = nf90_close(ncid)
  record = record + 1

end subroutine fv3ca_write 

end program standalone_ggca
