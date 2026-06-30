module cellular_automata_gg_mod

use mpp_mod, only : mpp_pe, mpp_root_pe
use kinddef, only : kind_io8, kind_phys, kind_qdt_prec
use constants_mod, only : radius
use spectral_transforms, only : stochy_la2ga
use netcdf

implicit none
private

public :: ggca_grid_t
public :: ggca_init
public :: ggca_write
public :: cellular_automata_gg

! Gaussian grid cellular autamata
type :: ggca_grid_t 
    integer :: ntrunc = 0
    integer :: nlat = 0
    integer :: nlon = 0
    real(kind=kind_io8) :: wlon = 0.0_kind_io8
    real(kind=kind_io8) :: rnlat = 0.0_kind_io8
    real(kind=kind_io8), allocatable :: lats(:)
    real(kind=kind_io8), allocatable :: lons(:)
end type ggca_grid_t 

contains

    subroutine cellular_automata_gg(         &
                    kstep,                   &
                    restart,                 &
                    first_time_step,         &
                    ca1_cpl,ca2_cpl,ca3_cpl, & ! Output coarse-grained to FV3
                    nca,                     &
                    ncells,                  &
                    nlives,                  &
                    nfracseed,               &
                    nseed,                   &
                    iseed_ca,                &
                    ca_smooth,               &
                    nspinup,                 &
                    nsmooth,                 &
                    ca_amplitude )
    
    use kinddef,           only: kind_dbl_prec, kind_phys
    use update_ca,         only: update_cells_gg
    use random_numbers,    only: random_01_CB
    use stochy_internal_state_mod, only: stochy_internal_state
    
    implicit none
    
    !L.Bengtsson, 2017-06
    !P.Pegion, 2021-09
    ! swtich to new random number generator and improve computational efficiency
    ! and remove unsued code. Also add restart capability ca_global
    !
    !K.Barton, 2026-06 : added gaussian grid
    
    !This program evolves a cellular automaton uniform over the globe 
    
    integer, intent(in) :: kstep,ncells,nca,nlives,nseed,nspinup,nsmooth
    integer(kind=kind_dbl_prec), intent(in) :: iseed_ca
    real(kind=kind_phys), intent(in) :: nfracseed,ca_amplitude
    logical, intent(in) :: ca_smooth,first_time_step, restart
    real(kind=kind_phys), intent(out) :: ca1_cpl(:,:),ca2_cpl(:,:),ca3_cpl(:,:)
    integer, parameter :: nxc=10,nyc=10
    integer :: i,j,k,nf, my_pe, root_pe
    integer :: iini_g(nxc,nyc,nca), ilives_g(nxc,nyc)
   
    ! For GG initialization
    type(ggca_grid_t)                     :: ggrid
    real(kind=kind_dbl_prec)              :: wlon, rnlat
    real(kind=kind_dbl_prec), parameter   :: l_min = 10000000
    real(kind=kind_dbl_prec), allocatable :: gg_lons(:), gg_lats(:), ca_field(:, :)
    integer                               :: ntrunc

    my_pe = mpp_pe()
    root_pe = mpp_root_pe()
    
    
    ! Initialize gaussian grid
    print *, "Initializing grid"
    call ggca_init(l_min, ggrid)

    print *, "Allocating CA field"
    allocate(ca_field(ggrid%nlon, ggrid%nlat))

    print *, "Populating CA field"
    do i=1,ggrid%nlon
        do j=1,ggrid%nlat
            ca_field(i,j) = i+j
        end do 
    end do

    print *, "Outputting CA field"
    call ggca_write("ggca_write.nc", ggrid, ca_field, my_pe, root_pe)

!    do i=1,nxc
!        do j=1,nyc
!        ilives_g(i,j)=i*j
!            do k=1,nca
!                iini_g(i,j,k)=i*j+k
!            end do
!        end do
!    end do
!    
!    print *, "call update_cells_gg"
!    
!    ! Update cells does the following:
!    ! Allocate board, lives, and board_halo
!    ! First time step: board&lives = iini&ilives
!    ! Seed cells with new active cells each nseed time-step
!    ! Perform spinup (?)
!    ! Evolve CA
!    ! Coarse-grain back to NWP Grid
!    
!    nf = 1
!    call update_cells_gg(       &
!            kstep,              &
!            first_time_step,    &
!            iseed_ca,           &
!            restart,            &
!            nca,                &
!            nxc,nyc,            &
!            ca1_cpl,            &
!            iini_g,             &
!            ilives_g,           &
!            nf)
!    
    
    end subroutine cellular_automata_gg
    
    subroutine ggca_init(l_min, grid)
        ! Input / Output
        real(kind=kind_phys), intent(in) :: l_min
        type(ggca_grid_t), intent(inout) :: grid

        ! Local
        real(kind=kind_io8) :: circ 

        circ = 2.0_kind_io8 * acos(-1.0_kind_io8) * real(radius, kind_io8)

        ! Generate truncation value from minimum resolved length scale
        grid%ntrunc = int(circ / l_min)
        grid%ntrunc = ((grid%ntrunc + 1)/4)*4 + 2

        ! Check logic here
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
    
    subroutine ggca_write(outfile, grid, field, mype, root_pe)
        ! In/Out
        character(len=*), intent(in) :: outfile
        type(ggca_grid_t), intent(in) :: grid
        real(kind=kind_io8), intent(in) :: field(grid%nlon, grid%nlat)
        integer, intent(in), optional :: mype, root_pe

        ! Local
        integer :: ncid, dlon, dlat, vlon, vlat, vfield, ierr

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

        ierr = nf90_create(trim(outfile), NF90_CLOBBER, ncid)
        call check_nf90(ierr, 'nf90_create('//trim(outfile)//')')
        ierr = nf90_def_dim(ncid, "lon", grid%nlon, dlon)
        call check_nf90(ierr, 'nf90_def_dim(lon)')
        ierr = nf90_def_dim(ncid, "lat", grid%nlat, dlat)
        call check_nf90(ierr, 'nf90_def_dim(lat)')

        ierr = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/dlon/), vlon)
        call check_nf90(ierr, 'nf90_def_var(lon)')
        ierr = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/dlat/), vlat)
        call check_nf90(ierr, 'nf90_def_var(lat)')
        ierr = nf90_def_var(ncid, "field", NF90_DOUBLE, (/dlon, dlat/), vfield)
        call check_nf90(ierr, 'nf90_def_var(field)')

        ierr = nf90_enddef(ncid)
        call check_nf90(ierr, 'nf90_enddef')

        ierr = nf90_put_var(ncid, vlon, grid%lons)
        call check_nf90(ierr, 'nf90_put_var(lon)')
        ierr = nf90_put_var(ncid, vlat, grid%lats)
        call check_nf90(ierr, 'nf90_put_var(lat)')
        ierr = nf90_put_var(ncid, vfield, field)
        call check_nf90(ierr, 'nf90_put_var(field)')

        ierr = nf90_close(ncid)
        call check_nf90(ierr, 'nf90_close')
    end subroutine ggca_write 

    subroutine build_gaussian_latitudes(nlat, lats)
        integer, intent(in) :: nlat
        real(kind=kind_io8), intent(out) :: lats(nlat)

        integer :: k, lgghaf
        real(kind=kind_io8), allocatable :: colrad(:), wgt(:), rcs2(:)
        real(kind=kind_io8), parameter   :: rad2deg = 180.0_kind_io8 / acos(-1.0_kind_io8)

        integer :: iter, k1, l2
        real(kind=kind_qdt_prec) :: drad, dradz, p1, p2, pi, rad, rc
        real(kind=kind_qdt_prec) :: rl2, scale, si, w, x

        real(kind=kind_io8), parameter :: cons1 = 1.0_kind_io8
        real(kind=kind_io8), parameter :: cons2 = 2.0_kind_io8
        real(kind=kind_io8), parameter :: cons4 = 4.0_kind_io8
        real(kind=kind_io8), parameter :: cons0p25 = 0.25_kind_io8

#ifdef NO_QUAD_PRECISION
        real(kind=kind_qdt_prec), parameter :: eps = 1.0e-12_kind_qdt_prec
#else
        real(kind=kind_qdt_prec), parameter :: eps = 1.0e-20_kind_qdt_prec
#endif
        
        print *, "building gaussain lats"

        lgghaf = nlat / 2
        allocate(colrad(lgghaf), wgt(lgghaf), rcs2(lgghaf))

        ! Based on glats_stochy
        si = 1.0_kind_qdt_prec
        l2 = 2*lgghaf
        rl2 = real(l2, kind_qdt_prec)
        scale = 2.0_kind_qdt_prec / (rl2 * rl2)
        k1 = l2 - 1
        pi = atan(si) * real(cons4, kind_qdt_prec)

        dradz = pi / real(lgghaf, kind_qdt_prec) / 200.0_kind_qdt_prec
        rad = 0.0_kind_qdt_prec

        do k = 1, lgghaf
            iter = 0
            drad = dradz
1           call poly(l2, rad, p2)
2           p1 = p2
            rad = rad + drad
            call poly(l2, rad, p2)
            if (sign(si, p1) == sign(si, p2)) go to 2
            if (drad < eps) go to 3
            rad = rad - drad
            drad = drad * real(cons0p25, kind_qdt_prec)
            go to 1
3           continue

            colrad(k) = real(rad, kind_io8)
            call poly(k1, rad, p1)
            x = cos(rad)
            w = scale * (1.0_kind_qdt_prec - x*x) / (p1*p1)
            wgt(k) = real(w, kind_io8)
            rc = 1.0_kind_qdt_prec / (sin(rad) * sin(rad))
            rcs2(k) = real(rc, kind_io8)
        end do
        ! end glats_stochy 

        do k = 1, lgghaf
            lats(k) = -1.0_kind_io8 * colrad(lgghaf - k + 1) * rad2deg
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

    subroutine check_nf90(status, where)
        use netcdf
        integer, intent(in) :: status
        character(len=*), intent(in) :: where
        if (status /= nf90_noerr) then
          write(0,*) 'NetCDF error in ', trim(where), ': ', trim(nf90_strerror(status))
          stop 2
        end if
    end subroutine check_nf90


end module cellular_automata_gg_mod
