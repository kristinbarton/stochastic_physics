! Remapping module for Fv3 to Gaussian Grid remapping

module fv3_to_gaussian_remap_mod
  use ESMF
  use mpi_f08
  use kinddef, only: kind_phys
  use cellular_automata_mod, only: ca_grid_t

  implicit none
  private
  public :: remap_fv3_to_gaussian

contains

  subroutine remap_fv3_to_gaussian(fv3_grid, fv3_data, gaussian_grid, gaussian_data, comm)
    type(ESMF_Grid), intent(in) :: fv3_grid
    real(kind_phys), intent(in) :: fv3_data(:,:)
    type(ca_grid_t), intent(in) :: gaussian_grid
    type(MPI_Comm),  intent(in) :: comm
    real(ESMF_KIND_R8), allocatable, intent(out) :: gaussian_data(:,:) 

    type(ESMF_Grid) :: dst_grid
    type(ESMF_Field) :: src_field, dst_field
    type(ESMF_RouteHandle) :: route
    real(ESMF_KIND_R8), pointer :: src_ptr(:,:), dst_ptr(:,:)
    real(ESMF_KIND_R8), pointer :: lon_ptr(:,:), lat_ptr(:,:)
    integer :: rc, mpierr, rank
    integer :: dims(2), min_dims(2), max_dims(2)
    integer :: local_de_count, ownership(2), global_ownership(2), owner
    integer :: lower(2), upper(2), i, j

    ! Retrieve the communicator
    call MPI_Comm_rank(comm, rank, mpierr)
    call check_mpi(mpierr, 'MPI_Comm_rank')

    ! Get gaussian grid dimensions
    dims = (/gaussian_grid%nlon, gaussian_grid%nlat/)

    ! Create FV3 grid's ESMF Field
    src_field = ESMF_FieldCreate(fv3_grid, typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 name='fv3_scalar', rc=rc)
    call check_esmf(rc, 'creating the FV3 field')
    call ESMF_FieldGet(src_field, localDeCount=local_de_count, rc=rc)
    call check_esmf(rc, 'querying FV3 local DEs')
    call ESMF_FieldGet(src_field, localDe=0, farrayPtr=src_ptr, rc=rc)
    call check_esmf(rc, 'getting FV3 field storage')
    if (any(shape(src_ptr) /= shape(fv3_data))) call fail('FV3 data does not match local ESMF field shape')
    src_ptr = real(fv3_data, ESMF_KIND_R8)

    ! GG must have periodic longitudes
    dst_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), maxIndex=dims, &
                                       regDecomp=(/1,1/), indexflag=ESMF_INDEX_GLOBAL, &
                                       coordSys=ESMF_COORDSYS_SPH_DEG, &
                                       coordTypeKind=ESMF_TYPEKIND_R8, &
                                       name='gaussian_grid', rc=rc)
    call check_esmf(rc, 'creating the gaussian grid')
    call ESMF_GridAddCoord(dst_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check_esmf(rc, 'allocationg gaussian coords')
    call ESMF_GridGet(dst_grid, localDECount=local_de_count, rc=rc)
    call check_esmf(rc, 'querying gaussian local DEs')

    ! Initially remap on single process --> send out after remap
    ownership = (/local_de_count, local_de_count* (rank+1)/)
    call MPI_Allreduce(ownership, global_ownership, 2, MPI_INTEGER, MPI_SUM, comm, mpierr)
    call check_mpi(mpierr, 'finding gaussian DE owner')
    if (global_ownership(1) /= 1) call fail('Gaussian grid must have exactly one DE')
    owner = global_ownership(2)-1

    if (local_de_count==1) then
      call ESMF_GridGetCoord(dst_grid, coordDim=1, localDe=0, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lon_ptr, &
                             computationalLBound=lower, computationalUBound=upper, rc=rc)
      call check_esmf(rc, 'getting gaussian longitude storage')
      if (any(lower /= (/1,1/)) .or. any(upper /= dims)) call fail('Gaussian de does not cover global index space')
      call ESMF_GridGetCoord(dst_grid, coordDim=2, localDe=0, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lat_ptr, rc=rc)
      call check_esmf(rc, 'getting gaussian lattiude storage')

      do j=1,dims(2)
        do i=1,dims(1)
          lon_ptr(i,j) = real(gaussian_grid%lons(i), ESMF_KIND_R8)
          lat_ptr(i,j) = real(gaussian_grid%lats(j), ESMF_KIND_R8)
        end do
      end do
    
    end if ! (local_de_count==1)

    dst_field = ESMF_FieldCreate(dst_grid, typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 name='gaussian_scalar', rc=rc)
    call check_esmf(rc, 'creating gaussian field')
    call ESMF_FieldGet(dst_field, localDeCount=local_de_count, rc=rc)
    call check_esmf(rc, 'querying gaussian field local DEs')

    ! Construct the regridding route
    call ESMF_FieldRegridStore(srcField=src_field, dstField=dst_field, &
                               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
                               unmappedaction=ESMF_UNMAPPEDACTION_ERROR, &
                               routehandle=route, rc=rc)
    call check_esmf(rc, 'building the FV3-to-gaussian remapping')
    call ESMF_FieldRegrid(src_field, dst_field, routehandle=route, rc=rc)
    call check_esmf(rc, 'applying the fv3-to-gaussian remapping')

    ! Copy ESMF field data back into regular array
    allocate(gaussian_data(dims(1), dims(2)))
    if (local_de_count == 1) then
      call ESMF_FieldGet(dst_field, localDe=0, farrayPtr=dst_ptr, rc=rc)
      call check_esmf(rc, 'getting gaussian field storage')
      gaussian_data = dst_ptr
    endif

    ! Broadcast array back to other PETs
    call MPI_Bcast(gaussian_data, dims(1)*dims(2), MPI_REAL8, owner, comm, mpierr)
    call check_mpi(mpierr, 'broadcasting gaussian field')

    ! Finalize by releasing fields
    call ESMF_FieldRegridRelease(route, rc=rc)
    call check_esmf(rc, 'releasing regridding route')
    call ESMF_FieldDestroy(src_field, rc=rc)
    call check_esmf(rc, 'destroying fv3 field')
    call ESMF_FieldDestroy(dst_field, rc=rc)
    call check_esmf(rc, 'destroying gg field')
    call ESMF_GridDestroy(dst_grid, rc=rc)
    call check_esmf(rc, 'destroying gg grid')
  end subroutine remap_fv3_to_gaussian

  ! Error handling
  subroutine check_esmf(rc, operation)
    integer, intent(in) :: rc
    character(len=*), intent(in) :: operation
    if (rc /= ESMF_SUCCESS) then
      write(0, '(a, a, a, i0)') 'ESMF failure while ', operation, ': rc=', rc
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      error stop 'ESMF abort returned unexpectedly'
    end if
  end subroutine check_esmf

  subroutine check_mpi(mpierr, operation)
    integer, intent(in) :: mpierr
    character(len=*), intent(in) :: operation
    if (mpierr /= MPI_SUCCESS) then
      write(0, '(a, a, a, i0)') 'MPI failure while ', operation, ': error=', mpierr
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      error stop 'ESMF abort returned unexpectedly'
    end if
  end subroutine check_mpi

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(0, '(a, a)') 'FV3-to-Gaussian remap: ', message
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
    error stop 'ESMF abort returned unexpectedly'
  end subroutine fail

end module fv3_to_gaussian_remap_mod
