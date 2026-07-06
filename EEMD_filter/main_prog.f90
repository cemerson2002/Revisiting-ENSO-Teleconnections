program read_era5_netcdf
  use netcdf_reader
  use netcdf_writer_1x1
  use spatiotemporal_filter
  implicit none
 
  character(len=256) :: in_dir, in_fname, ivarname, tname, levname
  character(len=256) :: out_dir, out_fname, ovarname, var_lname, var_units
  character(len=256) :: time_units, time_calendar
  character(len=256) :: in_file, out_file
  real :: level_val, fill_value

  namelist /netcdf_infile/ in_dir, in_fname, ivarname, tname, levname, level_val
  namelist /netcdf_outfile/ out_dir, out_fname, ovarname, var_lname, var_units, time_units, time_calendar, fill_value

  real, allocatable :: idata(:,:,:)
  real, allocatable :: iidata(:,:,:)
  real, allocatable :: odata(:,:,:,:)
  real, allocatable :: latitude(:), longitude(:)
  integer, allocatable :: time(:), modes(:), iskip
  integer :: i, nlat, nlon, ntime, nmodes, i1, j1, nx, ny, j, k
  real :: min_val, max_val

  ! Read namelist ! 
  open(unit=10, file='input.nml', status='old', action='read')
  read(10, nml=netcdf_infile)
  rewind(10)
  read(10, nml=netcdf_outfile)
  close(10)

  in_file = trim(in_dir) // trim(in_fname)
  out_file = trim(out_dir) // trim(out_fname)

  ! Read NetCDF data
  call read_netcdf(in_file, ivarname, tname, levname, level_val, idata, latitude, longitude, nlat, nlon, ntime)
  print*, 'ntime, nlat, nlon ', ntime, nlat, nlon
  nmodes = nint(log(float(ntime))/log(2.0))+2
  allocate(modes(nmodes))
  modes(1) = 1

!  ntime = 365
  allocate(time(ntime))
  do i = 1, ntime
    time(i) = i
  end do
! coarsening the input data
  iskip = 2 ! changed 4 to 1
  nx = nlon/iskip
  ny = nlat/iskip + 1
!
  print*, 'nx = ', nlon, nx
  print*, 'ny = ', nlat, ny
  print*, 'nmodes = ', nmodes
!
  allocate(iidata(nx,ny,ntime))
!
  do k=1,ntime 
   i1 = 0
   do i = 1, nlon, iskip 
    i1 = i1 +1
    j1 = 0
    do j= 1, nlat, iskip 
     j1 = j1 + 1 
     iidata(i1,j1,k) = idata(i,j,k)
    end do
   end do
  end do

  allocate(odata(nx,ny,ntime,nmodes))
 
  call filter(iidata,odata,nx,ny,nx*ny,ntime,nmodes)

! Find min/max values
  print*,'after filter '
  min_val = minval(odata)
  max_val = maxval(odata)

  print *, "Min " // trim(ivarname) // " height:", min_val
  print *, "Max " // trim(ivarname) // " height:", max_val

! Write subset to new NetCDF file
  call write_netcdf(out_file, ovarname, var_lname, var_units, time_units, time_calendar, &
          fill_value, odata, time, modes, ntime, nmodes)

  ! Deallocate memory
  deallocate(idata, odata, latitude, longitude)

end program read_era5_netcdf
