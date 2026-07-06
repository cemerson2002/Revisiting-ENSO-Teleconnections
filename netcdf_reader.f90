module netcdf_reader
  ! This module is more generic to handle 3D and 4D data !
  ! Check the dimension order of the variable in the NetCDF file !
  ! Read short, float, or double NetCDF variable types, converting them to real(kind=4) !
  ! Allow selecting a pressure level (e.g., 500 hPa) through a namelist file !
  ! Replace the fill/missing values to -9999.0 !
  ! This is compatible with NaN values !
  use netcdf
  implicit none
  private
  public :: read_netcdf
  contains
  
  subroutine read_netcdf(filename, varname, tname, levname, level_val, var, latitude, longitude, nlat, nlon, ntime)
    character(len=*), intent(in) :: filename, varname
    character(len=*), intent(in) :: tname
    character(len=*), intent(in), optional :: levname
    character(len=NF90_MAX_NAME) :: dimnames(4)
    character(len=NF90_MAX_NAME) :: dimname
    real, intent(in), optional :: level_val
    real, allocatable, intent(out) :: var(:,:,:)
    real, allocatable, intent(out) :: latitude(:), longitude(:)
    integer, intent(out) :: nlat, nlon

    integer :: i, ndims, nvars, natts, unlimdimid, dimlen, nlev, ntime, ilev
    integer :: ncid, varid, varid_lat, varid_lon, varid_lev, status
    integer :: dimid_time, dimid_lat, dimid_lon, dimid_lev, dimids(NF90_MAX_VAR_DIMS)
    real, allocatable :: level(:)
    integer :: start3(3), count3(3)
    integer :: start4(4), count4(4)
    integer :: vartype
    character(len=32) :: xtype_name
    real :: scale_factor, add_offset
    logical :: has_scale, has_offset
    logical :: has_level
    real, parameter :: nfill_value = -9999.0 
    real :: fill_val, min_val, max_val
    logical :: has_fillval, has_nans

    ! Temporary variables for holding raw data
    integer*2, allocatable :: tmp_short3(:,:,:),tmp_short4(:,:,:,:)
    real, allocatable :: tmp_real3(:,:,:),tmp_real4(:,:,:,:)
    double precision, allocatable :: tmp_double3(:,:,:),tmp_double4(:,:,:,:)

    ! Open the NetCDF file
    status = nf90_open(filename, NF90_NOWRITE, ncid)
    status = nf90_inquire(ncid, ndims, nvars, natts, unlimdimid)

    print *, "======== Dimensions ========"
    has_level = .false.
    do i = 1, ndims
        status = nf90_inquire_dimension(ncid, i, dimname, dimlen)
        print *, "Dimension: ", trim(dimname), " Size: ", dimlen
        if (trim(dimname) .eq. trim(levname)) then
                has_level = .true.
        end if
    end do
    print *, "============================"

    ! Check Variable
    status = nf90_inq_varid(ncid, trim(varname), varid)
    status = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids, xtype=vartype)

    ! Get and print dimension names
    do i = 1, ndims
       status = nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i))
       if (status /= nf90_noerr) stop "Error inquiring dimension"
       print *, "Dimension ", i, ": ", trim(dimnames(i))
    end do

    ! Check if it's either lon-lat-time or lon-lat-lev-time
    print *, ndims
    if (ndims == 3) then
       if (.not. ((trim(dimnames(1)) == "lon" .or. trim(dimnames(1)) == "longitude") .and. &
          (trim(dimnames(2)) == "lat" .or. trim(dimnames(2)) == "latitude") .and. &
          (trim(dimnames(3)) == trim(tname)))) then
          print *, "This module won't work: expected dimension order [time,lat,lon]"
          stop
       else
          print *, "Netcdf variable dimension order : [time,lat,lon]"
       end if
    elseif (ndims == 4) then
       if (.not. ((trim(dimnames(1)) == "lon" .or. trim(dimnames(1)) == "longitude") .and. &
          (trim(dimnames(2)) == "lat" .or. trim(dimnames(2)) == "latitude") .and. &
          (trim(dimnames(3)) == trim(levname) ) .and. (trim(dimnames(4)) == trim(tname)))) then
          print *, "This module won't work: expected dimension order [time,lev,lat,lon]"
          stop
       else
          print *, "Netcdf variable dimension order : [time,lev,lat,lon]"
       end if
    else
       print *, "This module won't work: unsupported number of dimensions"
       stop
    end if
    print *, "============================"

    ! Get dimension IDs and sizes
    status = nf90_inq_dimid(ncid, trim(tname), dimid_time)
    status = nf90_inquire_dimension(ncid, dimid_time, len=ntime)

    status = nf90_inq_dimid(ncid, "latitude", dimid_lat)
    if (status /= nf90_noerr) then
        status = nf90_inq_dimid(ncid, "lat", dimid_lat)
    end if
    status = nf90_inquire_dimension(ncid, dimid_lat, len=nlat)

    status = nf90_inq_dimid(ncid, "longitude", dimid_lon)
    if (status /= nf90_noerr) then
        status = nf90_inq_dimid(ncid, "lon", dimid_lon)
    end if
    status = nf90_inquire_dimension(ncid, dimid_lon, len=nlon)

    if (has_level) then
            status = nf90_inq_dimid(ncid, trim(levname), dimid_lev)
            status = nf90_inquire_dimension(ncid, dimid_lev, len=nlev)
    end if

    ! Allocate memory
    allocate(latitude(nlat))
    allocate(longitude(nlon))
    allocate(level(nlev))
    allocate(var(nlon,nlat,ntime))

    ! Get Lats and Lons
    status = nf90_inq_varid(ncid, 'latitude', varid_lat)
    if (status /= nf90_noerr) then
        status = nf90_inq_varid(ncid, "lat", varid_lat)
    end if
    status = nf90_get_var(ncid, varid_lat, latitude)

    status = nf90_inq_varid(ncid, 'longitude', varid_lon)
    if (status /= nf90_noerr) then
        status = nf90_inq_varid(ncid, "lon", varid_lon)
    end if
    status = nf90_get_var(ncid, varid_lon, longitude)
    
    if (has_level) then
       status = nf90_inq_varid(ncid, trim(levname), varid_lev)
       status = nf90_get_var(ncid, varid_lev, level)
       print *, level
       ilev = 1
       if (present(level_val)) then
          ilev = minloc([ (abs(level(i) - level_val), i=1, nlev) ], 1)
       end if
    print *, "Selecting level : ", level_val
    end if

    ! Check for _FillValue or missing_value
    has_fillval = .false.
    has_nans = .false.
    status = nf90_get_att(ncid, varid, "_FillValue", fill_val)
    if (status == nf90_noerr) then
       has_fillval = .true.
       if (fill_val - (-1*fill_val) == fill_val) has_nans = .true.
    else
       status = nf90_get_att(ncid, varid, "missing_value", fill_val)
       if (status == nf90_noerr) has_fillval = .true.
       if (fill_val - (-1*fill_val) == fill_val) has_nans = .true.
    end if

    print *, "Actual fill/missing value from file: ", fill_val
    print *, "Replacing fill/missing value : ", nfill_value
    print *, "============================"

    scale_factor = 1.0
    add_offset = 0.0
    has_scale = .false.
    has_offset = .false.

    ! Check for attributes
    if (nf90_get_att(ncid, varid, "scale_factor", scale_factor) == nf90_noerr) has_scale = .true.
    if (nf90_get_att(ncid, varid, "add_offset", add_offset) == nf90_noerr) has_offset = .true.

    if (.not. has_level) then
       start3 = (/1, 1, 1/)
       count3 = (/nlon, nlat, ntime/)

       select case (vartype)
       case (nf90_short)
          allocate(tmp_short3(nlon, nlat, ntime))
          status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
          status = nf90_get_att(ncid, varid, "add_offset", add_offset)
          status = nf90_get_var(ncid, varid, tmp_short3, start=start3, count=count3)
          var = scale_factor * real(tmp_short3) + add_offset
          deallocate(tmp_short3)

       case (nf90_real)
          allocate(tmp_real3(nlon, nlat, ntime))
          status = nf90_get_var(ncid, varid, tmp_real3, start=start3, count=count3)
          var = tmp_real3
          deallocate(tmp_real3)

       case (nf90_double)
          allocate(tmp_double3(nlon, nlat, ntime))
          status = nf90_get_var(ncid, varid, tmp_double3, start=start3, count=count3)
          var = real(tmp_double3)
          deallocate(tmp_double3)

       case default
          print *, "Unsupported variable type."
          stop
       end select
    else
       start4 = (/1, 1, ilev, 1/)
       count4 = (/nlon, nlat, 1, ntime/)
       select case (vartype)
       case (nf90_short)
          allocate(tmp_short4(nlon, nlat, 1, ntime))
          status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
          status = nf90_get_att(ncid, varid, "add_offset", add_offset)
          status = nf90_get_var(ncid, varid, tmp_short4, start=start4, count=count4)
          var = scale_factor * real(tmp_short4(:,:,1,:)) + add_offset
          deallocate(tmp_short4)

       case (nf90_real)
          allocate(tmp_real4(nlon, nlat, 1, ntime))
          status = nf90_get_var(ncid, varid, tmp_real4, start=start4, count=count4)
          var = tmp_real4(:,:,1,:)
          deallocate(tmp_real4)

       case (nf90_double)
          allocate(tmp_double4(nlon, nlat, 1, ntime))
          status = nf90_get_var(ncid, varid, tmp_double4, start=start4, count=count4)
          var = real(tmp_double4(:,:,1,:))
          deallocate(tmp_double4)

       case default
          print *, "Unsupported variable type."
          stop
       end select
    end if
   
    call fill_missing_values(var, fill_val, nfill_value, has_nans)

    min_val = minval(var)
    max_val = maxval(var)

    print *, "Min " // trim(varname) // " :", min_val
    print *, "Max " // trim(varname) // " :", max_val
    print *, "======================"

    if (status /= nf90_noerr) then
        print *, "Error opening file: ", trim(nf90_strerror(status))
        stop
    end if

    ! Close the file
    status = nf90_close(ncid)
  end subroutine read_netcdf

  subroutine fill_missing_values(var, fill_val, nfill_value, has_nans)
    implicit none

    real, intent(inout) :: var(:,:,:)
    real, intent(in)    :: fill_val
    real, intent(in)    :: nfill_value
    logical, intent(in) :: has_nans

    integer :: i, j, k
    logical :: is_nan

    ! Loop through all elements in the 3D array
    do i = 1, size(var, 1)
       do j = 1, size(var, 2)
          do k = 1, size(var, 3)
             ! Check if the current value is NaN or equal to the fill value
             if (has_nans) then
                is_nan = (var(i,j,k) + 1.0 == var(i,j,k))
                if (is_nan) then
                   var(i,j,k) = nfill_value 
                end if
             else
                if (var(i,j,k) == fill_val) then
                   var(i,j,k) = nfill_value 
                end if
             end if
          end do
       end do
    end do
    end subroutine fill_missing_values
end module netcdf_reader
