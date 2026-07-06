module netcdf_writer_1x1
  ! This module is tailor made for writing 4D geopotential !
  use netcdf
  implicit none
  contains

  subroutine write_netcdf(output_filename, ovarname, vlname, vunits, tunits, tcalendar, &
             fill_value, var, time, modes, ntime, nmodes)
    character(len=*), intent(in) :: output_filename
    character(len=*), intent(in) :: vlname, vunits, tunits, tcalendar, ovarname
    real, intent(in) :: var(:, :, :, :)   ! (nlon, nlat,  ntime,nmodes)
    real, intent(in) :: fill_value
    integer, intent(in) :: time(ntime), modes(nmodes)
    integer, intent(in) :: ntime, nmodes

    ! Hardcoded dimensions
    integer, parameter :: nlat = 181, nlon = 360 ! changed from 360 to 288
    real :: latitude(nlat), longitude(nlon)

    integer :: i, new_ncid, status
    integer :: dimids(4), varid, varid_lat, varid_lon, varid_time, varid_modes

    ! Generate latitude and longitude
    do i = 1, nlat
       latitude(i) = -90.0 + real(i-1) * 1.0
    end do
    do i = 1, nlon
       longitude(i) = 0.0 + real(i-1) * 1.0 ! changed from 1.0 to 1.25
    end do

    print *, latitude
    print *, longitude

    ! Create new NetCDF file
    status = nf90_create(output_filename, nf90_clobber, new_ncid)

    ! Define dimensions
    status = nf90_def_dim(new_ncid, 'longitude', nlon, dimids(1))
    status = nf90_def_dim(new_ncid, 'latitude', nlat, dimids(2))
    status = nf90_def_dim(new_ncid, 'time', ntime, dimids(3))
    status = nf90_def_dim(new_ncid, 'modes', nmodes, dimids(4))

    ! Define variables
    status = nf90_def_var(new_ncid, 'longitude', nf90_real, dimids(1), varid_lon)
    status = nf90_def_var(new_ncid, 'latitude', nf90_real, dimids(2), varid_lat)
    status = nf90_def_var(new_ncid, 'time', nf90_int, dimids(3), varid_time)
    status = nf90_def_var(new_ncid, 'modes', nf90_int, dimids(4), varid_modes)
    status = nf90_def_var(new_ncid, trim(ovarname), nf90_real, dimids, varid)

    status = nf90_put_att(new_ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(new_ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(new_ncid, varid_modes, 'units', 'unitless')
    status = nf90_put_att(new_ncid, varid_time, 'units', trim(tunits))
    status = nf90_put_att(new_ncid, varid_time, 'calendar', trim(tcalendar))
    status = nf90_put_att(new_ncid, varid, 'units', trim(vunits))
    status = nf90_put_att(new_ncid, varid, 'long_name', trim(vlname))
    !status = nf90_put_att(new_ncid, varid, 'level', "500")
    status = nf90_put_att(new_ncid, varid, '_FillValue', fill_value)

    ! End define mode
    status = nf90_enddef(new_ncid)

    ! Write data
    status = nf90_put_var(new_ncid, varid_lat, latitude)
    status = nf90_put_var(new_ncid, varid_lon, longitude)
    status = nf90_put_var(new_ncid, varid_time, time)
    status = nf90_put_var(new_ncid, varid_modes, modes)
    status = nf90_put_var(new_ncid, varid, var)
    if (status /= nf90_noerr) then
        print *, "Error writing file: ", trim(nf90_strerror(status))
        stop
    end if

    ! Close the file
    status = nf90_close(new_ncid)
  end subroutine write_netcdf

end module netcdf_writer_1x1
