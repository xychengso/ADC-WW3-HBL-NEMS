------------
mpifort -c  -O3 -assume byterecl -g -traceback -ip -module /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/WW3/model/mod -convert big_endian -I/home/xychen/ADC-WW3-Coupling-Dev/SrcCode/WW3/model/ftn -I/projects/ees/dhs-crc/xychen/ESMFlib_openmpi_xyc/mod/modg/Linux.intel.64.openmpi.default -I/projects/ees/dhs-crc/xychen/ESMFlib_openmpi_xyc/include -I/usr/share/Modules/software/CentOS-7/netcdf-Fortran/4.4.0_intel-18.0.0/include 
------------
------------
w3arrymd.F90(2091): remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'.
               '  Spectrum : ',A,'  Units : ',E8.3,1X,A,        &
-----------------------------------------------^
w3arrymd.F90(2092): remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'.
               '  Maximum value : ',E8.3,1X,A/)
-------------------------------------^
w3arrymd.F90(2089): remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'.
               '  Maximum value : ',E8.3,1X,A/)
-------------------------------------^
------------