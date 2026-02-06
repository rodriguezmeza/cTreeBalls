
  program dat2fits
  use fitstools, only : output_map
  use head_fits 
  ! CFITS & HEALPIX libraries are necessary to compile this code

  implicit none
  integer(4) :: nside
  integer(8) :: npix
  real(4), allocatable :: wl_data(:,:)
  character(len=200) :: infile,outfile
  character(len=80), dimension(1:64) :: mapheader

  infile='***.mag.dat'     ! input  mag.dat file
  outfile='***.fits'       ! output fits file
 
  open(12,file=infile,status="old",form='unformatted')
  read(12) nside,npix

  allocate(wl_data(0:npix-1,1:3))

  read(12) wl_data(0:npix-1,1)  ! convergence
  read(12) wl_data(0:npix-1,2)  ! shear1
  read(12) wl_data(0:npix-1,3)  ! shear2
!  read(12) wl_data(0:npix-1,4)  ! rotation
  close(12)

  call system('rm  '//outfile)

!!!! output fits 
  call write_minimal_header(mapheader,'MAP',nside=nside,ordering='RING')
  call add_card(mapheader,'TTYPE1','convergence','label for field 1')
  call add_card(mapheader,'TTYPE2','shear1','label for field 2')
  call add_card(mapheader,'TTYPE3','shear2','label for field 3')
!  call add_card(mapheader,'TTYPE4','rotation','label for field 4')

  call output_map(wl_data,mapheader,outfile)

  deallocate(wl_data)

  end program

       


