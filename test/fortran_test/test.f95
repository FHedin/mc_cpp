program fortran_interface

  use iso_c_binding
  implicit none
  
!   interface
!   
!     subroutine set_PBC(pbtype,a,b,c,alpha,beta,gamma) bind(c,name="set_pbc_")
!       use iso_c_binding
!       implicit none
!       character(kind=c_char) :: pbtype(128)
!       real(kind=c_double) :: a,b,c,alpha,beta,gamma
!     end subroutine
!     
!     subroutine set_ENS(natom,temperature) bind(c,name="set_ens_")
!       use iso_c_binding
!       implicit none
!       integer(kind=c_int) :: natom
!       real(kind=c_double) :: temperature
!     end subroutine
!     
!     subroutine set_FF_and_COR(cutoffmode,ctoff,cton,dcut,FFfile,CORfile) bind(c,name="set_ff_and_cor_")
!       use iso_c_binding
!       implicit none
!       character(kind=c_char) :: cutoffmode(128)
!       real(kind=c_double) :: ctoff, cton, dcut
!       character(kind=c_char) :: FFfile(128), CORfile(128)
!     end subroutine
!     
!     subroutine clean() bind(c,name="clean_")
!       use iso_c_binding
!       implicit none
!     end subroutine
!   
!   end interface
  
  character(len=128) :: pbtype
  real(kind=c_double) :: a,b,c,alpha,beta,gamma
  
  integer(kind=c_int) :: natom
  real(kind=c_double) :: temperature
  
  character(len=128) :: cutoffmode
  real(kind=c_double) :: ctoff, cton, dcut
  
  character(len=128) :: FFfile, CORfile
  
  real(kind=c_double),dimension(10) :: ener

  !------------------
  
  ! First define periodic boundary conditions
  ! use 'none' or 'cubic' or 'orbic' (orthorombic)
  pbtype='cubic'//C_NULL_CHAR
  a=30.0
  b=30.0
  c=30.0
  alpha=90.0
  beta=90.0
  gamma=90.0
  
  call set_pbc(pbtype,a,b,c,alpha,beta,gamma)
  
  ! Then define the ensemble
  natom=2658
  temperature=300.0 ! in K
  
  call set_ens(natom,temperature)
  
  ! set the forcefield and load parameters + coordinates files
  ! cuttof mode is 'full' (full non bonded energy evaluation) or 'switch'
  cutoffmode='full'//C_NULL_CHAR
  ctoff=12.0
  cton=10.0
  dcut=2.0
  
  ! FFfile is in the MDBAS format (see provided conversion tool)
  FFfile='./nma_solvated/forfield.dat'//C_NULL_CHAR
  ! a standard charmm cor file
  CORfile='./nma_solvated/nma.cor'//C_NULL_CHAR
  
  call set_ff_and_cor(cutoffmode,ctoff,cton,dcut,FFfile,CORfile)
  
  ! get the energy
  call get_energy(ener)
  write(6,*) "TOTAL POTENTIAL KINETIC ELEC VDW"
  write(6,*) ener(1:5)
  write(6,*) "BOND ANGLES UREY-BRADLEY DIHE IMPR"
  write(6,*) ener(6:10)
  
  
  ! clean things (memory deallocation of c++)
  call clean()
  
  return
  
end program