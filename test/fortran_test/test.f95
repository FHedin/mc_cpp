program fortran_interface

  ! this is required for compatibility with C/C++
  use iso_c_binding
  
  implicit none

  character(len=128) :: pbtype
  real(kind=c_double) :: a,b,c,alpha,beta,gamma
  
  integer(kind=c_int) :: natom
  real(kind=c_double) :: temperature
  
  character(len=128) :: cutoffmode
  real(kind=c_double) :: ctoff, cton, dcut
  
  character(len=128) :: FFfile, CORfile
  
  real(kind=c_double),dimension(10) :: ener
  
  real(kind=c_double),allocatable,dimension(:) :: x,y,z,mass
  
  character(len=128) :: move_type, move_mode, sel_mode, selection
  integer(kind=c_int) :: nstep
  
  real(kind=c_double) :: maxTrn, maxRot
  
  !------------------
  
  !uncomment this if you want to remove standard text output from the c++ : it will be silent, i.e. only warnings and fatal errors will cause text to appear to the terminal
  call disable_text_output()
  ! the following restores normal test output to console
  ! call restore_text_output()
  
  ! First define periodic boundary conditions
  ! use 'none' or 'cubic' or 'orbic' (orthorombic)
  ! for passing a character string from fortran to C/C++ this '//C_NULL_CHAR' is required at the end : it concatenates the symbol '\0' to the end of the string
  pbtype='none'//C_NULL_CHAR
  
  ! if no periodic conditions you can remove the next 6 lines or just keep everything to 0
  !a=65.1907
  !b=58.9821
  !c=46.5648
  !alpha=90.0
  !beta=90.0
  !gamma=90.0
  
  call set_pbc(pbtype,a,b,c,alpha,beta,gamma)
  
  ! Then define the ensemble
  ! This 'ensemble' internal representation was useful for my Monte-Carlo code but is not really useful when just energy evaluation is desired
  ! the natom is required in the code so always set it
  ! if this natom is not the same as the number of atoms in the coordinate file the code will fail later
  natom=3
  
  ! temperature is not used for the moment so set it or not
  temperature=300.0
  
  call set_ens(natom,temperature)
  
  ! set the forcefield and load parameters + coordinates files
  ! cuttof mode is 'full' (full non-bonded energy evaluation) or 'switch' (using switching function for non-bonded interactions)
  ! if not 'full' then a verlet list will be constructed
  ! the exclusion list for LJ 1,4 interactions is always constructed
  cutoffmode='full'//C_NULL_CHAR
  ! if switch mode enable we need a cutoff and cuton 
  !ctoff=12.0 ! default 12.0 in most of simulations
  !cton=10.0  ! default 10.0 in most of simulations
  ! dcut is an extra distance used for building the verlet list, usually 2, it is added to cutoff so the verlet list is built up 
  ! to a distance of (cutoff+dcut) ; thus the list contains more atoms than required but it will avoid updating the list to regularly
  !dcut=2.0
  
  ! FFfile is in the MDBAS format (see provided conversion tool)
  FFfile='clcn/ffield.dat'//C_NULL_CHAR
  ! a standard charmm cor file
  CORfile='clcn/clcn-opt.cor'//C_NULL_CHAR
  
  call set_ff_and_cor(cutoffmode,ctoff,cton,dcut,FFfile,CORfile)
  
  ! get the initial energy of the system
  ! an array of 10 double precision real numbers
  !   ener(1) = total energy
  !   ener(2) = potential energy
  !   ener(3) = kinetic energy
  !   ener(4) = electrostatic energy
  !   ener(5) = vdw energy
  !   ener(6) = bond energy
  !   ener(7) = angles energy
  !   ener(8) = urey bradley energy
  !   ener(9) = dihedrals energy
  !   ener(10) = impropers energy
  call get_energy(ener)
  
  write(6,*) "TOTAL POTENTIAL KINETIC ELEC VDW"
  write(6,*) ener(1:5)
  write(6,*) "BOND ANGLES UREY-BRADLEY DIHE IMPR"
  write(6,*) ener(6:10)

!   write(6,*) "TOTAL E : ", ener(1)

!   allocate(x(natom))
!   allocate(y(natom))
!   allocate(z(natom))
! updating the coordinates ; please be sure that the size of the arrays is really natom otherwise it may crash or produce weird things
!   x=1.d0
!   y=2.d0
!   z=3.d0
!   call set_coords(x,y,z)
  
  ! check the energy again
!   call get_energy(ener)
!   write(6,*) "TOTAL POTENTIAL KINETIC ELEC VDW"
!   write(6,*) ener(1:5)
!   write(6,*) "BOND ANGLES UREY-BRADLEY DIHE IMPR"
!   write(6,*) ener(6:10)
  
  ! get the atomic masses
  ! same warning as before i.e. allocation size should be the natom provided to set_ens
  !allocate(mass(natom))
  !call get_mass(mass)
  !write(6,*) "Dump of atomic mass for each atom"
  !write(6,*) mass
  
  ! define allowed types of random moves
  !
  !   * move_type is either "rot" for rotations or "trn" for translation
  !   * move_mode is "all", "residue" or "atom" : when applying a move it decides if we apply this to all
  !     atoms of the selection at once, or only to the residue containing the atom selection, or only to the atom selection.
  ! 
  !   * note that move_type="rot" and  move_mode="atom" is impossible as we can't rotate just one atom
  ! 
  !   * Then the atom selection is given using this selection node
  !   * sel_mode are :
  !     "all" or "none" which are explicit enough...
  !     "residue_id" which uses the column 2 of cor files, "residue_name" using column 3
  !     "atom_idx" which uses the column 1 of the cor file
  ! 
  !   * sele="..." is simply the string containing the selection (unused for sel_mode "all" or "none")
  !
  
  ! translation allowed for all atoms, and at once, by the same random distance : won't change energy
!   move_type="trn"//C_NULL_CHAR
!   move_mode="all"//C_NULL_CHAR
!   sel_mode="all"//C_NULL_CHAR
!   selection = ""//C_NULL_CHAR
!   call add_move_type(move_type,move_mode,sel_mode,selection)
  
  ! translation allowed for all atoms, by atom, so only one atom translated for a given move attempt
!   move_type = "trn"//C_NULL_CHAR
!   move_mode = "atom"//C_NULL_CHAR
!   sel_mode  = "all"//C_NULL_CHAR
!   selection = ""//C_NULL_CHAR
!   call add_move_type(move_type,move_mode,sel_mode,selection)
    
  ! rotation allowed for all atoms, and at once, by the same random angle : won't change energy
!   move_type="rot"//C_NULL_CHAR
!   move_mode="all"//C_NULL_CHAR
!   sel_mode="all"//C_NULL_CHAR
!   selection = ""//C_NULL_CHAR
!   call add_move_type(move_type,move_mode,sel_mode,selection)

  ! Per residue based rotation of atoms defined in the selection (all)
!   move_type="rot"//C_NULL_CHAR
!   move_mode="residue"//C_NULL_CHAR
!   sel_mode="all"//C_NULL_CHAR
!   selection = ""//C_NULL_CHAR
!   call add_move_type(move_type,move_mode,sel_mode,selection)

  ! print initial coordinates to file
  ! if file exist coordinates appended to the end
  call write_xyz("init.xyz"//C_NULL_CHAR)
  
  !allow max translation of +/- 1 angstroem
!   maxTrn = 0.25
  !allow max angular rotation of +/- 15 degrees
!   maxRot = 15.0

  ! apply random move to the first registered move : energy should be the same but coordinates will change

!   maxTrn = 1.0
!   call random_move(1,maxTrn)
!   !get modified coordinates
!   call write_XYZ("move1.xyz"//C_NULL_CHAR)
!   call get_energy(ener)
!   write(6,*) "TOTAL E : ", ener(1)

  
!   ! apply random move to the second registered move several time : energy should change
!   do nstep=1,1000
!     call random_move(2,maxTrn)
!     call write_XYZ("move2.xyz"//C_NULL_CHAR)
!   enddo
!   call get_energy(ener)
!   write(6,*) "TOTAL E : ", ener(1)
!   
!   ! apply random move to the third registered move : energy should not change
!   call random_move(3,maxRot)
!   call write_XYZ("move3.xyz"//C_NULL_CHAR)
!   call get_energy(ener)
!   write(6,*) "TOTAL E : ", ener(1)
!   
!     ! apply random move to the fourth registered move several time : energy should change
!   do nstep=1,1000
!     call random_move(4,maxRot)
!     call write_XYZ("move4.xyz"//C_NULL_CHAR)
!   enddo
!   call get_energy(ener)
!   write(6,*) "TOTAL E : ", ener(1)
!   
!   ! if we want to store coordinates back from the c++ code to the fortran arrays
!   call get_coords(x,y,z)

  ! clean things (memory deallocation of c++) do this we you don't need anymore any of the 
  ! subroutine used before
  call clean()

  return

end program
