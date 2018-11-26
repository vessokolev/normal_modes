program create_cb_ha_db

use iso_c_binding,only:C_INT,C_FLOAT
use mod_pdb
use mod_process

implicit none

integer(C_INT)            :: i
integer(C_INT)            :: j
integer(C_INT)            :: num_atoms
integer(C_INT)            :: min_pdb_num
integer(C_INT)            :: max_pdb_num

real(C_FLOAT)             :: angle_2d
real(C_FLOAT)             :: angle_3d
real(C_FLOAT)             :: ca(3)
real(C_FLOAT)             :: axis(3)
real(C_FLOAT)             :: coord(3,6)
real(C_FLOAT),allocatable :: cb(:,:)
real(C_FLOAT),allocatable :: ha(:,:)

logical                   :: flag

character(len=4096)       :: src_folder
character(len=4096)       :: dst_folder
character(len=10)         :: str_tmp
character(len=3)          :: resname
character(len=1024)       :: pdb_prefix
character(len=4)          :: pdb_ext
character(len=4096)       :: file_name

if (command_argument_count()>=3) then
   !
   call get_command_argument(1,src_folder)
   call get_command_argument(2,resname)
   call get_command_argument(3,dst_folder)
   !
else
   !
   print *
   print *,'Invoke the program with the following arguments:'
   print *
   print *,'src_folder resname min_pdb_num max_pdb_num dst_folder'
   print *
   !
   call exit(1)
   !
end if
!
pdb_prefix='cg_'
pdb_ext='.pdb'
!
min_pdb_num=0
!
max_pdb_num=0
!
do
   !
   write(file_name,fmt='(I8.8)') max_pdb_num
   !
   file_name=trim(adjustl(src_folder))//'/'//&
             trim(adjustl(resname))//'/cg/'//&
             trim(adjustl(pdb_prefix))//&
             trim(adjustl(file_name))//&
             trim(adjustl(pdb_ext))
   !
   inquire(file=file_name,exist=flag)
   !
   if (.not.flag) then
      !
      exit
      !
   end if
   !
   max_pdb_num=max_pdb_num+1
   !
end do
!
print *
print *,'Found ',max_pdb_num,'files!'
print *
!
max_pdb_num=max_pdb_num-1
!
allocate(cb(3,max_pdb_num+1))
allocate(ha(3,max_pdb_num+1))
!
! Obtain the number of atoms per PDB, before running
! the loop on the PDB for analyzing the coordinates.
!
!
! Construct the file name based on the number of frame with
! leading zeroes:
!
write(file_name,fmt='(I8.8)') min_pdb_num
!
file_name=trim(adjustl(src_folder))//'/'//&
          trim(adjustl(resname))//'/cg/'//&
          trim(adjustl(pdb_prefix))//&
          trim(adjustl(file_name))//&
          trim(adjustl(pdb_ext))
!
! Count the number of atom records inside the PDB file:
!
call count_atom_entries(file_name,num_atoms)
!
j=0
!
do i=min_pdb_num,max_pdb_num
   !
   write(file_name,fmt='(I8.8)') i
   !
   file_name=trim(adjustl(src_folder))//'/'//&
             trim(adjustl(resname))//'/cg/'//&
             trim(adjustl(pdb_prefix))//&
             trim(adjustl(file_name))//&
             trim(adjustl(pdb_ext))
   !
   ! Reads the atom coordinates from the PDB file:
   !
   call read_cg_atom_entries(file_name,coord)
   !
   ! Canonize the atom coordinates of the residue. During that process
   ! the "N", "C", and "CA" atoms are placed on the plane perpendicular
   ! to the z-axis, and the coordinates are shifted so that those of "CA"
   ! to become (0,0,0).
   !
   call canonize_residue_coords(coord,ca,angle_2d,angle_3d,axis)
   !
   ! Obtain the records for the database:
   !
   j=j+1
   !
   call get_angles_bond_len_pres(coord,cb(:,j),ha(:,j))
   !
end do
!
print *
!
file_name=trim(adjustl(dst_folder))//'/'//trim(adjustl(resname))//'_CB.bin'
!
open(unit=666,file=file_name,form='UNFORMATTED',access='direct',&
     recl=3*j*C_FLOAT)
!
write(666,rec=1) cb
!
close(unit=666)
!
print *,'Created the binary file ',trim(adjustl(file_name))
!
file_name=trim(adjustl(dst_folder))//'/'//trim(adjustl(resname))//'_HA.bin'
!
open(unit=666,file=file_name,form='UNFORMATTED',access='direct',&
     recl=3*j*C_FLOAT)
!
write(666,rec=1) ha
!
close(unit=666)
!
print *,'Created the binary file ',trim(adjustl(file_name))
!
file_name=trim(adjustl(dst_folder))//'/'//trim(adjustl(resname))//'_ND.bin'
!
open(unit=666,file=file_name, form='UNFORMATTED',access='direct',&
     recl=C_INT)
!
write(666,rec=1) j
!
close(unit=666)
!
print *,'Created the binary file ',trim(adjustl(file_name))
!
print *

end program create_cb_ha_db
