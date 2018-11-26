program test_read_db
!
! Reading the records stored inside the binary database files.
! It is a test run. Its goal is to check if the records are
! readable.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2018112600
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
integer(C_INT)            :: i
integer(C_INT)            :: max_pdb_num
!
real(C_FLOAT),allocatable :: cb(:,:)
real(C_FLOAT),allocatable :: ha(:,:)
!
character(len=4096)       :: src_folder
character(len=3)          :: resname(20)
character(len=4096)       :: file_name

src_folder='/home/vesso/project/db/'
!
resname=['ARG','ASN','ASP','CYI','CYS','GLN','GLU','HIE','HIP','HIS',&
         'ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
!
do i=1,ubound(resname,1)
   !
   file_name=trim(adjustl(src_folder))//'/'//resname(i)//'_ND.bin'
   !
   open(unit=666,file=file_name, form='unformatted',access='direct',&
        status='old',recl=C_INT)
   !
   read(666,rec=1) max_pdb_num
   !
   close(unit=666)
   !
   allocate(cb(3,max_pdb_num))
   allocate(ha(3,max_pdb_num))
   !
   print *,resname(i),max_pdb_num
   !
   file_name=trim(adjustl(src_folder))//'/'//resname(i)//'_CB.bin'
   !
   open(unit=666,file=file_name,form='unformatted',access='direct',&
        status='old',recl=3*(max_pdb_num)*C_FLOAT)
   !
   read(666,rec=1) cb
   !
   close(unit=666)
   !
   file_name=trim(adjustl(src_folder))//'/'//resname(i)//'_HA.bin'
   !
   open(unit=666,file=file_name,form='UNFORMATTED',access='direct',&
        status='old',recl=3*(max_pdb_num)*C_FLOAT)
   !
   read(666,rec=1) ha
   !
   close(unit=666)
   !
   deallocate(cb)
   !
   deallocate(ha)
   !
end do

end program test_read_db
