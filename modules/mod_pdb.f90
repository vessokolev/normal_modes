module mod_pdb

   ! The module contains the subroutines for reading the PDB files
   ! containing MOLARIS type CG protein residue topology.
   !
   ! Author : Veselin Kolev <vesso.kolev@gmail.com>
   ! Version: 2018112600
   ! License: GPLv2

   use iso_c_binding,only: C_INT,C_FLOAT
   !
   implicit none


contains


    subroutine count_atom_entries(pdb_in,num_atoms)
    !
    ! Counts the number ot 'ATOM' or/and 'HETATM' records
    ! stored inside the examined PDB file.
    !
    ! Interface variables:
    !
    character(len=*),intent(in) :: pdb_in
    integer(C_INT),intent(out)  :: num_atoms
    !
    ! Local variables:
    !
    character(len=54)           :: line
    integer(C_INT)              :: resnum
    integer(C_INT)              :: stat

    num_atoms=0
    !
    open(666,file=pdb_in,status='old')
    !
    do
      !
      read(666,fmt='(A54)',iostat=stat) line
      !
      if (stat/=0) then
         !
         exit
         !
      else
         !
         if (line(1:6)=='ATOM  '.or.line(1:6)=='HETATM') then
            !
            num_atoms=num_atoms+1
            !
         end if
         !
      end if
      !
    end do
    !
    close(666)

    end subroutine count_atom_entries


    subroutine read_cg_atom_entries(pdb_in,coord)
    !
    ! Reads the coordinates of 'N, 'CA', 'HA', 'C', and 'CB'
    ! atoms from within the PDB file.
    !
    character(len=*),intent(in) :: pdb_in
    real(C_FLOAT),intent(out)   :: coord(3,6)
    !
    ! Local variables:
    !
    character(len=54) :: line
    character(len=4)  :: aname
    integer(C_INT)    :: stat
    integer(C_INT)    :: tmp
    logical           :: check(5)

    check=.false.
    !
    coord=0.0
    !
    open(666,file=pdb_in,status='old',iostat=stat)
    !
    if (stat/=0) then
      !
      write(*,*) 'Cannot read the input PDB file!'
      !
      call exit(1)
      !
    end if
    !
    do
      !
      read(666,fmt='(A54)',iostat=stat) line
      !
      if (stat/=0) then
         !
         exit
         !
      else
         !
         if (line(1:6)=='ATOM  '.or.line(1:6)=='HETATM') then
            !
            aname=trim(adjustl(line(13:16)))
            !
            if (aname=='HA  ') then
               !
               tmp=5
               !
             else if (aname=='CB  ') then
               !
               tmp=4
               !
            else if (aname=='CA  ') then
               !
               tmp=3
               !
            else if (aname=='N   ') then
               !
               tmp=2
               !
            else if (aname=='C   ') then
               !
               tmp=1
               !
            else
               !
               tmp=0
               !
            end if
            !
            if (tmp.gt.0) then
               !
               ! Reads the x, y, and z-coordinates.
               !
               read(line(31:38),fmt='(F8.3)') coord(1,tmp)
               !
               read(line(39:46),fmt='(F8.3)') coord(2,tmp)
               !
               read(line(47:54),fmt='(F8.3)') coord(3,tmp)
               !
            end if
            !
         end if
         !
         check(tmp)=.true.
         !
      end if
      !
    end do
    !
    close(666)
    !
    if (any(check==.false.)) then
      !
      print *
      print *,'FATAL ERROR: The PDB file:'
      print *
      print *,trim(adjustl(pdb_in))
      print *
      print *,'does not contain a description of one of the atoms '//&
              '"HA", "CA", "CB", "C", "N"!'
      print *
      !
    end if

    end subroutine read_cg_atom_entries


    function get_atom_type(atom_name) result(res)
    !
    ! Assigns a formal atom type as an integer number. That type helps
    ! to identify the atom.
    !
    ! Interface variables:
    !
    character(len=4),intent(in) :: atom_name
    integer(C_INT)              :: res
    !
    ! Local variables:
    !
    character(len=4)            :: tmp

    tmp=trim(adjustl(atom_name))
    !
    if (tmp=='HA  ') then
      !
      res=5
      !
    else if (tmp=='CB  ') then
      !
      res=4
      !
    else if (tmp=='CA  ') then
      !
      res=3
      !
    else if (tmp=='N   ') then
      !
      res=2
      !
    else if (tmp=='C   ') then
      !
      res=1
      !
    else
      !
      res=0
      !
    end if

    end function get_atom_type


end module mod_pdb

