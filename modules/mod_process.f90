module mod_process

   ! Processing the CG topology - obtains bond angles, bond lengths,
   ! rotates atoms, finds principle axes.
   !
   ! Author : Veselin Kolev <vesso.kolev@gmail.com>
   ! Version: 2017120900
   ! License: GPLv2

   use iso_c_binding,only:C_INT,C_FLOAT
   !
   implicit none

contains


    function s_norm02(vector) result(res)
    !
    ! Finds the norm of the vector "vector" in N-dimensional space.
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(in) :: vector(:)
    real(C_FLOAT)            :: res

    res=sqrt(sum(vector**2))

    end function s_norm02


    subroutine s_angle_between_vectors(v1,v2,angle)
    !
    ! Computes the angle between two vectors in 3D:
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(in) :: v1(3)
    real(C_FLOAT),intent(in) :: v2(3)
    real(C_FLOAT),intent(out) :: angle
    !
    ! Local variables:
    !
    real(C_FLOAT) :: vcp(3)
    !
    call s_vector_cross_product(v1,v2,vcp)
    !
    angle=atan2(sqrt(dot_product(vcp,vcp)),dot_product(v1,v2))
    !
    end subroutine s_angle_between_vectors


    subroutine s_rotate_around_axis(vect,rodrigues_matrix)
    !
    ! Rotates the coordinates of the vector "vect" in 3D space by
    ! implementing the Rodrugues' formula. The Rodrugues' rotation
    ! matrix should be computed in advance by calling the
    ! subroutine get_rodrigues_matrix.
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: vect(3)
    real(C_FLOAT),intent(in)    :: rodrigues_matrix(3,3)
    !
    ! Local variables:
    !
    real(C_FLOAT)  :: rotated(3)
    integer(C_INT) :: i
    integer(C_INT) :: j
    !
    do i=1,3
      !
      rotated(i)=0.0_C_FLOAT
      !
      do j=1,3
         !
         ! Note the index order when doing the matrix to vector dot
         ! multiplication. j is the first index of the matrix
         ! elements. This is beacause of the way Fortran access the
         ! elements of multidimensional arrays.
         !
         rotated(i)=rotated(i)+rodrigues_matrix(j,i)*vect(j)
         !
      end do
      !
    end do
    !
    vect=rotated
    !
    end subroutine s_rotate_around_axis


    subroutine s_vector_cross_product(v1,v2,prod)
    !
    ! Computes the cross product of two vectors in 3D.
    !
    ! Interface varirables:
    !
    real(C_FLOAT),intent(in) :: v1(3)
    real(C_FLOAT),intent(in) :: v2(3)
    real(C_FLOAT),intent(out) :: prod(3)
    !
    ! Local variables:
    !
    integer(C_INT),parameter :: factor1(3)=(/2,3,1/)
    integer(C_INT),parameter :: factor2(3)=(/3,1,2/)
    !
    ! Total vectorization:
    !
    prod(:)=v1(factor1(:))*v2(factor2(:))-v1(factor2(:))*v2(factor1(:))
    !
    end subroutine s_vector_cross_product


    subroutine get_dummy_atom_pos(coord)
    !
    ! Finds the position of dummy atom 'Q', located in the middle
    ! between 'N' and 'C'. Of the distances 'N'-'CA' and 'C'-'NA'
    ! are not equal, equalize those two distances only inside this
    ! subroutine (do not transform the coordinates of either 'N'
    ! or 'CA').
    !
    !
    !           CA
    !         .    .
    !        .      .
    !       .        .
    !      .          .
    !    N      Q      C
    !
    ! The coordinates of 'Q' will appear as coord(:,6).
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout),target :: coord(3,6)
    !
    ! Local variables:
    !
    real(C_FLOAT) :: l1
    real(C_FLOAT) :: l2
    real(C_FLOAT),pointer :: q1(:)
    real(C_FLOAT),pointer :: q2(:)

    l1=s_norm02(coord(:,2)-coord(:,3))
    !
    l2=s_norm02(coord(:,1)-coord(:,3))
    !
    if (l1>l2) then
      !
      l1=l1/l2
      !
      q1=>null()
      q2=>null()
      !
      q1=>coord(:,1)
      q2=>coord(:,2)
      !
    else
      !
      l1=l2/l1
      !
      q1=>null()
      q2=>null()
      !
      q1=>coord(:,2)
      q2=>coord(:,1)
      !
    end if
    !
    ! Compute the coordinates of the dummy atom 'Q':
    !
    coord(:,6)=coord(:,3)+l1*(q1-coord(:,3))
    !
    coord(:,6)=coord(:,6)+0.5*(q2-coord(:,6))

    end subroutine get_dummy_atom_pos


    subroutine s_get_equation_of_plane(v1,v2,v3,equ)
    !
    ! Computes the equation of a 3D plane as:
    !
    ! equ(1)*x+equ(2)*y+equ(3)*z+equ(4)
    !
    ! The parameters are stored as the elements of the array "equ".
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(in) :: v1(3)
    real(C_FLOAT),intent(in) :: v2(3)
    real(C_FLOAT),intent(in) :: v3(3)
    real(C_FLOAT),intent(out) :: equ(4)

    call s_vector_cross_product(v2-v1,v2-v3,equ(1:3))
    !
    equ(4)=v1(1)*equ(1)+v1(2)*equ(2)+v1(3)*equ(3)
    equ(4)=-equ(4)

    end subroutine s_get_equation_of_plane


    subroutine rotate_vect_around_axis(vect,rodrigues_matrix)
    !
    ! Rotates 3D vector around 3D axis.
    !
    ! NOTES: "vect" is the vector to rotate. The matrix named
    !        "rodrigues_matrix" can be computined by calling the
    !         subroutine "get_rodrigues_matrix".
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: vect(3)
    real(C_FLOAT),intent(in)    :: rodrigues_matrix(3,3)
    !
    ! Local variables:
    !
    real(C_FLOAT)  :: rotated(3)
    integer(C_INT) :: i

    do i=1,3
      !
      ! Total vectorization of the second loop inside the sum operator:
      !
      rotated(i)=sum(rodrigues_matrix(:,i)*vect(:))
      !
    end do
    !
    vect=rotated(:)

    end subroutine rotate_vect_around_axis


    subroutine get_rodrigues_matrix(vect,angle,matrix)
    !
    ! Calculates the rotation matrix elements that allow the
    ! application of the Rodrigues' formula. For more details see:
    ! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    !
    ! NOTES: "vect" defines the 3D equation of the axis line.
    !        It MUST begin at (0,0,0). Supply the angle "angle"
    !        in radians!
    !
    ! IMPORTANT!!!
    !
    ! The matrix indexing foolows the Fortran way, not the
    ! one of C. Do take this into account! See for example the
    ! subroutine rotate_around_axis.
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(in)  :: vect(3)
    real(C_FLOAT),intent(in)  :: angle
    real(C_FLOAT),intent(out) :: matrix(3,3)
    !
    ! Local variables:
    !
    real(C_FLOAT) :: norm(3)
    real(C_FLOAT) :: cos_
    real(C_FLOAT) :: sin_
    real(C_FLOAT) :: omcos_
    real(C_FLOAT) :: n12
    real(C_FLOAT) :: n13
    real(C_FLOAT) :: n23

    norm=vect/s_norm02(vect)
    !
    cos_=cos(angle)
    sin_=sin(angle)
    !
    omcos_=1-cos_
    !
    n12=norm(1)*norm(2)
    n13=norm(1)*norm(3)
    n23=norm(2)*norm(3)
    !
    matrix(1,1)=cos_+norm(1)*norm(1)*omcos_
    matrix(2,1)=n12*omcos_-norm(3)*sin_
    matrix(3,1)=n13*omcos_+norm(2)*sin_
    !
    matrix(1,2)=n12*omcos_+norm(3)*sin_
    matrix(2,2)=cos_+norm(2)*norm(2)*omcos_
    matrix(3,2)=n23*omcos_-norm(1)*sin_
    !
    matrix(1,3)=n13*omcos_-norm(2)*sin_
    matrix(2,3)=n23*omcos_+norm(1)*sin_
    matrix(3,3)=cos_+norm(3)*norm(3)*omcos_

    end subroutine get_rodrigues_matrix


    subroutine canonize_residue_coords(coord,ca,angle_2d,angle_3d,axis)
    !
    ! Finds the equation of the plane defined by the atoms "N", "C", and
    ! "CA".
    ! Then translates the atom coordinates to make the "CA" positioned
    ! at (0,0,0).
    ! After that it rotates the atoms until the vector normal to the
    ! plane
    ! becomes (0,0,1). Right after that the subroutine computes the
    ! middle point
    ! between the normalized coordinates of "N" and "C". That point
    ! coordinates
    ! are stored as coord(:,6). Finally, it rotates the atoms around
    ! (0,0,1)
    ! until coord(:,6) becomes (a,0,0), where a>0.
    !
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: coord(3,6)
    real(C_FLOAT),intent(out)   :: ca(3)
    real(C_FLOAT),intent(out)   :: angle_2d
    real(C_FLOAT),intent(out)   :: angle_3d
    real(C_FLOAT),intent(out)   :: axis(3)
    !
    ! Local variables:
    !
    integer(C_INT) :: i
    real(C_FLOAT)  :: l1
    real(C_FLOAT)  :: l2
    real(C_FLOAT)  :: tmp(4)
    real(C_FLOAT)  :: n(3)
    real(C_FLOAT)  :: matrix(3,3)

    ca=coord(:,3)
    !
    ! Center the structure at "CA" atom so that "CA" coordinates
    ! to become (0,0,0):
    !
    do i=1,5
      !
      coord(:,i)=coord(:,i)-ca
      !
    end do
    !
    ! Compute the equation of plane containing "C", "NA", anc "CA"
    ! atoms.
    ! The parameters of the equation are stored as components of q1.
    !
    call s_get_equation_of_plane(coord(:,2),coord(:,3),coord(:,1),tmp)
    !
    ! Get the vector normal to the plane:
    !
    n=tmp(1:3)/s_norm02(tmp(1:3))
    !
    ! Compute the angle between the normal vector and vector (0,0,1):
    !
    call s_angle_between_vectors(n,[0.0,0.0,1.0],angle_3d)
    !
    ! Get the axis for rotating the to bring the normal vector to
    ! (0,0,1):
    !
    call s_vector_cross_product(n,[0.0,0.0,1.0],axis)
    !
    ! Compute the components of the Rodrigues' matrix:
    !
    call get_rodrigues_matrix(axis,-angle_3d,matrix)
    !
    ! and use them to rotate the atom coordinates:
    !
    do i=1,5
      !
      call rotate_vect_around_axis(coord(:,i),matrix)
      !
    end do
    !
    ! Find the middle point - the one in the middle between the
    ! normalized vectors "CA" - "N" and "CA" - "C":
    !
    ! Normalize the vector "CA" - "C":
    !
    tmp(1:3)=coord(:,1)/s_norm02(coord(:,1))
    !
    ! Compute the coordinates of the middle point and implicitly
    ! normalize the vector "CA" - "N":
    !
    coord(:,6)=tmp(1:3)+0.5*(coord(:,2)/s_norm02(coord(:,2))-tmp(1:3))
    !
    ! Rotate the atoms around z-axis, so that the middle point
    ! y and z-coordinates are zeros:
    !
    angle_2d=atan2(coord(2,6),coord(1,6))
    !
    ! Compute the cosine and sinus functions on advance to save
    ! CPU time later when their values are repeatedly used:
    !
    l1=cos(-angle_2d)
    l2=sin(-angle_2d)
    !
    ! Do the rotation around axis z:
    !
    do i=1,6
      !
      ! Copy the coordinates of the atom before starting with the
      ! rotation:
      !
      tmp(1:3)=coord(:,i)
      !
      coord(1,i)=tmp(1)*l1-tmp(2)*l2
      coord(2,i)=tmp(1)*l2+tmp(2)*l1
      !
    end do

    end subroutine canonize_residue_coords


    subroutine canonize_back_residue_coords(coord,ca,angle_2d,angle_3d,axis)
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: coord(3,6)
    real(C_FLOAT),intent(in)    :: ca(3)
    real(C_FLOAT),intent(in)    :: angle_2d
    real(C_FLOAT),intent(in)    :: angle_3d
    real(C_FLOAT),intent(in)    :: axis(3)
    !
    ! Local variables:
    !
    integer(C_INT) :: i
    real(C_FLOAT)  :: l1
    real(C_FLOAT)  :: l2
    real(C_FLOAT)  :: tmp(3)
    real(C_FLOAT)  :: matrix(3,3)

    !
    ! Compute the cosine and sinus functions on advance to save
    ! CPU time later when their values are repeatedly used:
    !
    l1=cos(angle_2d)
    l2=sin(angle_2d)
    !
    ! Do the reverse rotation around axis z:
    !
    do i=1,6
      !
      ! Copy the coordinates of the atom before starting with the
      ! rotation:
      !
      tmp=coord(:,i)
      !
      coord(1,i)=tmp(1)*l1-tmp(2)*l2
      coord(2,i)=tmp(1)*l2+tmp(2)*l1
      !
    end do
    !
    ! Compute the components of the Rodrigues' matrix:
    !
    call get_rodrigues_matrix(axis,angle_3d,matrix)
    !
    ! and use them to rotate the atom coordinates:
    !
    do i=1,5
      !
      call rotate_vect_around_axis(coord(:,i),matrix)
      !
    end do
    !
    ! Translate back the coordinates with respect to the coordinates
    ! of the "CA" atom:
    !
    do i=1,5
      !
      coord(:,i)=coord(:,i)+ca
      !
    end do

    end subroutine canonize_back_residue_coords


    subroutine get_angles_bond_len_pres(coord,cb,ha)
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: coord(3,6)
    real(C_FLOAT),intent(out)   :: cb(3)
    real(C_FLOAT),intent(out)   :: ha(3)
    !
    ! Local variables:
    !
    integer(C_INT) :: i
    real(C_FLOAT)  :: l1
    real(C_FLOAT)  :: l2
    real(C_FLOAT)  :: tmp(3)

    !
    ! Start collecting the Euler angles for keeping the relative
    ! position of "CB":
    !
    ! 1. Get y=0.0 by rotating "CB" around z-axis.
    !
    ! 2. Get x=0.0 by rotating "CB" around x-axis.
    !
    ! 3. Get the bond length "CA" - "CB" from the z-position of
    !    "CB" atom.
    !
    ! Rotation around z-axis:
    !
    cb(1)=atan2(coord(2,4),coord(1,4))
    !
    ! Compute the cosine and sinus functions on advance to save
    ! CPU time later when their values are repeatedly used:
    !
    l1=cos(-cb(1))
    l2=sin(-cb(1))
    !
    do i=4,5
      !
      tmp=coord(:,i)
      !
      coord(1,i)=tmp(1)*l1-tmp(2)*l2
      coord(2,i)=tmp(1)*l2+tmp(2)*l1
      !
    end do
    !
    ! Rotation around y-axis:
    !
    cb(2)=atan2(coord(1,4),coord(3,4))
    !
    ! Compute the cosine and sinus functions on advance to save
    ! CPU time later when their values are repeatedly used:
    !
    l1=cos(-cb(2))
    l2=sin(-cb(2))
    !
    do i=4,5
      !
      tmp=coord(:,i)
      !
      coord(1,i)= tmp(1)*l1+tmp(3)*l2
      coord(3,i)=-tmp(1)*l2+tmp(3)*l1
      !
    end do
    !
    ! This is in fact the bond length "CA" - "CB":
    !
    cb(3)=coord(3,4)
    !
    ha=coord(:,5)

    end subroutine get_angles_bond_len_pres


    subroutine reconstruct_atom_pos(coord,cb,ha)
    !
    ! Interface variables:
    !
    real(C_FLOAT),intent(inout) :: coord(3,6)
    real(C_FLOAT),intent(in)    :: cb(3)
    real(C_FLOAT),intent(in)    :: ha(3)
    !
    ! Local variables:
    !
    integer(C_INT) :: i
    real(C_FLOAT)  :: l1
    real(C_FLOAT)  :: l2
    real(C_FLOAT)  :: tmp(3)
    real(C_FLOAT)  :: tmp_coord(3,2)

    coord(:,4)=(/0.0,0.0,cb(3)/)
    coord(:,5)=ha

    !
    ! Reconstruct the atom positions by rotating back
    ! around y-axis:
    !
    ! Compute the cosine and sinus functions on advance to save
    ! CPU time later when their values are repeatedly used:
    !
    l1=cos(cb(2))
    l2=sin(cb(2))
    !
    do i=4,5
      !
      tmp=coord(:,i)
      !
      coord(1,i)= tmp(1)*l1+tmp(3)*l2
      coord(3,i)=-tmp(1)*l2+tmp(3)*l1
      !
    end do
    !
    ! Reconstruct the atom positions by rotating back
    ! around y-axis:
    !
    l1=cos(cb(1))
    l2=sin(cb(1))
    !
    do i=4,5
      !
      tmp=coord(:,i)
      !
      coord(1,i)=tmp(1)*l1-tmp(2)*l2
      coord(2,i)=tmp(1)*l2+tmp(2)*l1
      !
    end do

    end subroutine reconstruct_atom_pos


end module mod_process

