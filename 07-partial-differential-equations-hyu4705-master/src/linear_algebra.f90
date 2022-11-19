!-----------------------------------------------------------------------
!Module: linear_algebra
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! These subroutines and functions are used to compute the simple linear algebra like ax=b and inverse of matrix a.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! ...
!!----------------------------------------------------------------------
module linear_algebra
use types
implicit none
private
public :: solve_linear_system, invert_matrix, solve_equations
contains

!-----------------------------------------------------------------------
!! Subroutine: solve_linear_system
!-----------------------------------------------------------------------
!! by: Heng Yu
!!
!! Take input of a matrix and b vector to compute x vector where: ax = b
!!----------------------------------------------------------------------
!! Input:
!!
!! a_matrinx        real        2D array containing the $a$ matrix
!! b_vector         real        1D array containing the $b$ vector
!-----------------------------------------------------------------------
!! Output:
!!
!! x_vector         real        1D array with the solution to the system of equations
!! a_inverse        real        2D array with the inverse matrix $a^{-1}$
!-----------------------------------------------------------------------
subroutine solve_linear_system(a_matrix, b_vector, x_vector, a_inverse)
    implicit none
    real(dp), intent(in) :: a_matrix(:,:), b_vector(:)
    real(dp), intent(out) ::  x_vector(:), a_inverse(:,:)
    ! The first thing is to make sure that all the arrays have proper sizes to
    ! make sure that  the matrix and vectors provided actually represent a
    ! system of equations. Otherwise the subroutines below will give errors or
    ! worst won't behave as expected but we won't notice.
    call test_array_sizes(a_matrix, b_vector, a_inverse, x_vector)
    call invert_matrix(a_matrix, a_inverse)
    
    ! Now that you have the inverse of a_matrix, how can you use a_inverse to 
    ! solve the system of equations? 
    
    ! Write a subroutine for that purpose and call it here (Don't forget to
    ! document the subroutine)
    ! call ...
    call solve_equations(b_vector, a_matrix, x_vector)

end subroutine solve_linear_system

subroutine solve_equations(b_vector, a_inverse, x_vector)
    implicit none
    real(dp), intent(in) :: b_vector(:), a_inverse(:,:)
    real(dp), intent(out) :: x_vector(:)
    real(dp) :: sum
    integer :: i, j

    do i=1, size(x_vector)
        sum = 0._dp

        do j=1, size(b_vector)
            sum = sum + a_inverse(i,j)*b_vector(j)
        enddo

        x_vector(i) = sum
    enddo
    
end subroutine solve_equations
!-----------------------------------------------------------------------
!! Subroutine: solve_linear_system
!-----------------------------------------------------------------------
!! By Heng Yu
!!
!! The subroutine is to test the size of a, b, x matrix/vector to make sure the processes are succesful.
!!----------------------------------------------------------------------
!! Input:
!!
!! a_matrinx        real        2D array containing the $a$ matrix
!! b_vector         real        1D array containing the $b$ vector
!! x_vector         real        1D array with the solution to the system of equations
!! a_inverse        real        2D array with the inverse matrix $a^{-1}$
!-----------------------------------------------------------------------
subroutine test_array_sizes(a_matrix, b_vector, a_inverse, x_vector)
    implicit none
    real(dp), intent(in) :: a_matrix(:,:), b_vector(:), a_inverse(:,:), x_vector(:)

    integer :: shape_a(1:2), b_size, shape_ainv(1:2), x_size
    
    shape_a = shape(a_matrix)
    shape_ainv = shape(a_inverse)
    b_size = size(b_vector)
    x_size = size(x_vector)

    ! Test that a_matrix and a_inverse are square and of the same size
    if (shape_a(1) /= shape_ainv(1) .OR. shape_a(2) /= shape_ainv(2)) then
        print *, 'Notice: the shape of a_matrix and a_inverse are not same'
        stop
    else if(shape_a(1) /= shape_a(2)) then
        print *, 'Notice: the shape of a matrix is not square'
        stop
    endif
    ! test that the number of columns in a_matrix is equal to the number of elements in x_vector
    if (shape_a(2) /= x_size) then
        print *, 'Notice: the a_matrix columns and size of x_vector are not same'
        stop
    endif
    ! test that the number of rows in a_matrix is equal to the number of elements in b_vector
    if (shape_a(1) /= b_size) then
        print *, 'Notice: the a_matrix rows and size of b_vector are not same'
        stop
    endif
end subroutine test_array_sizes

!-----------------------------------------------------------------------
!! Subroutine: invert_matrix
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Given a non singular matrix $a$, returns its inverse $a^{-1}$
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real    2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a_inv    real    2D array with the $a^{-1}$ matrix
!-----------------------------------------------------------------------
subroutine invert_matrix(a, a_inv)
    implicit none
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: a_inv(:,:)
    real(dp), allocatable :: a_work(:,:)
    integer :: shape_a(1:2), n, i
    real(dp) :: d
    integer, allocatable :: indx(:)

    ! This I'll give you for free. It's the LU decomposition we discussed in
    ! class with the following back-substitution
    
    allocate(a_work,mold=a)
    shape_a = shape(a)
    n = shape_a(1)
    allocate(indx(1:n))
    
    ! ludcmp destroys the input matrix a. In order to preserve a we will copy
    ! it into a work array that will be used in ludcmp
    a_work = a
    call ludcmp(a_work,indx,d)
    
    ! We construct a matrix that has orthogonal unit vectors as columns
    a_inv = 0._dp
    do i=1,n
        a_inv(i,i) = 1._dp
    enddo

    ! And then feed each column to the back-substitution routine
    do i = 1,n
        call lubksb(a_work,indx,a_inv(:,i))
    enddo
    ! This results in a_inv being the inverse of a
end subroutine invert_matrix

! The subroutines below were taken from numerical recipes and were slightly
! modified  to work with double precision reals.

! Notice how much harder it is to understand what a code does when 
! explicit informative names are not used for the different variables
! and processes 

!-----------------------------------------------------------------------
!! Subroutine: ludcmp
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Adapted from numerical recipes subroutine.
!! Performs LU decomposition on a non singular matrix $a$.
!! The original $a$ matrix is destroyed as the LU decomposition is returned
!! in the same array
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a        real        2D array with LU decomposition of the $a$ matrix
!! indx     integer     1D array that records the row permutation effected by the partial pivoting
!! d        real        +1 or -1 depending on whether the number of row interchanges was even or odd, respectively
!-----------------------------------------------------------------------
subroutine ludcmp(a, indx, d)
    implicit none
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: indx(:)
    real(dp), intent(out) :: d
    integer :: n,i,imax,j,k
    real(dp) aamax,dum,sum
    real(dp), allocatable :: vv(:)
    n = size(indx)
    allocate(vv(1:n))
    d=1._dp
    do i=1,n
        aamax=0._dp
        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0._dp) then
            print *, 'singular matrix in ludcmp'
            stop
        endif
        vv(i)=1._dp/aamax
    enddo

    do j=1,n
        do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
        enddo
        aamax=0._dp
        do i=j,n
            sum=a(i,j)
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
                imax=i
                aamax=dum
            endif
        enddo
        if (j.ne.imax)then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0._dp) a(j,j) = tiny(1._sp)
        if(j.ne.n)then
            dum=1._dp/a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            enddo
        endif
    enddo
end subroutine ludcmp

!-----------------------------------------------------------------------
!! Subroutine: lubksb
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Adapted from numerical recipes subroutine.
!!
!! Performs back-substitution after a LU decomposition in order to solve the
!! linear system of equations $a \cdot x = b$. The $b$ vector is given in the b
!! array (which is destroyed) and the solution $x$ is returned in its place
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the LU decomposition $a$ matrix (as returned by ludecomp)
!! indx     integer     1D array with the record of the row permutation effected by the partial pivoting (as returned by ludecomp)
!! b        real        1D array containing the $b$ vector
!!----------------------------------------------------------------------
!! Output:
!! b        real        1D array containing the $x$ vector
!-----------------------------------------------------------------------
subroutine lubksb(a, indx, b)
    implicit none
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: indx(:)
    real(dp), intent(inout) :: b(:)

    integer :: n
    integer :: i,ii,j,ll
    real(dp) :: sum

    n = size(b)
    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            enddo
        else if (sum.ne.0.) then
            ii=i
        endif
        b(i)=sum
    enddo
    do i=n,1,-1
        sum=b(i)
        do j=i+1,n
            sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
    enddo
end subroutine lubksb

end module linear_algebra