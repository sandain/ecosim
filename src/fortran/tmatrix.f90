!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
!
!    Copyright (C) 2013  Jason M. Wood, Montana State University
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for storing a triangular matrix as a vector.
!>
!> Reduces the storage space required to store a triangular matrix, from 
!> N^2 when using a standard 2D-array to N*(N+1)/2 when using a vector.
!>
!> #### Usage
!>
!> @code
!>   program tmatrixTest
!>     use tmatrix
!>     implicit none
!>     integer(kind=4)      :: i
!>     integer(kind=4)      :: j
!>     integer(kind=4)      :: k
!>     type(tmatrixInteger) :: integerMatrix
!>     ! Initialize the triangular matrix with an initial size of 10.
!>     call tmatrixInitialize(integerMatrix, 10)
!>     ! Fill the triangular matrix.
!>     k = 1
!>     do i = 1, 10
!>       do j = 1, i
!>         call tmatrixSet(integerMatrix, i, j, k)
!>         k = k + 1
!>       end do
!>     end do
!>     ! Write the triangular matrix to the screen.
!>     write(unit = *, fmt = '(5X10(A1I0.2A1))') ('[',j,']', j = 1, 10)
!>     do i = 1, 10
!>       write(unit = *, fmt ='(A1I0.2A1X10(XI2X))') '[', i, ']', &
!>         (tmatrixGet(integerMatrix, i, j), j = 1, i)
!>     end do
!>     ! Deallocate the triangular matrix.
!>     call tmatrixClose(integerMatrix)
!>   end program tmatrixTest
!> @endcode
!>
!> Expected output:
!>
!>      [01][02][03][04][05][06][07][08][09][10]
!> [01]   1
!> [02]   2   3
!> [03]   4   5   6
!> [04]   7   8   9  10
!> [05]  11  12  13  14  15
!> [06]  16  17  18  19  20  21
!> [07]  22  23  24  25  26  27  28
!> [08]  29  30  31  32  33  34  35  36
!> [09]  37  38  39  40  41  42  43  44  45
!> [10]  46  47  48  49  50  51  52  53  54  55
!>
!> @author Jason M. Wood
!> @copyright GNU General Public License
module tmatrix
  implicit none
  private

  ! Public matrix types.
  public :: tmatrixInteger
  public :: tmatrixLong
  public :: tmatrixReal
  public :: tmatrixDouble

  ! Public methods.
  public :: tmatrixInitialize
  public :: tmatrixClose
  public :: tmatrixGet
  public :: tmatrixSet

  !> The integer tmatrix: integer(kind=4).
  type :: tmatrixInteger
    integer(kind = 4), allocatable :: matrix(:)
    integer(kind = 4)              :: capacity
  end type tmatrixInteger

  !> The long tmatrix: integer(kind=8).
  type :: tmatrixLong
    integer(kind = 8), allocatable :: matrix(:)
    integer(kind = 4)              :: capacity
  end type tmatrixLong

  !> The real tmatrix: real(kind=4).
  type :: tmatrixReal
    real(kind = 4), allocatable    :: matrix(:)
    integer(kind = 4)              :: capacity
  end type tmatrixReal

  !> The double precision tmatrix: real(kind=8).
  type :: tmatrixDouble
    real(kind = 8), allocatable    :: matrix(:)
    integer(kind = 4)              :: capacity
  end type tmatrixDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface tmatrixInitialize
    module procedure tmatrixInitializeInteger
    module procedure tmatrixInitializeLong
    module procedure tmatrixInitializeReal
    module procedure tmatrixInitializeDouble
  end interface tmatrixInitialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface tmatrixClose
    module procedure tmatrixCloseInteger
    module procedure tmatrixCloseLong
    module procedure tmatrixCloseReal
    module procedure tmatrixCloseDouble
  end interface tmatrixClose

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element needed.
  !> @param[in]     j             The j index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface tmatrixGet
    module procedure tmatrixGetInteger
    module procedure tmatrixGetLong
    module procedure tmatrixGetReal
    module procedure tmatrixGetDouble
  end interface tmatrixGet

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element to change.
  !> @param[in]     j             The j index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface tmatrixSet
    module procedure tmatrixSetInteger
    module procedure tmatrixSetLong
    module procedure tmatrixSetReal
    module procedure tmatrixSetDouble
  end interface tmatrixSet

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixInitializeInteger (tmatrix, capacity)
    type(tmatrixInteger), intent(inout) :: tmatrix
    integer(kind=4), intent(in)         :: capacity
    ! Local variables
    integer :: allocateStatus
    integer :: matrixSize
    ! Set the capacity of the triangular matrix.
    tmatrix%capacity = capacity
    ! Allocate space for the triangular matrix.
    matrixSize = capacity * (capacity + 1) / 2
    allocate(tmatrix%matrix(matrixSize), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixInitializeInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixInitializeLong (tmatrix, capacity)
    type(tmatrixLong), intent(inout)    :: tmatrix
    integer(kind=4), intent(in)         :: capacity
    ! Local variables
    integer :: allocateStatus
    integer :: matrixSize
    ! Set the capacity of the triangular matrix.
    tmatrix%capacity = capacity
    ! Allocate space for the triangular matrix.
    matrixSize = capacity * (capacity + 1) / 2
    allocate(tmatrix%matrix(matrixSize), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixInitializeLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixInitializeReal (tmatrix, capacity)
    type(tmatrixReal), intent(inout)    :: tmatrix
    integer(kind=4), intent(in)         :: capacity
    ! Local variables
    integer :: allocateStatus
    integer :: matrixSize
    ! Set the capacity of the triangular matrix.
    tmatrix%capacity = capacity
    ! Allocate space for the triangular matrix.
    matrixSize = capacity * (capacity + 1) / 2
    allocate(tmatrix%matrix(matrixSize), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixInitializeReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixInitializeDouble (tmatrix, capacity)
    type(tmatrixDouble), intent(inout)  :: tmatrix
    integer(kind=4), intent(in)         :: capacity
    ! Local variables
    integer :: allocateStatus
    integer :: matrixSize
    ! Set the capacity of the triangular matrix.
    tmatrix%capacity = capacity
    ! Allocate space for the triangular matrix.
    matrixSize = capacity * (capacity + 1) / 2
    allocate(tmatrix%matrix(matrixSize), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixInitializeDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixCloseInteger (tmatrix)
    type(tmatrixInteger), intent(inout) :: tmatrix
    ! Local variables.
    integer :: allocateStatus
    ! Deallocate space for the triangular matrix.
    deallocate(tmatrix%matrix, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixCloseInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixCloseLong (tmatrix)
    type(tmatrixLong), intent(inout)    :: tmatrix
    ! Local variables.
    integer :: allocateStatus
    ! Deallocate space for the triangular matrix.
    deallocate(tmatrix%matrix, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixCloseLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixCloseReal (tmatrix)
    type(tmatrixReal), intent(inout)    :: tmatrix
    ! Local variables.
    integer :: allocateStatus
    ! Deallocate space for the triangular matrix.
    deallocate(tmatrix%matrix, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixCloseReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixCloseDouble (tmatrix)
    type(tmatrixDouble), intent(inout)  :: tmatrix
    ! Local variables.
    integer :: allocateStatus
    ! Deallocate space for the triangular matrix.
    deallocate(tmatrix%matrix, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the triangular matrix!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine tmatrixCloseDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Private method to calculate the index of the vector from the size and
  !> i, j indices of the triangular matrix.
  !>
  !> @param[in]     capacity      The capacity of the triangular matrix.
  !> @param[in]     i             The i index of the matrix element needed.
  !> @param[in]     j             The j index of the matrix element needed.
  !> @return                      The index of the vector.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tmatrixGetIndex (capacity, i, j) result (returnValue)
    integer(kind=4), intent(in)         :: capacity
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    integer(kind=4)                     :: returnValue
    ! Make sure the provided index is valid.
    if (i .lt. 1 .or. j .lt. 1 .or. &
      i .gt. capacity .or. j .gt. capacity) then
      write (unit = *, fmt = *) &
        "Invalid index to the triangular matrix!", i, j
      ! Error, exit the program.
      stop
    end if
    ! Calculate the index.
    if (i .lt. j) then
      returnValue = (i - 1) * capacity - (i - 1) * i / 2 + j
    else
      returnValue = (j - 1) * capacity - (j - 1) * j / 2 + i
    end if
    return
  end function tmatrixGetIndex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the triangular matrix.
  !>
  !> @param[in]     tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element needed.
  !> @param[in]     j             The j index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tmatrixGetInteger (tmatrix, i, j) result (returnValue)
    type(tmatrixInteger), intent(in)    :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    integer(kind=4)                     :: returnValue
    ! Grab the element at the provided index.
    returnValue = tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j))
    return
  end function tmatrixGetInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the triangular matrix.
  !>
  !> @param[in]     tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element needed.
  !> @param[in]     j             The j index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tmatrixGetLong (tmatrix, i, j) result (returnValue)
    type(tmatrixLong), intent(in)       :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    integer(kind=8)                     :: returnValue
    ! Grab the element at the provided index.
    returnValue = tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j))
    return
  end function tmatrixGetLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the triangular matrix.
  !>
  !> @param[in]     tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element needed.
  !> @param[in]     j             The j index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tmatrixGetReal (tmatrix, i, j) result (returnValue)
    type(tmatrixReal), intent(in)       :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    real(kind=4)                        :: returnValue
    ! Grab the element at the provided index.
    returnValue = tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j))
    return
  end function tmatrixGetReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the triangular matrix.
  !>
  !> @param[in]     tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element needed.
  !> @param[in]     j             The j index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tmatrixGetDouble (tmatrix, i, j) result (returnValue)
    type(tmatrixDouble), intent(in)     :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    real(kind=8)                        :: returnValue
    ! Grab the element at the provided index.
    returnValue = tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j))
    return
  end function tmatrixGetDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element to change.
  !> @param[in]     j             The j index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixSetInteger (tmatrix, i, j, input)
    type(tmatrixInteger), intent(inout) :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    integer(kind=4), intent(in)         :: input
    ! Set the element at the provided index.
    tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j)) = input
    return
  end subroutine tmatrixSetInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element to change.
  !> @param[in]     j             The j index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixSetLong (tmatrix, i, j, input)
    type(tmatrixLong), intent(inout)    :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    integer(kind=8), intent(in)         :: input
    ! Set the element at the provided index.
    tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j)) = input
    return
  end subroutine tmatrixSetLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element to change.
  !> @param[in]     j             The j index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixSetReal (tmatrix, i, j, input)
    type(tmatrixReal), intent(inout)    :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    real(kind=4), intent(in)            :: input
    ! Set the element at the provided index.
    tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j)) = input
    return
  end subroutine tmatrixSetReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the triangular matrix.
  !>
  !> @param[inout]  tmatrix       The triangular matrix.
  !> @param[in]     i             The i index of the element to change.
  !> @param[in]     j             The j index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tmatrixSetDouble (tmatrix, i, j, input)
    type(tmatrixDouble), intent(inout)  :: tmatrix
    integer(kind=4), intent(in)         :: i
    integer(kind=4), intent(in)         :: j
    real(kind=8), intent(in)            :: input
    ! Set the element at the provided index.
    tmatrix%matrix(tmatrixGetIndex(tmatrix%capacity, i, j)) = input
    return
  end subroutine tmatrixSetDouble

end module tmatrix
