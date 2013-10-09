!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
!
!    Copyright (C) 2013  Jason M. Wood, Montana State University
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
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
!> Module for dynamic arrays.
!>
!> #### Usage
!>
!> @code
!>   program darrayTest
!>     use darray
!>     implicit none
!>     integer(kind=4)     :: i
!>     type(darrayInteger) :: integerArray
!>     ! Initialize the dynamic array with an initial size of 5.
!>     call darrayInitialize(integerArray, 5)
!>     ! Write more than 5 integers to the dynamic array.
!>     do i = 1, 100
!>       call darraySet(integerArray, i, i * 23)
!>     end do
!>     ! Write the dynamic array to the screen.
!>     do i = 1, 100
!>       write(*,*) i, darrayGet(integerArray, i)
!>     end do
!>     ! Deallocate the dynamic array.
!>     call darrayClose(integerArray)
!>   end program darrayTest
!> @endcode
!>
!> @author Jason M. Wood
!> @copyright GNU General Public License
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module darray
  implicit none
  private

  ! Public array types.
  public :: darrayInteger
  public :: darrayLong
  public :: darrayReal
  public :: darrayDouble

  ! Public methods.
  public :: darrayInitialize
  public :: darrayClose
  public :: darrayGet
  public :: darraySet

  !> The integer darray: integer(kind=4).
  type :: darrayInteger
    integer(kind = 4), allocatable :: array(:)
    integer(kind = 4)              :: capacity
    integer(kind = 4)              :: initialCapacity
  end type darrayInteger

  !> The long darray: integer(kind=8).
  type :: darrayLong
    integer(kind = 8), allocatable :: array(:)
    integer(kind = 4)              :: capacity
    integer(kind = 4)              :: initialCapacity
  end type darrayLong

  !> The real darray: real(kind=4).
  type :: darrayReal
    real(kind = 4), allocatable    :: array(:)
    integer(kind = 4)              :: capacity
    integer(kind = 4)              :: initialCapacity
  end type darrayReal

  !> The double precision darray: real(kind=8).
  type :: darrayDouble
    real(kind = 8), allocatable    :: array(:)
    integer(kind = 4)              :: capacity
    integer(kind = 4)              :: initialCapacity
  end type darrayDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     capacity      The initial capacity of the dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface darrayInitialize
    module procedure darrayInitializeInteger
    module procedure darrayInitializeLong
    module procedure darrayInitializeReal
    module procedure darrayInitializeDouble
  end interface darrayInitialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface darrayClose
    module procedure darrayCloseInteger
    module procedure darrayCloseLong
    module procedure darrayCloseReal
    module procedure darrayCloseDouble
  end interface darrayClose

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface darrayGet
    module procedure darrayGetInteger
    module procedure darrayGetLong
    module procedure darrayGetReal
    module procedure darrayGetDouble
  end interface darrayGet

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface darraySet
    module procedure darraySetInteger
    module procedure darraySetLong
    module procedure darraySetReal
    module procedure darraySetDouble
  end interface darraySet

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     capacity      The initial capacity of the dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayInitializeInteger (darray, capacity)
    type(darrayInteger), intent(inout) :: darray
    integer, intent(in)                :: capacity
    ! Local variables
    integer :: allocate_status
    ! Set the capacity of the dynamic array.
    darray%capacity = capacity
    darray%initialCapacity = capacity
    ! Allocate space for the dynamic array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayInitializeInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     capacity      The initial capacity of the dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayInitializeLong (darray, capacity)
    type(darrayLong), intent(inout)    :: darray
    integer, intent(in)                :: capacity
    ! Local variables
    integer :: allocate_status
    ! Set the capacity of the dynamic array.
    darray%capacity = capacity
    darray%initialCapacity = capacity
    ! Allocate space for the dynamic array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayInitializeLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     capacity      The initial capacity of the dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayInitializeReal (darray, capacity)
    type(darrayReal), intent(inout)    :: darray
    integer, intent(in)                :: capacity
    ! Local variables
    integer :: allocate_status
    ! Set the capacity of the dynamic array.
    darray%capacity = capacity
    darray%initialCapacity = capacity
    ! Allocate space for the dynamic array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayInitializeReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize memory for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     capacity      The initial capacity of the dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayInitializeDouble (darray, capacity)
    type(darrayDouble), intent(inout)  :: darray
    integer, intent(in)                :: capacity
    ! Local variables
    integer :: allocate_status
    ! Set the capacity of the dynamic array.
    darray%capacity = capacity
    darray%initialCapacity = capacity
    ! Allocate space for the dynamic array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayInitializeDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for an array.
  !>
  !> @param[inout]  darray        The dynamic array to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayCloseInteger (darray)
    type(darrayInteger), intent(inout) :: darray
    ! Local variables.
    integer :: allocate_status
    ! Deallocate space for the dynamic array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayCloseInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for an array.
  !>
  !> @param[inout]  darray        The dynamic array to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayCloseLong (darray)
    type(darrayLong), intent(inout)    :: darray
    ! Local variables.
    integer :: allocate_status
    ! Deallocate space for the dynamic array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayCloseLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for a array.
  !>
  !> @param[inout]  darray        The dynamic array to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayCloseReal (darray)
    type(darrayReal), intent(inout)    :: darray
    ! Local variables.
    integer :: allocate_status
    ! Deallocate space for the dynamic array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayCloseReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate memory reserved for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array to deallocate.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayCloseDouble (darray)
    type(darrayDouble), intent(inout)  :: darray
    ! Local variables.
    integer :: allocate_status
    ! Deallocate space for the dynamic array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine darrayCloseDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the dynamic array.
  !>
  !> @param[in]     darray        The dynamic array.
  !> @param[in]     index         The index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function darrayGetInteger (darray, index) result (returnValue)
    type(darrayInteger), intent(in)    :: darray
    integer(kind=4), intent(in)        :: index
    integer(kind=4)                    :: returnValue
    ! Grab the element at the provided index.
    if (index .le. darray%capacity) then
      returnValue = darray%array(index)
    else
      write (unit = *, fmt = *) &
        "Invalid index to the dynamic array!", index
      ! Error, exit the program.
      stop
    end if
    return
  end function darrayGetInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the dynamic array.
  !>
  !> @param[in]     darray        The dynamic array.
  !> @param[in]     index         The index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function darrayGetLong (darray, index) result (returnValue)
    type(darrayLong), intent(in)       :: darray
    integer(kind=4), intent(in)        :: index
    integer(kind=8)                    :: returnValue
    ! Grab the element at the provided index.
    if (index .le. darray%capacity) then
      returnValue = darray%array(index)
    else
      write (unit = *, fmt = *) &
        "Invalid index to the dynamic array!", index
      ! Error, exit the program.
      stop
    end if
    return
  end function darrayGetLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the dynamic array.
  !>
  !> @param[in]     darray        The dynamic array.
  !> @param[in]     index         The index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function darrayGetReal (darray, index) result (returnValue)
    type(darrayReal), intent(in)       :: darray
    integer(kind=4), intent(in)        :: index
    real(kind=4)                       :: returnValue
    ! Grab the element at the provided index.
    if (index .le. darray%capacity) then
      returnValue = darray%array(index)
    else
      write (unit = *, fmt = *) &
        "Invalid index to the dynamic array!", index
      ! Error, exit the program.
      stop
    end if
    return
  end function darrayGetReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get an element at the provided index of the dynamic array.
  !>
  !> @param[in]     darray        The dynamic array.
  !> @param[in]     index         The index of the element needed.
  !> @return                      The element at the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function darrayGetDouble (darray, index) result (returnValue)
    type(darrayDouble), intent(in)     :: darray
    integer(kind=4), intent(in)        :: index
    real(kind=8)                       :: returnValue
    ! Grab the element at the provided index.
    if (index .le. darray%capacity) then
      returnValue = darray%array(index)
    else
      write (unit = *, fmt = *) &
        "Invalid index to the dynamic array!", index
      ! Error, exit the program.
      stop
    end if
    return
  end function darrayGetDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darraySetInteger (darray, index, input)
    type(darrayInteger), intent(inout) :: darray
    integer(kind=4), intent(in)        :: index
    integer(kind=4), intent(in)        :: input
    ! Make sure that the dynamic array capacity is large enough.
    do while (index .ge. darray%capacity)
      call darrayReallocateInteger(darray)
    end do
    darray%array(index) = input
    return
  end subroutine darraySetInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darraySetLong (darray, index, input)
    type(darrayLong), intent(inout)    :: darray
    integer(kind=4), intent(in)        :: index
    integer(kind=8), intent(in)        :: input
    ! Make sure that the dynamic array capacity is large enough.
    do while (index .ge. darray%capacity)
      call darrayReallocateLong(darray)
    end do
    darray%array(index) = input
    return
  end subroutine darraySetLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darraySetReal (darray, index, input)
    type(darrayReal), intent(inout)    :: darray
    integer(kind=4), intent(in)        :: index
    real(kind=4), intent(in)           :: input
    ! Make sure that the dynamic array capacity is large enough.
    do while (index .ge. darray%capacity)
      call darrayReallocateReal(darray)
    end do
    darray%array(index) = input
    return
  end subroutine darraySetReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set an element at the provided index of the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !> @param[in]     index         The index of the element to change.
  !> @param[in]     input         The element for the provided index.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darraySetDouble (darray, index, input)
    type(darrayDouble), intent(inout)  :: darray
    integer(kind=4), intent(in)        :: index
    real(kind=8), intent(in)           :: input
    ! Make sure that the dynamic array capacity is large enough.
    do while (index .ge. darray%capacity)
      call darrayReallocateDouble(darray)
    end do
    darray%array(index) = input
    return
  end subroutine darraySetDouble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Private subroutine to reallocate the space for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayReallocateInteger (darray)
    type(darrayInteger), intent(inout) :: darray
    ! Local variables
    integer(kind=4)              :: allocate_status
    integer(kind=4)              :: capacity
    integer(kind=4), allocatable :: tmp_array(:)
    ! Increase capacity.
    capacity = darray%capacity + darray%initialCapacity
    ! Allocate space for the temporary array.
    allocate(tmp_array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the data to the temporary array.
    tmp_array(1:darray%capacity) = darray%array
    ! Deallocate the old array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Allocate space for the new array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the temporary array into the newly resized array.
    darray%array = tmp_array
    ! Deallocate the temporary array.
    deallocate(tmp_array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Update the capacity information.
    darray%capacity = capacity
    return
  end subroutine darrayReallocateInteger

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Private subroutine to reallocate the space for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayReallocateLong (darray)
    type(darrayLong), intent(inout)    :: darray
    ! Local variables
    integer(kind=4)              :: allocate_status
    integer(kind=4)              :: capacity
    integer(kind=8), allocatable :: tmp_array(:)
    ! Increase capacity.
    capacity = darray%capacity + darray%initialCapacity
    ! Allocate space for the temporary array.
    allocate(tmp_array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the data to the temporary array.
    tmp_array(1:darray%capacity) = darray%array
    ! Deallocate the old array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Allocate space for the new array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the temporary array into the newly resized array.
    darray%array = tmp_array
    ! Deallocate the temporary array.
    deallocate(tmp_array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Update the capacity information.
    darray%capacity = capacity
    return
  end subroutine darrayReallocateLong

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Private subroutine to reallocate the space for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayReallocateReal (darray)
    type(darrayReal), intent(inout)    :: darray
    ! Local variables
    integer(kind=4)              :: allocate_status
    integer(kind=4)              :: capacity
    real(kind=4), allocatable    :: tmp_array(:)
    ! Increase capacity.
    capacity = darray%capacity + darray%initialCapacity
    ! Allocate space for the temporary array.
    allocate(tmp_array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the data to the temporary array.
    tmp_array(1:darray%capacity) = darray%array
    ! Deallocate the old array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Allocate space for the new array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the temporary array into the newly resized array.
    darray%array = tmp_array
    ! Deallocate the temporary array.
    deallocate(tmp_array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Update the capacity information.
    darray%capacity = capacity
    return
  end subroutine darrayReallocateReal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Private subroutine to reallocate the space for the dynamic array.
  !>
  !> @param[inout]  darray        The dynamic array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine darrayReallocateDouble (darray)
    type(darrayDouble), intent(inout)  :: darray
    ! Local variables
    integer(kind=4)              :: allocate_status
    integer(kind=4)              :: capacity
    real(kind=8), allocatable    :: tmp_array(:)
    ! Increase capacity.
    capacity = darray%capacity + darray%initialCapacity
    ! Allocate space for the temporary array.
    allocate(tmp_array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the data to the temporary array.
    tmp_array(1:darray%capacity) = darray%array
    ! Deallocate the old array.
    deallocate(darray%array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Allocate space for the new array.
    allocate(darray%array(capacity), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the dynamic array!"
      ! Error, exit the program.
      stop
    end if
    ! Copy the temporary array into the newly resized array.
    darray%array = tmp_array
    ! Deallocate the temporary array.
    deallocate(tmp_array, stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to allocate memory for the temporary array!"
      ! Error, exit the program.
      stop
    end if
    ! Update the capacity information.
    darray%capacity = capacity
    return
  end subroutine darrayReallocateDouble

end module darray
