!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
!
!    Copyright (C) 2000  George Marsaglia and Wai Wan Tsang
!    Copyright (C) 2013-2026  Jason M. Wood
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
!> The ziggurat algorithm is a rejection sampling algorithm for pseudo-random
!> number sampling. The original C version was created by George Marsaglia and
!> Wai Wan Tsang [1].  The Fortran 90 version compatible with OpenMP was
!> created by Jason M. Wood.
!>
!> #### Usage
!>
!> @code
!>   program zigguratTest
!>     use, intrinsic :: iso_fortran_env
!>     use :: ziggurat
!>     implicit none
!>     integer, parameter  :: iii = 362436069 !< The RNG seed.
!>     integer             :: i
!>     real(kind = real64) :: x
!>     type(ziggurat_t)    :: rng
!>     ! Seed the RNG.
!>     call ziggurat_seed (rng, iii)
!>     ! Generate some random numbers.
!>     do i = 1, 10
!>       x = ziggurat_uni (rng)
!>       write (*,*) i, x
!>     end do
!>   end program zigguratTest
!> @endcode
!>
!> #### Reference
!>  [1] George Marsaglia, Wai Wan Tsang.  2000.  The Ziggurat Method for
!>        Generating Random Variables.  'Journal of Statistical Software'.
!>        Vol. 5, Issue 8.  http://www.jstatsoft.org/v05/i08/
!>
!> @author George Marsaglia
!> @author Wai Wan Tsang
!> @author Jason M. Wood
!> @copyright GNU General Public License
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ziggurat
  ! Load intrinsic modules.
  use, intrinsic :: iso_fortran_env

  implicit none

  private

  ! Declare public methods.
  public :: ziggurat_t
  public :: ziggurat_rexp
  public :: ziggurat_rnor
  public :: ziggurat_seed
  public :: ziggurat_shr3
  public :: ziggurat_uni

  ! Mask selecting the low 32 bits of a 64-bit integer.
  integer(kind = int64), parameter :: mask32 = 4294967295_int64

  ! The private state variables for the ziggurat algorithm.
  type :: ziggurat_t
    private
    integer(kind = int32) :: jsr
    integer(kind = int64) :: ke(0:255)
    integer(kind = int32) :: kn(0:127)
    real(kind = real32)   :: fe(0:255)
    real(kind = real32)   :: fn(0:127)
    real(kind = real32)   :: we(0:255)
    real(kind = real32)   :: wn(0:127)
  end type ziggurat_t

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This procedure sets the seed and creates the RNOR and REXP tables for
  !> the ziggurat algorithm.
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @param iii The seed for the ziggurat algorithm. A seed of 123456789 gives
  !>            an all-zero state and must not be used.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ziggurat_seed (state, iii)
    type(ziggurat_t), intent(inout)   :: state
    integer(kind = int32), intent(in) :: iii
    ! Local parameters.
    real(kind = real64), parameter :: m1 = 2147483648.0d0
    real(kind = real64), parameter :: m2 = 4294967296.0d0
    real(kind = real64), parameter :: ve = 3.949659822581572d-3
    real(kind = real64), parameter :: vn = 9.91256303526217d-3
    ! Local variables.
    integer(kind = int32) :: i
    real(kind = real64)   :: de
    real(kind = real64)   :: dn
    real(kind = real64)   :: te
    real(kind = real64)   :: tn
    real(kind = real64)   :: q
    dn = 3.442619855899d0
    de = 7.697117470131487d0
    tn = dn
    te = de
    state%jsr = ieor (123456789_int32, iii)
    ! Tables for RNOR:
    q = vn / exp (-0.5 * dn * dn)
    state%kn(0) = int (dn / q * m1, kind = int32)
    state%kn(1) = 0
    state%wn(0) = real (q / m1, kind = real32)
    state%wn(127) = real (dn / m1, kind = real32)
    state%fn(0) = 1.0
    state%fn(127) = real (exp (-0.5 * dn * dn), kind = real32)
    do i = 126, 1, -1
      dn = sqrt (-2.0 * log (vn / dn + exp (-0.5 * dn * dn)))
      state%kn(i + 1) = int (dn / tn * m1, kind = int32)
      tn = dn
      state%fn(i) = real (exp (-0.5 * dn * dn), kind = real32)
      state%wn(i) = real (dn / m1, kind = real32)
    end do
    ! Tables for REXP:
    q = ve / exp (-de)
    state%ke(0) = int (de / q * m2, kind = int64)
    state%ke(1) = 0
    state%we(0) = real (q / m2, kind = real32)
    state%we(255) = real (de / m2, kind = real32)
    state%fe(0) = 1.0
    state%fe(255) = real (exp (-de), kind = real32)
    do i = 254, 1, -1
      de = -log (ve / de + exp (-de))
      state%ke(i + 1) = int (de / te * m2, kind = int64)
      te = de
      state%fe(i) = real (exp (-de), kind = real32)
      state%we(i) = real (de / m2, kind = real32)
    end do
  end subroutine ziggurat_seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates pseudo-random 32-bit integers.
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return A pseudo-random 32-bit integer.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_shr3 (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    integer(kind = int32)           :: return_value
    ! Local variables.
    integer(kind = int32) :: jz
    integer(kind = int64) :: s
    jz = state%jsr
    state%jsr = ieor (state%jsr, ishft (state%jsr, 13))
    state%jsr = ieor (state%jsr, ishft (state%jsr, -17))
    state%jsr = ieor (state%jsr, ishft (state%jsr, 5))
    ! Add modulo 2**32 in 64-bit arithmetic to avoid signed 32-bit overflow.
    s = iand (int (jz, kind = int64) + int (state%jsr, kind = int64), mask32)
    if (s .gt. int (huge (return_value), kind = int64)) then
      s = s - 4294967296_int64
    end if
    return_value = int (s, kind = int32)
  end function ziggurat_shr3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a uniformly distributed pseudo-random value in the range (0,1).
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return A uniformly distributed pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_uni (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real64)             :: return_value
    ! Local variables.
    real(kind = real64) :: x
    x = real (ziggurat_shr3 (state), kind = real64)
    return_value = 0.5d0 + x * 0.2328306d-9
  end function ziggurat_uni

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a normally distributed pseudo-random value with a mean of 0 and
  !> a variance of 1.
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return A normally distributed pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_rnor (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    ! Local parameters.
    real(kind = real32), parameter :: r = 3.442620  ! Start of the right tail.
    real(kind = real64), parameter :: rinv = 0.2904764d0  ! 0.2904764 is 1/r.
    ! Local variables.
    integer(kind = int32) :: hz
    integer(kind = int32) :: iz
    real(kind = real64) :: u
    real(kind = real32) :: x
    real(kind = real32) :: y
    real(kind = real64) :: z
    hz = ziggurat_shr3 (state)
    iz = iand (hz, 127)
    if (abs (int (hz, kind = int64)) .lt. state%kn(iz)) then
      return_value = hz * state%wn(iz)
    else
      ! RNOR rejection occurred, generate variates from the residue.
      do
        x = hz * state%wn(iz)
        ! Handle the base strip.
        if (iz .eq. 0) then
          do
            x = real (-log (ziggurat_uni (state)) * rinv, kind = real32)
            y = real (-log (ziggurat_uni (state)), kind = real32)
            if (y + y .ge. x * x) exit
          end do
          if (hz .gt. 0) then
            return_value = r + x
          else
            return_value = -r - x
          end if
          return
        end if
        ! Handle the wedges of other strips.
        u = ziggurat_uni (state)
        z = state%fn(iz) + u * (state%fn(iz - 1) - state%fn(iz))
        if (z .lt. exp (-0.5d0 * x * x)) then
          return_value = x
          return
        end if
        ! Try to exit do loop.
        hz = ziggurat_shr3 (state)
        iz = iand (hz, 127)
        if (abs (int (hz, kind = int64)) .lt. state%kn(iz)) then
          return_value = hz * state%wn(iz)
          return
        end if
      end do
    end if
  end function ziggurat_rnor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates exponentially distributed pseudo-random values with a rate (and
  !> mean) of 1. The SHR3 output is interpreted as an unsigned 32-bit integer.
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return An exponential pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_rexp (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    ! Local variables.
    integer(kind = int32) :: iz
    integer(kind = int32) :: jz
    integer(kind = int64) :: ujz
    real(kind = real64)   :: u
    real(kind = real32)   :: x
    real(kind = real64)   :: y
    jz = ziggurat_shr3 (state)
    ujz = iand (int (jz, kind = int64), mask32)
    iz = iand (jz, 255)
    if (ujz .lt. state%ke(iz)) then
      return_value = real (ujz, kind = real32) * state%we(iz)
    else
      ! REXP rejection occurred, generate variates from the residue.
      do
        ! Handles the base strip.
        if (iz .eq. 0) then
          u = ziggurat_uni (state)
          return_value = real (7.69711d0 - log (u), kind = real32)
          return
        end if
        ! Handle the wedges of other strips.
        x = real (ujz, kind = real32) * state%we(iz)
        u = ziggurat_uni (state)
        y = state%fe(iz) + u * (state%fe(iz - 1) - state%fe(iz))
        if (y .lt. exp (-real (x, kind = real64))) then
          return_value = x
          return
        end if
        ! Try to exit do loop.
        jz = ziggurat_shr3 (state)
        ujz = iand (int (jz, kind = int64), mask32)
        iz = iand (jz, 255)
        if (ujz .lt. state%ke(iz)) then
          return_value = real (ujz, kind = real32) * state%we(iz)
          return
        end if
      end do
    end if
  end function ziggurat_rexp

end module ziggurat
