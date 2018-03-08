!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
!
!    Copyright (C) 2000  George Marsaglia and Wai Wan Tsang
!    Copyright (C) 2013-2018  Jason M. Wood, Montana State University
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
!>     use :: ziggurat
!>     implicit none
!>     integer, parameter :: iii = 123456789 !< The RNG seed.
!>     integer            :: i
!>     real               :: x
!>     type(ziggurat_t)   :: rng
!>     ! Seed the RNG.
!>     call ziggurat_seed (rng, iii)
!>     ! Generate some random numbers.
!>     do i = 1, 10
!>       x = ziggurat_uni(rng)
!>       write(*,*) i, x
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

  ! The state variables for the ziggurat algorithm.
  type :: ziggurat_t
    integer(kind = int32) :: hz
    integer(kind = int64) :: iz
    integer(kind = int64) :: jz
    integer(kind = int64) :: jsr
    integer(kind = int64) :: ke(0:255)
    integer(kind = int64) :: kn(0:127)
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
  !> @param iii The seed for the ziggurat algorithm.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ziggurat_seed (state, iii)
    type(ziggurat_t), intent(inout)   :: state
    integer(kind = int32), intent(in) :: iii
    ! Local parameters.
    integer(kind = int64), parameter :: m1 = 2147483648_int64
    integer(kind = int64), parameter :: m2 = 4294967296_int64
    real(kind = real64), parameter   :: ve = 3.949659822581572d-3
    real(kind = real64), parameter   :: vn = 9.91256303526217d-3
    ! Local variables.
    integer(kind = int32) :: i
    real(kind = real64)   :: de
    real(kind = real64)   :: dn
    real(kind = real64)   :: te
    real(kind = real64)   :: tn
    real(kind = real64)   :: q
    dn = 3.442619855899
    de = 7.697117470131487
    tn = dn
    te = de
    state%jsr = ieor (state%jsr, int (iii, kind = int64))
    ! Tables for RNOR:
    q = vn / exp (-0.5 * dn * dn)
    state%kn(0) = int (dn / q * m1, kind = int64)
    state%kn(1) = 0
    state%wn(0) = real (q / m1, kind = real32)
    state%wn(127) = real (dn / m1, kind = real32)
    state%fn(0) = 1.0
    state%fn(127) = real (exp (-0.5 * dn * dn), kind = real32)
    do i = 126, 1, -1
      dn = sqrt (-2.0 * log (vn / dn + exp (-0.5 * dn * dn)))
      state%kn(i + 1) = int (dn / tn * m1, kind = int64)
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
    state%jz = state%jsr
    state%jsr = ieor (state%jsr, ishft (state%jsr, 13))
    state%jsr = ieor (state%jsr, ishft (state%jsr, -17))
    state%jsr = ieor (state%jsr, ishft (state%jsr, 5))
    return_value = int (state%jz + state%jsr, kind = int32)
    return
  end function ziggurat_shr3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a uniformly distributed pseudo-random value in the range [0,1).
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return A uniformly distributed pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_uni (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    return_value = 0.5 + ziggurat_shr3(state) * 0.2328306e-9
    return
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
    real(kind = real32), parameter :: r = 3.442620      ! Start of the right tail.
    real(kind = real32), parameter :: rinv = 0.2904764  ! 0.2904764 is 1/r.
    ! Local variables.
    real(kind = real32) :: x
    real(kind = real32) :: y
    real(kind = real32) :: z
    state%hz = ziggurat_shr3(state)
    state%iz = iand (state%hz, 127)
    if (abs (state%hz) .lt. state%kn(state%iz)) then
      return_value = state%hz * state%wn(state%iz)
    else
      ! RNOR rejection occurred, generate variates from the residue.
      do
        x = state%hz * state%wn(state%iz)
        ! Handle the base strip.
        if (state%iz .eq. 0) then
          do
            x = -log (ziggurat_uni(state)) * rinv
            y = -log (ziggurat_uni(state))
            if (y + y .ge. x * x) exit
          end do
          if (state%hz .gt. 0) then
            return_value = r + x
          else
            return_value = -r - x
          end if
          return
        end if
        ! Handle the wedges of other strips.
        z = state%fn(state%iz) + ziggurat_uni(state) * &
          (state%fn(state%iz - 1) - state%fn(state%iz))
        if (z .lt. exp (-0.5 * x * x)) then
          return_value = x
          return
        end if
        ! Try to exit do loop.
        state%hz = ziggurat_shr3(state)
        state%iz = iand (state%hz, 127)
        if (abs (state%hz) .lt. state%kn(state%iz)) then
          return_value = state%hz * state%wn(state%iz)
          return
        end if
      end do
    end if
    return
  end function ziggurat_rnor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates exponentially distributed pseudo-random values.
  !>
  !> @param state The state variables for the ziggurat algorithm.
  !> @return A exponential pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_rexp (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    ! Local variables.
    real(kind = real64) :: x
    real(kind = real64) :: y
    state%jz = ziggurat_shr3(state)
    state%iz = iand (state%jz, 255)
    if (abs (state%jz) .lt. state%ke(state%iz)) then
      x = state%jz * real (state%we(state%iz), kind = real64)
      return_value = real (x, kind = real32)
    else
      ! REXP rejection occurred, generate variates from the residue.
      do
        ! Handles the base strip.
        if (state%iz .eq. 0) then
          return_value = 7.69711 - log (ziggurat_uni(state))
          exit
        end if
        ! Handle the wedges of other strips.
        x = state%jz * real (state%we(state%iz), kind = real64)
        y = state%fe(state%iz) + ziggurat_uni(state) * &
          (state%fe(state%iz - 1) - state%fe(state%iz))
        if (y .lt. exp (-x)) then
          return_value = real (x, kind = real32)
          exit
        end if
        ! Try to exit do loop.
        state%jz = ziggurat_shr3(state)
        state%iz = iand (state%jz, 255)
        if (state%jz .lt. state%ke(state%iz)) then
          x = state%jz * real (state%we(state%iz), kind = real64)
          return_value = real (x, kind = real32)
          exit
        end if
      end do
    end if
    return
  end function ziggurat_rexp

end module ziggurat
