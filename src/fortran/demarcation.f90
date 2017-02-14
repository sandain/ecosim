!    Ecotype Simulation models the sequence diversity within a bacterial clade as
!    the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
! 
!    Copyright (C) 2009  Fred Cohan, Wesleyan University
!                        Carlo Francisco, Wesleyan University
!                        Danny Krizanc, Wesleyan University
!                        Andrew Warner, Wesleyan University
!                        Jason Wood, Montana State University
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hill-climbing version of non-ultrametric design
! This is a non-ultrametric divergence design, with drift but no
! recombination. September 21, 2003
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program demarcation
  use methods
  implicit none
  ! Local variables
  character(len = 256)                :: inputFile
  character(len = 256)                :: outputFile
  logical                             :: fileExists
  integer                             :: i
  integer                             :: istep
  integer                             :: indexnpop
  integer, parameter                  :: nparams = 4
  integer, parameter                  :: outputUnit = 2
  double precision                    :: params(nparams)
  double precision                    :: yvalue
  type(acinas_data)                   :: acinas
  type(number_increments)             :: numincs
  type(parameters_data)               :: bottom
  type(parameters_data)               :: top
  type(parameters_data)               :: parameters
  ! Provide default file names to use.
  inputFile = 'demarcationIn.dat'
  outputFile = 'demarcationOut.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, inputFile)
      case (2)
        call getArgument (2, outputFile)
      case (3)
        call getArgument (3, numberThreads)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) "Expected: demarcationIn.dat demarcationOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The demarcationIn.dat file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read in the input file.
  call readinput (trim (inputFile), acinas, bottom, top, numincs)
  ! Open the output file.
  open(unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  parameters%omega = bottom%omega
  parameters%sigma = bottom%sigma
  parameters%xn = bottom%xn
  istep = numincs%npop
  if (istep .le. 0) istep = 1
  do indexnpop = bottom%npop, top%npop, numincs%npop
    parameters%npop = indexnpop
    params(1) = log (parameters%omega)
    params(2) = log (parameters%sigma)
    params(3) = parameters%npop
    params(4) = log (parameters%xn)
    yvalue = 0.0
    ! yvalue is the negative of the likelihood
    call fredprogram (nparams, params, yvalue)
    write (unit = outputUnit, fmt = *) params(3), ',', yvalue
  end do
  ! Close the random number generator.
  call randomClose ()
  ! Close the output file.
  close (unit = outputUnit)
  stop

  contains

  subroutine fredprogram (nparams, params, yvalue)
    ! Subroutine parameters
    integer, intent(in)             :: nparams
    double precision, intent(inout) :: params(nparams)
    double precision, intent(out)   :: yvalue
    ! Local variables
    real             :: avgsuccess(6)
    ! Make sure that npop is in the right range.
    if (params(3) .lt. 1.0d0) params(3) = 1.0d0
    if (params(3) .gt. acinas%nu) params(3) = acinas%nu
    parameters%omega = exp (params(1))
    parameters%sigma = exp (params(2))
    parameters%npop = nint (params(3))
    parameters%xn = exp (params(4))
    if (debug) then
      write (unit = *, fmt = *) 'omega= ', parameters%omega
      write (unit = *, fmt = *) 'sigma= ', parameters%sigma
      write (unit = *, fmt = *) 'npop= ', parameters%npop
      write (unit = *, fmt = *) 'xn= ', parameters%xn
    end if
    ! avgsuccess is a count of the number of results that are within X% tolerance
    ! for a particular set of parameter values.
    call simulation (acinas, parameters, avgsuccess)
    yvalue = -1.0d0 * avgsuccess(acinas%whichavg)
    if (debug) then
      write (unit = *, fmt = *) 'yvalue= ', yvalue
      write (unit = *, fmt = *)
    end if
    return
  end subroutine fredprogram

end program demarcation
