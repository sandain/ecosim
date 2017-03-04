!    Ecotype Simulation models the sequence diversity within a bacterial clade as
!    the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
! 
!    Copyright (C) 2009  Fred Cohan, Wesleyan University
!                        Danny Krizanc, Wesleyan University
!                        Jason M. Wood, Montana State University
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a hillclimbing version of finding the drift confidence interval.  For each value of npop tested,
! xn stays fixed at infinity, and the optimal values of omega and sigma are found for that npop.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program driftCI
  use methods
  use simplexmethod
  implicit none
  ! Local variables
  character(len = 256)                :: inputFile
  character(len = 256)                :: outputFile
  logical                             :: fileExists
  integer                             :: i
  integer                             :: ier
  integer                             :: iprint
  integer                             :: iquad
  integer                             :: maxf
  integer                             :: nloop
  integer                             :: indexxn
  integer                             :: xincs
  integer, parameter                  :: nparams = 3
  integer, parameter                  :: outputUnit = 2
  double precision                    :: simp
  double precision                    :: stopcr
  double precision                    :: yvalue
  double precision                    :: params(nparams)
  double precision                    :: step(nparams)
  double precision                    :: var(nparams)
  double precision                    :: diffxn
  type(acinas_data)                   :: acinas
  type(number_increments)             :: numincs
  type(parameters_data)               :: bottom
  type(parameters_data)               :: top
  type(parameters_data)               :: parameters
  ! Nelder Mead function.
  procedure(nelmeadFunction), pointer :: functn
  functn => fredprogram
  ! Provide default file names to use.
  inputFile = 'driftIn.dat'
  outputFile = 'driftOut.dat'
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
        write (unit = *, fmt = *) "Expected: driftIn.dat driftOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The driftIn.dat file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read in the input file.
  call readinput (trim (inputFile), acinas, bottom, top, numincs)
  ! Open the output file.
  open(unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  ! Set max. no. of function evaluations = MAXF, print every IPRINT.
  maxf = 100
  iprint = -1
  ! Set value for stopping criterion.   Stopping occurs when the
  ! standard deviation of the values of the objective function at
  ! the points of the current simplex < stopcr.
  stopcr = 1.0d-1
  nloop = 8
  !  Fit a quadratic surface to be sure a minimum has been found.
  iquad = 0
  ! As function value is being evaluated in DOUBLE PRECISION, it
  ! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
  ! we should get about 9 dec. digits accuracy in fitting the surface.
  simp = 1.0d-6
  ! Fixed factor for drift:
  parameters%xn = bottom%xn
  ! Starting values for omega, sigma, and npop:
  parameters%omega = bottom%omega
  parameters%sigma = bottom%sigma
  parameters%npop = bottom%npop
  ! if xincs=0, it means that only one set of omega values will be used
  xincs = 1
  diffxn = log (top%xn) - log (bottom%xn)
  if (numincs%xn .eq. 0) then
    top%xn = bottom%xn * 10.0d0
    diffxn = log (top%xn) - log (bottom%xn)
    numincs%xn = 1
    xincs = 0
  endif
  do indexxn = 0, numincs%xn
    parameters%xn = exp (log (bottom%xn) + indexxn * (diffxn / numincs%xn) * 0.999d0)
    params(1) = log (parameters%omega)
    params(2) = log (parameters%sigma)
    params(3) = parameters%npop
    step(1) = log (parameters%omega) / 2.0d0
    step(2) = log (parameters%sigma) / 2.0d0
    step(3) = parameters%npop / 2.0
    if (log (parameters%omega) .lt. 0.30d0 .and. log (parameters%omega) .gt. -0.30d0) step(1) = 0.15d0
    if (log (parameters%sigma) .lt. 0.30d0 .and. log (parameters%sigma) .gt. -0.30d0) step(2) = 0.15d0
    yvalue = 0.0
    ! yvalue is the negative of the likelihood
    call nelmead (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
     outputUnit, acinas%probthreshold)
    write (unit = outputUnit, fmt = *) parameters%omega, parameters%sigma, parameters%npop, parameters%xn, yvalue
    if (-yvalue .lt. acinas%probthreshold .or. xincs .eq. 0) exit
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
    parameters%omega = exp (params(1))
    parameters%sigma = exp (params(2))
    parameters%npop = nint (params(3))
    if (parameters%omega .lt. 1.0d-7) parameters%omega = 1.0d-7
    if (parameters%sigma .lt. 1.0d-7) parameters%sigma = 1.0d-7
    if (parameters%npop .lt. 1) parameters%npop = 1
    if (parameters%npop .gt. acinas%nu) parameters%npop = acinas%nu
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

end program driftCI
