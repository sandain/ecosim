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
! This is a hillclimbing version of finding the npop confidence interval.  For each value of npop tested,
! xn stays fixed at infinity, and the optimal values of omega and sigma are found for that npop.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program npopCI
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
  integer                             :: indexnpop
  integer                             :: npop
  integer                             :: npoprange(2)
  integer                             :: numincsomega
  integer                             :: numincssigma
  integer                             :: numincsnpop
  integer                             :: numincsxn
  integer, parameter                  :: outputUnit = 2
  integer, parameter                  :: nparams = 2
  double precision                    :: simp
  double precision                    :: stopcr
  double precision                    :: yvalue
  double precision                    :: params(nparams)
  double precision                    :: step(nparams)
  double precision                    :: var(nparams)
  double precision                    :: omega
  double precision                    :: omegarange(2)
  double precision                    :: sigma  
  double precision                    :: sigmarange(2)
  double precision                    :: xn
  double precision                    :: xnrange(2)
  real                                :: probthreshold
  ! bldanny common block
  integer :: numcrit
  integer :: nu
  integer :: nrep
  integer :: lengthseq
  integer :: realdata(1000)
  integer :: whichavg
  real    :: crit(1000)
  common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,whichavg
  ! Nelder Mead function.
  procedure(nelmeadFunction), pointer :: functn
  functn => fredprogram
  ! Provide default file names to use.
  inputFile = 'npopIn.dat'
  outputFile = 'npopOut.dat'
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
        write (unit = *, fmt = *) "Expected: npopIn.dat npopOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The npopIn.dat file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read in the input file.
  call readinput (trim (inputFile), omegarange, sigmarange, npoprange, &
    xnrange, numincsomega, numincssigma, numincsnpop, numincsxn, &
    numcrit, nu, nrep, lengthseq, realdata, crit, whichavg, probthreshold)
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
  ! Fit a quadratic surface to be sure a minimum has been found.
  iquad = 0
  ! As function value is being evaluated in DOUBLE PRECISION, it
  ! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
  ! we should get about 9 dec. digits accuracy in fitting the surface.
  simp = 1.0d-6
  xn = xnRange(1)
  omega = omegaRange(1)
  sigma = sigmaRange(1)
  do indexnpop = npopRange(1), npopRange(2), numincsnpop
    npop = indexnpop
    ! Note, for npop values besides the original one, 
    ! we will start with the omega and sigma values
    ! calculated for the previous npop value. 
    params(1) = log (omega)
    params(2) = log (sigma)
    step(1) = log (omega) / 2.0d0
    step(2) = log (sigma) / 2.0d0
    if (log (omega) .lt. 0.3d0 .and. log (omega) .gt. -0.3d0) step(1) = 0.15d0
    if (log (sigma) .lt. 0.3d0 .and. log (sigma) .gt. -0.3d0) step(2) = 0.15d0
    yvalue = 0.0d0
    call nelmead (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
     outputUnit, probthreshold)
    write (unit = outputUnit, fmt = *) omega, sigma, npop, xn, yvalue
    if (-yvalue .lt. probthreshold) exit
  end do
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
    integer          :: npop
    double precision :: omega
    double precision :: sigma
    double precision :: xn
    real             :: avgsuccess(6)
    ! bldanny common block
    integer :: numcrit
    integer :: nu
    integer :: nrep
    integer :: lengthseq
    integer :: realdata(1000)
    integer :: whichavg
    real    :: crit(1000)
    common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,whichavg
    omega = exp (params(1))
    sigma = exp (params(2))
    if (omega .lt. 1.0d-7) omega = 1.0d-7
    if (sigma .lt. 1.0d-7) sigma = 1.0d-7
    if (debug) then
      write (unit = *, fmt = *) 'omega= ', omega
      write (unit = *, fmt = *) 'sigma= ', sigma
      write (unit = *, fmt = *) 'npop= ', npop
      write (unit = *, fmt = *) 'xn= ', xn
    end if
    ! avgsuccess is a count of the number of results that are within X% tolerance
    ! for a particular set of parameter values.
    call runProgram (omega, sigma, npop, xn, numcrit, nu, nrep, &
      lengthseq, realdata, crit, avgsuccess)
    yvalue = -1.0d0 * avgsuccess(whichavg)
    if (debug) then
      write (unit = *, fmt = *) 'yvalue= ', yvalue
      write (unit = *, fmt = *)
    end if
    return
  end subroutine fredprogram

end program npopCI
