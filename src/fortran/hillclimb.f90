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
! note, this will soon be replaced by a quicker method that Danny
! is writing, which will not require the n-square binning algorithm,
! but rather an nlogn sorting routine
! note, the programs I've changed for nonultra are:
! main, whicheventandwhen, 
! divergenonultra (and deleted binning), binningdanny, coalescence,
! and startpops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program hillclimb
  use methods
  use simplexmethod
  implicit none
  character(len = 256)                :: inputFile
  character(len = 256)                :: outputFile
  logical                             :: fileExists
  integer                             :: i
  integer                             :: ier
  integer                             :: iprint
  integer                             :: iquad
  integer                             :: maxf
  integer                             :: nloop
  integer                             :: npoprange(2)
  integer                             :: numincsomega
  integer                             :: numincssigma
  integer                             :: numincsnpop
  integer                             :: numincsxn
  integer, parameter                  :: nparams = 4
  integer, parameter                  :: outputUnit = 3
  double precision                    :: simp
  double precision                    :: step(nparams)
  double precision                    :: stopcr
  double precision                    :: var(nparams)
  double precision                    :: params(nparams)
  double precision                    :: yvalue
  double precision                    :: omegarange(2)
  double precision                    :: sigmarange(2)
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
  inputFile = 'hClimbIn.dat'
  outputFile = 'hClimbOut.dat'
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
        write (unit = *, fmt = *) "Expected: hClimbIn.dat hClimbOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The hclimbIn.dat file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read in the input file.
  call readinput (trim (inputFile), omegarange, sigmarange, npoprange, &
    xnrange, numincsomega, numincssigma, numincsnpop, numincsxn, &
    numcrit, nu, nrep, lengthseq, realdata, crit, whichavg, probthreshold)
  ! Open the output file.
  open (unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  ! Set max. no. of function evaluations = MAXF, print every IPRINT.
  maxf = 100
  iprint = 1
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
  params(1) = log (omegaRange(1))
  params(2) = log (sigmaRange(1))
  params(3) = npopRange(1)
  params(4) = log (xnRange(1))
  step(1) = log (omegaRange(1)) / 2.0d0
  step(2) = log (sigmaRange(1)) / 2.0d0
  step(3) = npopRange(1) / 2.0d0
  step(4) = log (xnRange(1)) / 2.0d0
  yvalue = 0.0
  ! Danny, yvalue is just a dummy variable that gets passed to my
  ! program in the CALL statement, but the meaningful passing of 
  ! variables occurs through the COMMON block
  ! yvalue is the value returned for success level determined by jwhichxavg
  call nelmead (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
    outputUnit)
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
    ! Make sure that npop is in the right range.
    if (params(3) .lt. 2) params(3) = 2
    if (params(3) .gt. nu) params(3) = nu
    omega = exp (params(1))
    sigma = exp (params(2))
    npop = nint (params(3))
    xn = exp (params(4))
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

end program hillclimb
