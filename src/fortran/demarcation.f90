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
  integer                             :: npop
  integer                             :: npoprange(2)
  integer                             :: numincsomega
  integer                             :: numincssigma
  integer                             :: numincsnpop
  integer                             :: numincsxn
  integer, parameter                  :: nparams = 4
  integer, parameter                  :: outputUnit = 2
  double precision                    :: params(nparams)
  double precision                    :: yvalue
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
  call readinput (trim (inputFile), omegarange, sigmarange, npoprange, &
    xnrange, numincsomega, numincssigma, numincsnpop, numincsxn, &
    numcrit, nu, nrep, lengthseq, realdata, crit, whichavg, probthreshold)
  ! Open the output file.
  open(unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  omega = omegarange(1)
  sigma = sigmarange(1)
  xn = xnrange(1)
  istep = numincsnpop
  if (istep .le. 0) istep = 1
  do indexnpop = npoprange(1), npoprange(1), numincsnpop
    npop = indexnpop
    params(1) = log (omega)
    params(2) = log (sigma)
    params(3) = npop
    params(4) = log (xn)
    yvalue = 0.0
    ! yvalue is the negative of the likelihood
    call fredprogram (nparams, params, yvalue)
    write (unit = outputUnit, fmt = *) params(3), ',', yvalue
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
    ! Make sure that npop is in the right range.
    if (params(3) .lt. 1.0d0) params(3) = 1.0d0
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

end program demarcation
