!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation and periodic
!    selection, yielding a certain number of ecotypes.
!
!    Copyright (C) 2009-2013  Fred Cohan, Wesleyan University
!                             Danny Krizanc, Wesleyan University
!                             Jason M. Wood, Montana State University
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
!> The hillclimb program uses the Nelder-Mead simplex method for function
!> minimization to optimize the solution provided by the brute force search
!> method.
!>
!> @pre  Requires that the hclimbIn.dat file exists and that it contains all
!>         of the variables necessary to perform the hillclimbing.
!>
!> @post Generates the hillclimbOut.dat file that contains the optimized
!>         solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program hillclimb
  use methods
  use simplexmethod
  implicit none
  ! Local variables
  character(len = 256) :: inputFile
  character(len = 256) :: outputFile
  logical              :: fileExists
  integer              :: npop
  integer              :: i
  integer              :: ier
  integer              :: iprint
  integer              :: iquad
  integer              :: lout
  integer              :: maxf
  integer              :: nloop
  integer, parameter   :: nparams = 3
  integer, parameter   :: outputUnit = 4
  double precision     :: omega
  double precision     :: sigma
  double precision     :: simp
  double precision     :: step(nparams)
  double precision     :: stopcr
  double precision     :: var(nparams)
  double precision     :: params(nparams)
  double precision     :: yvalue
  ! bldanny common block
  integer :: numcrit
  integer :: nu
  integer :: nrep
  integer :: lengthseq
  integer :: realdata(1000)
  integer :: jwhichxavg
  real    :: crit(1000)
  common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
  ! The function to be used by the Nelder-Mead minimization function.
  procedure(nelmeadFunction), pointer :: functn
  functn => callfredprogram
  ! Provide default file names to use.
  inputFile = 'hillclimbIn.dat'
  outputFile = 'hillclimbOut.dat'
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
        write (unit = *, fmt = *) "Expected: hillclimbIn.dat hillclimbOut.dat"
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
  ! Read the input file.
  call readinput (trim (inputFile), omega, sigma, npop)
  ! Open the output file.
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! Set max. no. of function evaluations = maxf, print every iprint.
  maxf = 100
  if (debug) then
    iprint = 1
  else
    iprint = -1
  end if
  ! Send output to stdout (usually unit 6)
  lout = 6
  ! Set value for stopping criterion.  Stopping occurs when the
  ! standard deviation of the values of the objective function at
  ! the points of the current simplex < stopcr.
  stopcr = 1.0d-1
  nloop = 8
  ! Fit a quadratic surface to be sure a minimum has been found.
  iquad = 0
  ! As function value is being evaluated in double precision, it
  ! should be accurate to about 15 decimals.  If we set simp = 1.0d-6,
  ! we should get about 9 dec. digits accuracy in fitting the surface.
  simp = 1.0d-6
  ! Return value starts off at zero.
  yvalue = 0.0
  ! Make sure omega and sigma are greater than zero.
  if (omega .gt. 0.0 .and. sigma .gt. 0.0) then
    ! Setup the parameters for Nelder-Mead.
    params(1) = log (omega)
    step(1) = log (omega) / 2.0
    params(2) = log (sigma)
    step(2) = log (sigma) / 2.0
    params(3) = npop
    step(3) = npop / 2.0
    call nelmead (params, step, nparams, yvalue, maxf, iprint, stopcr, &
      nloop, iquad, simp, var, functn, ier, lout)
    omega = exp (params(1))
    sigma = exp (params(2))
    npop = nint (params(3))
  end if
  ! Output the answer.
  write (unit = outputUnit, fmt = *) omega, sigma, npop, yvalue
  ! Close the random number generator.
  call randomClose ()
  ! Close the output file.
  close (unit = outputUnit)
  ! Successful termination of program.
  stop

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This subroutine is called by the Nelder-Mead simplex method, using the
  !> runprogram subroutine to calculate the yvalue.
  !>
  !> @param[in]     nparams       The number of parameters.
  !> @param[in,out] params        The parameters (omega, sigma, npop).
  !> @param[out]    yvalue        The resulting value of the function.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine callfredprogram (nparams, params, yvalue)
    integer, intent(in)             :: nparams
    double precision, intent(inout) :: params(nparams)
    double precision, intent(out)   :: yvalue
    ! Local variables
    double precision   :: avgsuccess(6)
    double precision   :: omega
    double precision   :: sigma
    integer            :: npop
    ! bldanny common blocknparams
    integer :: numcrit
    integer :: nu
    integer :: nrep
    integer :: lengthseq
    integer :: realdata(1000)
    integer :: jwhichxavg
    real    :: crit(1000)
    common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
    if (params(3) .lt. 2.0) params(3) = 2
    if (params(3) .gt. nu)  params(3) = nu
    omega = exp (params(1))
    sigma = exp (params(2))
    npop = nint (params(3))
    if (debug) then
      write (unit = *, fmt = *) 'omega= ', omega
      write (unit = *, fmt = *) 'sigma= ', sigma
      write (unit = *, fmt = *) 'npop= ', npop
    end if
    call runprogram (omega, sigma, npop, numcrit, nu, nrep, &
      lengthseq, realdata, crit, avgsuccess)
    yvalue = -1.0 * avgsuccess(jwhichxavg)
    if (debug) then
      write (unit = *, fmt = *) 'yvalue= ', yvalue
    end if
    return
  end subroutine callfredprogram

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read the variables contained in the input file.
  !>
  !> @param[in]     fname         The path and file name of the input file.
  !> @param[out]    omega         The omega value to be tested.
  !> @param[out]    sigma         The sigma value to be tested.
  !> @param[out]    npop          The npop value to be tested.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readinput (fname, omega, sigma, npop)
    character(len = *), intent(in) :: fname
    integer, intent(out)           :: npop
    double precision, intent(out)  :: omega
    double precision, intent(out)  :: sigma
    ! Local variables
    integer            :: iii
    integer            :: jcrit
    integer, parameter :: input_unit = 1
    ! bldanny common block
    integer :: numcrit
    integer :: nu
    integer :: nrep
    integer :: lengthseq
    integer :: realdata(1000)
    integer :: jwhichxavg
    real    :: crit(1000)
    common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
    ! Open the input file.
    open (unit = input_unit, file = fname, action = 'read', &
      access = 'sequential', form = 'formatted')
    ! numcrit is the number of criteria for making cluster bins
    read (unit = input_unit, fmt = *) numcrit
    do jcrit = 1, numcrit
      ! realdata(jcrit) is the number of bins in the real data at criterion
      ! jcrit
      read (unit = input_unit, fmt = *) realdata(jcrit)
    end do
    do jcrit = 1, numcrit
      ! crit(jcrit) is the jcrit th criterion value
      read (unit = input_unit, fmt = *) crit(jcrit)
    end do
    ! omega is the rate of niche invasion per eligible parental population
    ! measured as niche invasions per nucleotide substitution in a given gene
    read (unit = input_unit, fmt = *) omega
    ! sigma is the rate of periodic selection per eligible population,
    ! measured as periodic selection events per population per nucleotide
    ! substitution in a given gene
    read (unit = input_unit, fmt = *) sigma
    ! npop is the number of ecotypes assumed to be in the environmental DNA
    ! sample
    read (unit = input_unit, fmt = *) npop
    ! nu is the number of homologous gene sequences in the environmental
    ! sample following Acinas et al., this should be in the thousands.
    read (unit = input_unit, fmt = *) nu
    ! nrep is the number of replicate simulations for a given set of sigma,
    ! omega, and npop
    read (unit = input_unit, fmt = *) nrep
    ! iii is the odd random number seed (up to nine digits)
    read (unit = input_unit, fmt = *) iii
    ! Initialize the random number generator.
    call randomInitialize (iii)
    ! lengthseq is the length in nucleotides of the sequence analyzed
    read (unit = input_unit, fmt = *) lengthseq
    ! jwhichxavg gives the precision:
    ! 1=5x, 2=2x, 3=1.5x, 4=1.25x, 5=1.1x, 6=1.05x
    read (unit = input_unit, fmt = *) jwhichxavg
    ! The highest sequence identity criterion cannot be 1.0 but should be
    ! 1-(1/(2*lengthseq))
    crit(numcrit) = 1.0 - 1.0 / (2.0 * lengthseq)
    ! Close the input file.
    close (unit = input_unit)
    return
  end subroutine readinput

end program hillclimb
