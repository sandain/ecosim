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
!> The bruteforce program uses the runprogram subroutine to search the range
!> of values spanned by omegabot to omegatop, sigmabot to sigmatop, and
!> npopbot to npoptop for results that fall within 500%, 200%, 150%, 125%,
!> 110%, and 105% of the expected values.
!>
!> @pre  Requires that the bruteforceIn.dat file exist and that it contains
!>         all of the variables necessary to perform the brute force search.
!>
!> @post Generates the bruteforceOut.dat file that contains the output from
!>         the run program subroutine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program bruteforce
  use methods
  implicit none
  ! Local variables
  character(len = 256)          :: inputFile
  character(len = 256)          :: outputFile
  character(len = *), parameter :: outputFormat = &
    "(f15.7,',',f15.7,',',i10,6(',',1e18.10))"
  logical                       :: fileExists
  integer                       :: npop
  integer                       :: npopbot
  integer                       :: npoptop
  integer                       :: i
  integer                       :: indexsigma
  integer                       :: indexomega
  integer                       :: indexnpop
  integer                       :: numincsomega
  integer                       :: numincssigma
  integer                       :: numincsnpop
  integer, parameter            :: outputUnit = 2
  double precision              :: omega
  double precision              :: omegabot
  double precision              :: omegatop
  double precision              :: sigma
  double precision              :: sigmabot
  double precision              :: sigmatop
  double precision              :: avgsuccess(6)
  double precision              :: diffomega
  double precision              :: diffnpop
  double precision              :: diffsigma
  ! bldanny common block
  integer :: numcrit
  integer :: nu
  integer :: nrep
  integer :: lengthseq
  integer :: realdata(1000)
  integer :: jwhichxavg
  real    :: crit(1000)
  common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
  ! Provide default file names to use.
  inputFile = 'bruteforceIn.dat'
  outputFile = 'bruteforceOut.dat'
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
        write (unit = *, fmt = *) &
          "Expected: bruteforceIn.dat bruteforceOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) &
      "The bruteforceIn.dat file was not found at: ", trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read the input file.
  call readinput (trim (inputFile), omegabot, omegatop, sigmabot, &
    sigmatop, npopbot, npoptop, numincsomega, numincssigma, numincsnpop)
  ! Open the output file.
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! For brute force main program, we will run subroutine "runprogram"
  ! over the range of values spanned by omegabot to omegatop, sigmabot
  ! to sigmatop, and npopbot to npoptop
  diffomega = log10 (omegatop) - log10 (omegabot)
  diffsigma = log10 (sigmatop) - log10 (sigmabot)
  diffnpop = log10 (real (npoptop)) - log10 (real (npopbot))
  if (numincsomega .eq. 0) then
    omegatop = omegabot + 1
    diffomega = log10 (omegatop) - log10 (omegabot)
    diffomega = 2.0 * diffomega
    numincsomega = 1
  endif
  if (numincssigma .eq. 0) then
    sigmatop = sigmabot + 1
    diffsigma = log10 (sigmatop) - log10 (sigmabot)
    diffsigma = 2.0 * diffsigma
    numincssigma = 1
  endif
  if (numincsnpop .eq. 0) then
    npoptop = npopbot + 1
    diffnpop = log10 (real (npoptop)) - log10 (real (npopbot))
    diffnpop = 2.0 * diffnpop
    numincsnpop = 1
  endif
  indexomega = 0
  do while (indexomega .lt. numincsomega)
    omega = 10.0 ** (log10 (omegabot) + indexomega * &
      (diffomega / numincsomega) * 0.999)
    indexsigma = 0
    do while (indexsigma .lt. numincssigma)
      sigma = 10.0 ** (log10 (sigmabot) + indexsigma * &
        (diffsigma / numincssigma) * 0.999)
      indexnpop = 0
      do while (indexnpop .lt. numincsnpop)
        npop = nint (10.0 ** (log10 (real (npopbot)) + indexnpop * &
          (diffnpop / numincsnpop) * 0.999))
        if (npop .gt. nu) npop = nu
        call runprogram (omega, sigma, npop, numcrit, nu, nrep, &
          lengthseq, realdata, crit, avgsuccess)
        ! Output the result if the 500% average success level is greater
        ! than zero.
        if (avgsuccess(1) .gt. 0.0) then
          write (unit = outputUnit, fmt = outputFormat) omega, sigma, npop, &
            (avgsuccess(i), i = 1, 6)
          if (debug) then
            write (unit = *, fmt = outputFormat) omega, sigma, npop, &
              (avgsuccess(i), i = 1, 6)
          end if
        end if
        indexnpop = indexnpop + 1
      end do
      indexsigma = indexsigma + 1
    end do
    indexomega = indexomega + 1
  end do
  ! Close the random number generator.
  call randomClose ()
  ! Close the output file.
  close (unit = outputUnit)
  ! Successful termination of program.
  stop

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read the variables contained in the input file.
  !>
  !> @param[in]     fname         The path and file name of the input file.
  !> @param[out]    omegabot      The bottom of the omega values to be tested.
  !> @param[out]    omegatop      The top of the omega values to be tested.
  !> @param[out]    sigmabot      The bottom of the sigma values to be tested.
  !> @param[out]    sigmatop      The top of the sigma values to be tested.
  !> @param[out]    npopbot       The bottom of the npop values to be tested.
  !> @param[out]    npoptop       The top of the npop values to be tested.
  !> @param[out]    numincsomega  The number of increments to be used for
  !>                                omega values.
  !> @param[out]    numincssigma  The number of increments to be used for
  !>                                sigma values.
  !> @param[out]    numincsnpop   The number of increments to be used for
  !>                                npop values.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readinput(fname, omegabot, omegatop, sigmabot, sigmatop, &
    npopbot, npoptop, numincsomega, numincssigma, numincsnpop)
    character(len = *), intent(in) :: fname
    integer, intent(out)           :: npopbot
    integer, intent(out)           :: npoptop
    double precision, intent(out)  :: omegabot
    double precision, intent(out)  :: omegatop
    double precision, intent(out)  :: sigmabot
    double precision, intent(out)  :: sigmatop
    integer, intent(out)           :: numincsomega
    integer, intent(out)           :: numincssigma
    integer, intent(out)           :: numincsnpop
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
      ! realdata(jcrit) is the number of bins in the real data at
      ! criterion jcrit
      read (unit = input_unit, fmt = *) realdata(jcrit)
    end do
    do jcrit = 1, numcrit
      ! crit(jcrit) is the jcrit th criterion value
      read (unit = input_unit, fmt = *) crit(jcrit)
    end do
    ! omega is the rate of niche invasion per eligible parental population
    ! measured as niche invasions per nucleotide substitution in a given gene
    ! omegabot and omegatop are the range
    read (unit = input_unit, fmt = *) omegabot, omegatop
    ! sigma is the rate of periodic selection per eligible population,
    ! measured as periodic selection events per population per nucleotide
    ! substitution in a given gene
    read (unit = input_unit, fmt = *) sigmabot, sigmatop
    ! npop is the number of ecotypes assumed to be in the environmental DNA
    ! sample
    read (unit = input_unit, fmt = *) npopbot, npoptop
    ! numincs values are the numbers of increments to be investigated
    read (unit = input_unit, fmt = *) numincsomega, numincssigma, &
      numincsnpop
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
    ! The highest sequence identity criterion cannot be 1.0 but should be
    ! 1-(1/(2*lengthseq))
    crit(numcrit) = 1.0 - 1.0 / (2.0 * lengthseq)
    ! Close the input file.
    close (unit = input_unit)
    return
  end subroutine readinput

end program bruteforce
