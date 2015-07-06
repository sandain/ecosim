!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation and periodic
!    selection, yielding a certain number of ecotypes.
!
!    Copyright (C) 2009-2014  Fred Cohan, Wesleyan University
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
!> The methods module contains all of the common subroutines used in the
!> various programs of the Ecotype Simulation application.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module methods
  use ISO_FORTRAN_ENV
  use darray
  use ziggurat
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  ! Declare public methods.
  public :: getArgument
  public :: gtest
  public :: randomClose
  public :: randomInitialize
  public :: randomNumber
  public :: runFredProgram

  ! Declare public variables.
  logical, public :: debug = .false.       !< Display debug information.
  integer(kind = int32), public :: &
    numberThreads = 1                      !< The number of threads to start.

  ! Declare private global parameters.
  integer(kind = int32), parameter :: &
    EVENT_NICHE_INVASION     = 1001        !< Niche invasion event.
  integer(kind = int32), parameter :: &
    EVENT_PERIODIC_SELECTION = 1002        !< Periodic selection event.

  type(ziggurat_t), allocatable :: rng(:)  !< The state variables for the RNG.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets an argument from the command line.
  !>
  !> @param[in]     arg           The argument number to get.
  !> @param[out]    out           The value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface getArgument
    module procedure getLogicalArgument
    module procedure getIntegerArgument
    module procedure getRealArgument
    module procedure getStringArgument
  end interface getArgument

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Binning
  !>
  !> @param[out]    bin           The number of clusters at the identity
  !>                                level sought.
  !> @param[in]     div           The number of substitutions per site for
  !>                                each ancestor.
  !> @param[in]     crit          The list of bin identity levels (between
  !>                                0.0 and 1.0).
  !> @param[in]     ncoalesce     The strains that will be coalesced.
  !> @param[in]     numcrit       The number of criteria for making cluster
  !>                                bins.
  !> @param[in]     numanctot     The total number of ancestors, over all
  !>                                strains.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine binning (bin, div, crit, ncoalesce, numcrit, numanctot)
    integer(kind = int32), intent(out)   :: bin(:)
    integer(kind = int32), intent(in)    :: numcrit
    integer(kind = int32), intent(in)    :: numanctot
    real(kind = real32), intent(in)      :: crit(:)
    type(darrayInteger), intent(in)      :: ncoalesce
    type(darrayReal), intent(in)         :: div
    ! Local variables.
    integer(kind = int32) :: lused
    integer(kind = int32) :: janc
    integer(kind = int32) :: jcrit
    integer(kind = int32) :: nbins
    real(kind = real32)   :: criterion
    real(kind = real32)   :: div2
    nbins = 1
    lused = 0
    do jcrit = 1, numcrit
      criterion = crit(jcrit)
      janc = numanctot - lused
      do while (janc .ge. 1)
        div2 = -1.5 * (exp ((-4.0 / 3.0) * darrayGet (div, janc)) - 1.0)
        if (div2 .lt. 1.0 - criterion) then
          bin(jcrit) = nbins
          exit
        end if
        nbins = nbins + darrayGet (ncoalesce, janc) - 1
        lused = lused + 1
        janc = janc - 1
      end do
    end do
    return
  end subroutine binning

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Canonical
  !>
  !> @param[in]     npop          The number of ecotypes assumed to be in the
  !>                                environmental DNA sample.
  !> @param[in]     nu            The number of homologous gene sequences in
  !>                                the environmental sample.
  !> @param[in,out] numstrain     An array of the number of strains in each
  !>                                ecotype.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine canonical (npop, nu, numstrain)
    integer(kind = int32), intent(in)    :: npop
    integer(kind = int32), intent(in)    :: nu
    integer(kind = int32), intent(out)   :: numstrain(:)
    ! Local variables.
    integer(kind = int32) :: upto
    integer(kind = int32) :: ipop
    integer(kind = int32) :: istrain
    integer(kind = int32) :: numberofstrains
    real(kind = real32)   :: x
    numstrain(1) = nu
    if (npop .gt. 1) then
      do upto = 2, npop
        ! for the upto - 1 population, choose a random population (ipop)
        do
          call randomNumber (x)
          ipop = int (x * (upto - 1), kind = int32) + 1
          numberofstrains = numstrain(ipop)
          if (numberofstrains .gt. 1) exit
        end do
        ! next, pick a random breakpoint
        call randomNumber (x)
        istrain = int (x * (numberofstrains - 1), kind = int32) + 1
        numstrain(ipop) = istrain
        numstrain(upto) = numberofstrains - istrain
      end do
    end if
    return
  end subroutine canonical

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs the niche invasion for the nascent population.
  !>
  !> @param[in,out] activepop     The number of active populations.
  !> @param[in,out] numanctot     The total number of ancestors, over all
  !>                                strains.
  !> @param[in,out] numstrain     An array of the number of strains in each
  !>                                ecotype.
  !> @param[in,out] ncoalesce     The strains that will be coalesced.
  !> @param[in]     lengthseq     The length of the sequences.

  !> @param[in]     time          The total time (in substitutions) before
  !>                                the present for a coalescence event.
  !> @param[in,out] div           The number of substitutions per site for
  !>                                each ancestor.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doNicheInvasion (activepop, numanctot, numstrain, ncoalesce, &
    lengthseq, time, div)
    integer(kind = int32), intent(inout) :: activepop
    integer(kind = int32), intent(inout) :: numanctot
    integer(kind = int32), intent(inout) :: numstrain(:)
    integer(kind = int32), intent(in)    :: lengthseq
    real(kind = real32), intent(in)      :: time
    type(darrayInteger), intent(inout)   :: ncoalesce
    type(darrayReal), intent(inout)      :: div
    ! Local variables.
    integer(kind = int32) :: popfornascent
    real(kind = real32)   :: x
    if (activepop .eq. 0) return
    ! Choose the population for the nascent population
    call randomNumber (x)
    popfornascent = int (x * activepop, kind = int32) + 1
    numanctot = numanctot + 1
    call darraySet (ncoalesce, numanctot, numstrain(popfornascent) + 1)
    call darraySet (div, numanctot, time / lengthseq)
    ! next, delete the nascent population
    numstrain(popfornascent) = numstrain(activepop)
    activepop = activepop - 1
    return
  end subroutine doNicheInvasion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs the periodic selection for the given population.
  !>
  !> @param[in]     activepop     The number of active populations.
  !> @param[in,out] numanctot     The total number of ancestors, over all
  !>                                strains.
  !> @param[in,out] numstrain     An array of the number of strains in each
  !>                                ecotype.
  !> @param[in,out] ncoalesce     The strains that will be coalesced.
  !> @param[in]     lengthseq     The length of the sequences.
  !> @param[in]     time          The total time (in substitutions) before
  !>                                the present for a coalescence event.
  !> @param[in,out] div           The number of substitutions per site for
  !>                                each ancestor.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doPeriodicSelection (activepop, numanctot, numstrain, &
    ncoalesce, lengthseq, time, div)
    integer(kind = int32), intent(in)    :: activepop
    integer(kind = int32), intent(inout) :: numstrain(:)
    integer(kind = int32), intent(inout) :: numanctot
    integer(kind = int32), intent(in)    :: lengthseq
    real(kind = real32), intent(in)      :: time
    type(darrayInteger), intent(inout)   :: ncoalesce
    type(darrayReal), intent(inout)      :: div
    ! Local variables.
    integer(kind = int32)              :: allocateStatus
    integer(kind = int32)              :: chosen
    integer(kind = int32), allocatable :: eligiblepspop(:)
    integer(kind = int32)              :: numcoalesce
    integer(kind = int32)              :: popforps
    integer(kind = int32)              :: jpop
    integer(kind = int32)              :: numeligible
    real(kind = real32)                   :: x
    ! Allocate memory.
    allocate (eligiblepspop(activepop), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory for eligiblepspop!"
      ! Error, exit the program.
      stop
    end if
    numeligible = 0
    do jpop = 1, activepop
      if (numstrain(jpop) .le. 1) cycle
      numeligible = numeligible + 1
      eligiblepspop(numeligible) = jpop
    end do
    if (numeligible .gt. 0) then
      call randomNumber (x)
      chosen = int (x * numeligible, kind = int32) + 1
      popforps = eligiblepspop(chosen)
      numanctot = numanctot + 1
      numcoalesce = numstrain(popforps)
      call darraySet (ncoalesce, numanctot, numcoalesce)
      call darraySet (div, numanctot, time / lengthseq)
      numstrain(popforps) = 1
    end if
    ! Deallocate memory.
    deallocate (eligiblepspop, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) &
        "Failed to deallocate memory for eligiblepspop!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine doPeriodicSelection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Eligible returns the number of populations eligible to be the nascent
  !> and parent ecotypes, as well as eligible for periodic selection.
  !>
  !> @param[in]     numstrain     An array of the number of strains in each
  !>                                ecotype.
  !> @param[in]     activepop     The number of active populations.
  !> @param[out]    eligibleNI    The number of populations eligible to be
  !>                                the parent population.
  !> @param[out]    eligiblePS    The number of populations eligible for
  !>                                periodic selection.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eligible (numstrain, activepop, eligibleNI, eligiblePS)
    integer(kind = int32), intent(in)    :: numstrain(:)
    integer(kind = int32), intent(in)    :: activepop
    integer(kind = int32), intent(out)   :: eligibleNI
    integer(kind = int32), intent(out)   :: eligiblePS
    ! Local variables.
    integer(kind = int32) :: jpop
    eligibleNI = activepop
    if (activepop .eq. 1) then
      eligibleNI = 0
    end if
    eligiblePS = 0
    do jpop = 1, activepop
      if (numstrain(jpop) .ne. 1) eligiblePS = eligiblePS + 1
    end do
    return
  end subroutine eligible

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets an integer argument from the command line.
  !>
  !> @param[in]     arg           The argument number to get.
  !> @param[out]    out           The value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getIntegerArgument (arg, out)
    integer(kind = int32), intent(in)    :: arg
    integer(kind = int32), intent(out)   :: out
    ! Local variables.
    character(len = 100)  :: buffer
    integer(kind = int32) :: error
    call getarg (arg, buffer)
    read (unit = buffer, fmt = *, iostat = error) out
    if (error .ne. 0) then
      write (unit = *, fmt = *) &
        "ERROR: Integer number not supplied with argument ", arg
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine getIntegerArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets a logical argument from the command line.
  !>
  !> @param[in]     arg           The argument number to get.
  !> @param[out]    out           The value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getLogicalArgument (arg, out)
    integer(kind = int32), intent(in)    :: arg
    logical, intent(out)                 :: out
    ! Local variables.
    character(len = 100) :: buffer
    call getarg (arg, buffer)
    ! Unless we receive an appropriate true string, return false.
    out = .false.
    ! Catch lower case true string.
    if (trim(buffer) .eq. 'true' .or. trim(buffer) .eq. 't') then
      out = .true.
    ! Catch upper case true string.
    else if (trim(buffer) .eq. 'TRUE' .or. trim(buffer) .eq. 'T') then
      out = .true.
    ! Catch mixed case true string.
    else if (trim(buffer) .eq. 'True') then
      out = .true.
    ! Catch binary true string.
    else if (trim(buffer) .eq. '1') then
      out = .true.
    end if
    return
  end subroutine getLogicalArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets a real argument from the command line.
  !>
  !> @param[in]     arg           The argument number to get.
  !> @param[out]    out           The value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getRealArgument (arg, out)
    integer(kind = int32), intent(in)    :: arg
    real(kind = real32), intent(out)     :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer(kind = int32)              :: error
    call getarg (arg, buffer)
    read (unit = buffer, fmt = *, iostat = error) out
    if (error .ne. 0) then
      write (unit = *, fmt = *) &
        "ERROR: Real number not supplied with argument ", arg
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine getRealArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets a string argument from the command line.
  !>
  !> @param[in]     arg           The argument number to get.
  !> @param[out]    out           The value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getStringArgument (arg, out)
    integer(kind = int32), intent(in)    :: arg
    character(len = *), intent(out)      :: out
    call getarg (arg, out)
    return
  end subroutine getStringArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate the G statistic for a contingency table.
  !>
  !> @param[in]     nrows         The number of rows in the table.
  !> @param[in]     ncols         The number of columns in the table.
  !> @param[in]     x             The contingency table.
  !> @return                      The G statistic.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function gtest (nrows, ncols, x)  result (return_value)
    integer(kind = int32), intent(in)    :: nrows
    integer(kind = int32), intent(in)    :: ncols
    real(kind = real64), intent(in)      :: x(nrows, ncols)
    real(kind = real64)                  :: return_value
    ! Local variables.
    integer(kind = int32) :: i
    integer(kind = int32) :: j
    real(kind = real64)   :: g
    real(kind = real64)   :: n
    real(kind = real64)   :: sr(nrows)
    real(kind = real64)   :: sc(ncols)
    real(kind = real64)   :: e(nrows, ncols)
    ! Calculate the number of observations.
    n = sum (x)
    ! Make sure there is at least one positive observation.
    if (.not. n .gt. 0.0) then
      write (unit = *, fmt = *) "At least one entry of x must be positive."
      stop
    end if
    ! Make sure there are no negative observations.
    if (any (x .lt. 0.0)) then
      write (unit = *, fmt = *) "All entrys must be nonnegative."
      stop
    end if
    ! Calculate the row and column totals.
    sr = (/ (sum (x(i,:)), i = 1, nrows) /)
    sc = (/ (sum (x(:,i)), i = 1, ncols) /)
    ! Calculate the expected matrix.
    do i = 1, nrows
      do j = 1, ncols
        e(i,j) = sr(i) * sc(j) / n
      end do
    end do
    ! Calculate the G statistic.
    g = 0.0
    do i = 1, nrows
      do j = 1, ncols
        if (x(i,j) .gt. 0.0) g = g + x(i,j) * log (x(i,j) / e(i,j))
      end do
    end do
    g = 2.0 * g
    ! Return the G statistic value.
    return_value = g
    return
  end function gtest

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Takes the expected number of substitutions, and gives back the actual
  !> number of substitutions using the Poisson distribution.
  !>
  !> @param         xmean         The expected number of substitutions.
  !> @param         nsubs         The actual number of substitutions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson (xmean, nsubs)
    real(kind = real64), intent(in)      :: xmean
    integer(kind = int64), intent(out)   :: nsubs
    ! Local variables.
    integer(kind = int32) :: jmut
    real(kind = real64)   :: accumprob
    real(kind = real64)   :: expect
    real(kind = real64)   :: prob
    real(kind = real32)   :: x
    expect = xmean
    call randomNumber (x)
    accumprob = 0.0
    prob = exp (-1.0 * expect)
    do jmut = 0, 100
      if (jmut .ne. 0) then
        prob = (prob * expect) / jmut
      end if
      accumprob = accumprob + prob
      if (x .lt. accumprob) then
        nsubs = jmut
        return
      end if
    end do
    nsubs = int (expect, kind = int64)
    return
  end subroutine poisson

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deallocate memory reserved for the random number generator.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomClose ()
    ! Local variables.
    integer(kind = int32) :: allocateStatus
    ! Deallocate space for the RNG.
    deallocate (rng, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) "Failed to deallocate memory for the RNG!"
      ! Error, exit the program.
      stop
    end if
    return
  end subroutine randomClose

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialize the random number generator.
  !>
  !> @param[in]     iii           The integer used to seed the random number
  !>                                generator.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomInitialize (iii)
    integer(kind = int32), intent(in)    :: iii
    ! Local variables.
    integer(kind = int32) :: i
    integer(kind = int32) :: allocateStatus
    real(kind = real32)   :: x
    ! Allocate space for each thread's RNG.
    allocate (rng(numberThreads), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory for the RNG!"
      ! Error, exit the program.
      stop
    end if
    ! Seed the RNG of each thread.
    do i = 1, numberThreads
      call random_number (x)
      call ziggurat_seed (rng(i), int (iii * x, kind = int32))
    end do
  end subroutine randomInitialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generate a random number in the range [0,1).
  !>
  !> @param[out]    x             The random number generated
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomNumber (x)
    real(kind = real32), intent(out)     :: x
    ! Local variables.
    integer(kind = int32) :: n
#ifdef _OPENMP
    n = omp_get_thread_num() + 1
#else
    n = 1
#endif
    x = ziggurat_uni(rng(n))
    return
  end subroutine randomNumber

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This subroutine is set up to be the general purpose subprogram for
  !> various main programs:  brute force, hillclimb, and three confidence
  !> interval programs.
  !>
  !> @param[in]     omega         The rate of niche invasion per eligible
  !>                                parental population, measured as niche
  !>                                invasions per nucleotide substitution
  !>                                in a given gene.
  !> @param[in]     sigma         The rate of periodic selection per eligible
  !>                                population, measured as periodic selection
  !>                                events per population per nucleotide
  !>                                substitution in a given gene.
  !> @param[in]     npop          The number of ecotypes assumed to be in the
  !>                                environmental DNA sample.
  !> @param[in]     numcrit       The number of criteria for making cluster
  !>                                bins.
  !> @param[in]     nu            The number of homologous gene sequences in
  !>                                the environmental sample.
  !> @param[in]     nrep          The number of replicate simulations for a
  !>                                given set of values.
  !> @param[in]     lengthseq     The length of the sequences.
  !> @param[in]     realdata      The number of bins in the real data.
  !> @param[in]     crit          The list of bin identity levels (between 0.0
  !>                                and 1.0).
  !> @param[out]    avgsuccess    The average number of results within xxx%
  !>                                tolerance.
  !>                                1: 500% tolerance.
  !>                                2: 200% tolerance.
  !>                                3: 150% tolerance.
  !>                                4: 125% tolerance.
  !>                                5: 110% tolerance.
  !>                                6: 105% tolerance.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine runFredProgram (omega, sigma, npop, numcrit, nu, nrep, &
    lengthseq, realdata, crit, avgsuccess)
    real(kind = real64), intent(in)      :: omega
    real(kind = real64), intent(in)      :: sigma
    real(kind = real64), intent(out)     :: avgsuccess(6)
    integer(kind = int32), intent(in)    :: lengthseq
    integer(kind = int32), intent(in)    :: npop
    integer(kind = int32), intent(in)    :: nu
    integer(kind = int32), intent(in)    :: numcrit
    integer(kind = int32), intent(in)    :: nrep
    integer(kind = int32), intent(in)    :: realdata(:)
    real(kind = real32), intent(in)      :: crit(:)
    ! Local variables.
    integer(kind = int32) :: i
    integer(kind = int32) :: irep
    logical               :: success(6, nrep)
    ! Initialize the avgsuccess array with values of 0.
    avgsuccess = (/ (0.0d0, i = 1, 6) /)
    ! Make sure that omega has valid values.
    if (isnan (omega)) return ! XXX Should use ieee_is_nan when supported.
    if (omega .lt. epsilon (0.0d0)) return
    if (omega .gt. huge (0.0d0) - 1.0d0) return
    ! Make sure that sigma has valid values.
    if (isnan (sigma)) return ! XXX Should use ieee_is_nan when supported.
    if (sigma .lt. epsilon (0.0d0)) return
    if (sigma .gt. huge (0.0d0) - 1.0d0) return
    ! Make sure that npop does not exceed the number of homologous gene
    ! sequences in the environmental sample, and is greater than zero.
    if (npop .gt. nu .or. npop .le. 0) return
#ifdef _OPENMP
    call omp_set_num_threads (numberThreads)
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(irep)
    !$OMP DO
#endif
    ! Run the simulation nrep number of times.
    do irep = 1, nrep
      ! Run the simulation.
      call simulation (omega, sigma, npop, numcrit, nu, lengthseq, &
        realdata, crit, success(:, irep))
    end do
#ifdef _OPENMP
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
#endif
    ! Calculate the average success level for each tolerance range.
    do irep = 1, nrep
      do i = 1, 6
        if (success(i, irep)) avgsuccess(i) = avgsuccess(i) + 1.0d0
      end do
    end do
    avgsuccess = (/ (avgsuccess(i) / real (nrep, kind = real64), i = 1, 6) /)
    return
  end subroutine runFredProgram

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Run the population simulation.
  !>
  !> @param[in]     omega         The rate of niche invasion per eligible
  !>                                parental population, measured as niche
  !>                                invasions per nucleotide substitution
  !>                                in a given gene.
  !> @param[in]     sigma         The rate of periodic selection per eligible
  !>                                population, measured as periodic selection
  !>                                events per population per nucleotide
  !>                                substitution in a given gene.
  !> @param[in]     npop          The number of ecotypes assumed to be in the
  !>                                environmental DNA sample.
  !> @param[in]     numcrit       The number of criteria for making cluster
  !>                                bins.
  !> @param[in]     nu            The number of homologous gene sequences in
  !>                                the environmental sample.
  !> @param[in]     lengthseq     The length of the sequences.
  !> @param[in]     realdata      The number of bins in the real data.
  !> @param[in]     crit          The list of bin identity levels (between 0.0
  !>                                and 1.0).
  !> @param[out]    success       The results are within xxx% tolerance.
  !>                                1: 500% tolerance.
  !>                                2: 200% tolerance.
  !>                                3: 150% tolerance.
  !>                                4: 125% tolerance.
  !>                                5: 110% tolerance.
  !>                                6: 105% tolerance.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation (omega, sigma, npop, numcrit, nu, lengthseq, &
    realdata, crit, success)
    real(kind = real64), intent(in)      :: omega
    real(kind = real64), intent(in)      :: sigma
    integer(kind = int32), intent(in)    :: npop
    integer(kind = int32), intent(in)    :: numcrit
    integer(kind = int32), intent(in)    :: nu
    integer(kind = int32), intent(in)    :: lengthseq
    integer(kind = int32), intent(in)    :: realdata(:)
    real(kind = real32), intent(in)      :: crit(:)
    logical, intent(out)                 :: success(6)
    ! Local variables.
    integer(kind = int32)                :: activepop
    integer(kind = int32)                :: allocateStatus
    integer(kind = int32)                :: bin(numcrit)
    integer(kind = int32)                :: event
    integer(kind = int32)                :: jpop
    integer(kind = int32)                :: ntotalpop
    integer(kind = int32)                :: numanctot
    integer(kind = int32), allocatable   :: numstrain(:)
    real(kind = real32)                  :: time
    type(darrayInteger) :: ncoalesce
    type(darrayReal)    :: div
    ! Initialize the dynamic arrays.
    call darrayInitialize(ncoalesce, 4 * nu)
    call darrayInitialize(div, 4 * nu)
    ! Allocate memory for the static array.
    allocate (numstrain(nu), stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory for numstrain!"
      ! Error, exit the program.
      stop
    end if
    ! Initialize the population simulation.
    call startpops (npop, nu, activepop, numstrain, div, numanctot, time)
    do
      ! subroutine whichEventAndWhen returns the key event, and gives the
      ! waiting time until the event (since the previous key event).
      call whichEventAndWhen (omega, sigma, numstrain, activepop, time, event)
      ! Run the chosen event.
      select case (event)
        case (EVENT_NICHE_INVASION)
          call doNicheInvasion (activepop, numanctot, numstrain, ncoalesce, &
            lengthseq, time, div)
        case (EVENT_PERIODIC_SELECTION)
          call doPeriodicSelection (activepop, numanctot, numstrain, &
            ncoalesce, lengthseq, time, div)
        case default
          write (unit = *, fmt = *) "Unknown event!"
          ! Error, exit the program.
          stop
      end select
      ! Calculate the total number of organisms as ntotalpop
      ntotalpop = 0
      do jpop = 1, activepop
        ntotalpop = ntotalpop + numstrain(jpop)
      end do
      if (ntotalpop .le. 1) exit
    end do
    call binning (bin, div, crit, ncoalesce, numcrit, numanctot)
    ! Check whether a replicate is within a factor of (5.00, 2.00, 1.50, 1.25,
    ! 1.10, 1.05) for all bins
    call testForSuccessFit (bin, numcrit, realdata, success)
    ! Close the dynamic arrays.
    call darrayClose(ncoalesce)
    call darrayClose(div)
    ! Deallocate memory for the static array.
    deallocate (numstrain, stat = allocateStatus)
    if (allocateStatus .gt. 0) then
      write (unit = *, fmt = *) "Failed to deallocate memory for numstrain!"
      ! Error, exit the program.
      stop
    end if
  end subroutine simulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Start the population simulation.
  !>
  !> @param[in]     npop           The number of ecotypes assumed to be in
  !>                                the environmental DNA sample.
  !> @param[in]     nu             The number of homologous gene sequences
  !>                                in the environmental sample.
  !> @param[out]    activepop      The number of active populations.
  !> @param[in,out] numstrain      An array of the number of strains in each
  !>                                ecotype.
  !> @param[in,out] div            The number of substitutions per site for
  !>                                each ancestor.
  !> @param[out]    numanctot      The total number of ancestors, over all
  !>                                strains.
  !> @param[out]    time           The time of the last ancestor described.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine startpops (npop, nu, activepop, numstrain, div, numanctot, time)
    integer(kind = int32), intent(in)    :: npop
    integer(kind = int32), intent(in)    :: nu
    integer(kind = int32), intent(out)   :: activepop
    integer(kind = int32), intent(out)   :: numstrain(:)
    integer(kind = int32), intent(out)   :: numanctot
    real(kind = real32), intent(out)     :: time
    type(darrayReal), intent(inout)      :: div
    ! Local variables.
    integer(kind = int32) :: i
    ! activepop is the number of active populations (note some will disappear
    ! in the backwards simulation as they become invented in reverse)
    activepop = npop
    numanctot = nu
    call canonical (npop, nu, numstrain)
    ! numstrain(i) is the number of strains in ecotype i; npop is the number
    ! of ecotypes
    ! Time is the time of the last ancestor described; it starts at 0
    time = 0.0
    ! div(i) is the number of substitutions per site down to ancestor i from
    ! the present the contemporary nu strains are considered ancestors 1
    ! through nu, and their divergence time is set to 0
    do i = 1, nu
      call darraySet (div, i, 0.0)
    end do
    return
  end subroutine startpops

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Adds one to each succcess level when a replicate is within a factor of
  !> xxx% for all bins.
  !>
  !> @param[in]     bin           The number of clusters at the identity level
  !>                                sought.
  !> @param[in]     numcrit       The number of criteria for making cluster
  !>                                bins.
  !> @param[in]     realdata      The number of bins in the real data.
  !> @param[out]    success       The results within xxx% tolerance.
  !>                                1: 500% tolerance.
  !>                                2: 200% tolerance.
  !>                                3: 150% tolerance.
  !>                                4: 125% tolerance.
  !>                                5: 110% tolerance.
  !>                                6: 105% tolerance.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testForSuccessFit (bin, numcrit, realdata, success)
    integer(kind = int32), intent(in)    :: bin(:)
    integer(kind = int32), intent(in)    :: numcrit
    integer(kind = int32), intent(in)    :: realdata(:)
    logical, intent(out)                 :: success(6)
    ! Local variables.
    integer(kind = int32) :: i
    integer(kind = int32) :: jcrit
    real(kind = real32)   :: xbin
    real(kind = real32)   :: xreal
    ! Initialize the success array.
    success = (/ (.false., i = 1, 6) /)
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 5.0 .or. xbin / xreal .gt. 5.0) return
    end do
    ! Meets the 500% criterion.
    success(1) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 2.0 .or. xbin / xreal .gt. 2.0) return
    end do
    ! Meets the 200% criterion.
    success(2) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.5 .or. xbin / xreal .gt. 1.5) return
    end do
    ! Meets the 150% criterion.
    success(3) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.25 .or. xbin / xreal .gt. 1.25) return
    end do
    ! Meets the 125% criterion.
    success(4) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.1 .or. xbin / xreal .gt. 1.1) return
    end do
    ! Meets the 110% criterion.
    success(5) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.05 .or. xbin / xreal .gt. 1.05) return
    end do
    ! Meets the 105% criterion.
    success(6) = .true.
    return
  end subroutine testForSuccessFit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns the type of the key event, and gives the waiting time until the
  !> event (since the previous key event).
  !>
  !> @param[in]     omega         The rate of niche invasion per eligible
  !>                                parental population, measured as niche
  !>                                invasions per nucleotide substitution in a
  !>                                given gene.
  !> @param[in]     sigma         The rate of periodic selection per eligible
  !>                                population, measured as periodic selection
  !>                                events per population per nucleotide
  !>                                substitution in a given gene.
  !> @param[in]     numstrain     An array of the number of strains in each
  !>                                ecotype.
  !> @param[in]     activepop     The number of active populations.
  !> @param[in,out] time          The total time (in substitutions) before the
  !>                                present for a coalescence event.
  !> @param[out]    event         The type of key event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine whichEventAndWhen (omega, sigma, numstrain, activepop, time, &
    event)
    real(kind = real64), intent(in)      :: omega
    real(kind = real64), intent(in)      :: sigma
    integer(kind = int32), intent(in)    :: numstrain(:)
    integer(kind = int32), intent(in)    :: activepop
    integer(kind = int32), intent(out)   :: event
    real(kind = real32), intent(inout)   :: time
    ! Local variables.
    integer(kind = int32) :: eligibleNI
    integer(kind = int32) :: eligiblePS
    integer(kind = int64) :: nsubs
    real(kind = real64)   :: effectiveOmega
    real(kind = real64)   :: effectiveSigma
    real(kind = real64)   :: rateKey
    real(kind = real64)   :: timeWait
    real(kind = real32)   :: x
    ! subroutine eligible returns the number of populations eligible for
    ! niche invasion and periodic selection.
    call eligible (numstrain, activepop, eligibleNI, eligiblePS)
    ! note that omega is a per population rate of niche invasion,
    ! so we multiply by the number of eligible parent populations
    effectiveOmega = eligibleNI * omega
    effectiveSigma = eligiblePS * sigma
    ! rateKey is the total rate of all key events.
    rateKey = effectiveOmega + effectiveSigma
    ! the value of timewait yielded here is the expected number of
    ! nucleotide substitutions
    call randomNumber (x)
    if (x .lt. 1.0e-6) x = 1.0e-6
    timeWait = -1 * log (x) / rateKey
    ! we next get the actual number of substitutions nsubs using poisson
    call poisson (timeWait, nsubs)
    ! time is the total time (in substitutions) before the present for a
    ! coalescence event
    time = time + nsubs
    ! now, what kind of key event
    call randomNumber (x)
    if (x .lt. effectiveOmega / rateKey) then
      event = EVENT_NICHE_INVASION
    else
      event = EVENT_PERIODIC_SELECTION
    end if
    return
  end subroutine whichEventAndWhen

end module methods
