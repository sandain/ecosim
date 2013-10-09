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
!> The methods module contains all of the common subroutines used in the
!> various programs of the Ecotype Simulation application.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module methods
  use darray
  use tmatrix
  use ziggurat
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  ! Declare public methods.
  public :: diverge
  public :: getArgument
  public :: randomClose
  public :: randomInitialize
  public :: randomNumber
  public :: runProgram

  ! Declare public variables.
  logical, public :: debug = .false.       !< Display debug information.
  integer, public :: numberThreads = 1     !< The number of threads to start.

  ! Declare private global parameters.
  integer, parameter :: &
    EVENT_NICHE_INVASION     = 1001        !< Niche invasion event.
  integer, parameter :: &
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
    integer, intent(out)               :: bin(:)
    integer, intent(in)                :: numcrit
    integer, intent(in)                :: numanctot
    real, intent(in)                   :: crit(:)
    type(darrayInteger), intent(in)    :: ncoalesce
    type(darrayReal), intent(in)       :: div
    ! Local variables.
    integer :: lused
    integer :: janc
    integer :: jcrit
    integer :: nbins
    real    :: criterion
    real    :: div2
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
    integer, intent(in)                :: npop
    integer, intent(in)                :: nu
    integer, intent(out)               :: numstrain(:)
    ! Local variables.
    integer :: upto
    integer :: ipop
    integer :: istrain
    integer :: numberofstrains
    real    :: x
    numstrain(1) = nu
    if (npop .gt. 1) then
      do upto = 2, npop
        ! for the upto - 1 population, choose a random population (ipop)
        do
          call randomNumber (x)
          ipop = int (x * (upto - 1)) + 1
          numberofstrains = numstrain(ipop)
          if (numberofstrains .gt. 1) exit
        end do
        ! next, pick a random breakpoint
        call randomNumber (x)
        istrain = int (x * (numberofstrains - 1)) + 1
        numstrain(ipop) = istrain
        numstrain(upto) = numberofstrains - istrain
      end do
    end if
    return
  end subroutine canonical

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a divergence matrix from the provided sequences.
  !>
  !> @param[in]     numstrains    The number of strains in the sample.
  !> @param[in]     length        The length of the sequences in the sample.
  !> @param[in]     seq           An array containing the sequences in the
  !>                                sample.
  !> @param[in,out] divergematrix The resulting divergence matrix.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diverge (numstrains, length, seq, divergematrix)
    character(len = 1), intent(in)     :: seq(:,:)
    integer, intent(in)                :: numstrains
    integer, intent(in)                :: length
    type(tmatrixReal), intent(inout)   :: divergematrix
    ! Local variables.
    integer :: jstrain, kstrain, knuc
    real    :: diff
    do jstrain = 1, numstrains
      do kstrain = 1, jstrain - 1
        diff = 0.0
        do knuc = 1, length
          if (seq(jstrain, knuc) .ne. seq(kstrain, knuc)) then
            diff = diff + 1.0
          end if
        end do
        call tmatrixSet(divergematrix, jstrain, kstrain, diff / length)
      end do
      call tmatrixSet(divergematrix, jstrain, jstrain, 0.0)
    end do
    return
  end subroutine diverge

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
    integer, intent(inout)             :: activepop
    integer, intent(inout)             :: numanctot
    integer, intent(inout)             :: numstrain(:)
    integer, intent(in)                :: lengthseq
    real, intent(in)                   :: time
    type(darrayInteger), intent(inout) :: ncoalesce
    type(darrayReal), intent(inout)    :: div
    ! Local variables.
    integer :: popfornascent
    real    :: x
    if (activepop .eq. 0) return
    ! Choose the population for the nascent population
    call randomNumber (x)
    popfornascent = int (x * activepop) + 1
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
    integer, intent(in)                :: activepop
    integer, intent(inout)             :: numstrain(:)
    integer, intent(inout)             :: numanctot
    integer, intent(in)                :: lengthseq
    real, intent(in)                   :: time
    type(darrayInteger), intent(inout) :: ncoalesce
    type(darrayReal), intent(inout)    :: div
    ! Local variables.
    integer              :: allocateStatus
    integer              :: chosen
    integer, allocatable :: eligiblepspop(:)
    integer              :: numcoalesce
    integer              :: popforps
    integer              :: jpop
    integer              :: numeligible
    real                 :: x
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
      chosen = int (x * numeligible) + 1
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
    integer, intent(in)                :: numstrain(:)
    integer, intent(in)                :: activepop
    integer, intent(out)               :: eligibleNI
    integer, intent(out)               :: eligiblePS
    ! Local variables.
    integer :: jpop
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
    integer, intent(in)               :: arg
    integer, intent(out)              :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer              :: error
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
    integer, intent(in)                :: arg
    logical, intent(out)               :: out
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
    integer, intent(in)                :: arg
    real, intent(out)                  :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer              :: error
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
    integer, intent(in)                :: arg
    character(len = *), intent(out)    :: out
    call getarg (arg, out)
    return
  end subroutine getStringArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Takes the expected number of substitutions, and gives back the actual
  !> number of substitutions using the Poisson distribution.
  !>
  !> @param         xmean         The expected number of substitutions.
  !> @param         nsubs         The actual number of substitutions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson (xmean, nsubs)
    double precision, intent(in)       :: xmean
    integer, intent(out)               :: nsubs
    ! Local variables.
    integer          :: jmut
    double precision :: accumprob
    double precision :: expect
    double precision :: prob
    real             :: x
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
    nsubs = int (expect)
    return
  end subroutine poisson

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deallocate memory reserved for the random number generator.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomClose ()
    ! Local variables.
    integer :: allocateStatus
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
    integer, intent(in)                :: iii
    ! Local variables.
    integer :: i
    integer :: allocateStatus
    real    :: x
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
      call ziggurat_seed (rng(i), int (iii * x))
    end do
  end subroutine randomInitialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generate a random number in the range [0,1).
  !>
  !> @param[out]    x             The random number generated
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomNumber (x)
    real, intent(out)                  :: x
    ! Local variables.
    integer :: n
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
  subroutine runProgram (omega, sigma, npop, numcrit, nu, nrep, &
    lengthseq, realdata, crit, avgsuccess)
    double precision, intent(in)       :: omega
    double precision, intent(in)       :: sigma
    double precision, intent(out)      :: avgsuccess(6)
    integer, intent(in)                :: lengthseq
    integer, intent(in)                :: npop
    integer, intent(in)                :: nu
    integer, intent(in)                :: numcrit
    integer, intent(in)                :: nrep
    integer, intent(in)                :: realdata(:)
    real, intent(in)                   :: crit(:)
    ! Local variables.
    integer :: i
    integer :: irep
    logical :: success(6, nrep)
    ! Initialize the avgsuccess array with values of 0.
    avgsuccess = (/ (0.0d0, i = 1, 6) /)
    ! Make sure that omega and sigma have valid values.
    if (omega .ne. omega .or. omega .le. epsilon (0.0d0) .or. &
      omega .ge. huge (0.0d0) - 1.0d0) return
    if (sigma .ne. sigma .or. sigma .le. epsilon (0.0d0) .or. &
      sigma .ge. huge (0.0d0) - 1.0d0) return
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
    avgsuccess = (/ (avgsuccess(i) / real (nrep), i = 1, 6) /)
    return
  end subroutine runProgram

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
    double precision, intent(in)       :: omega
    double precision, intent(in)       :: sigma
    integer, intent(in)                :: npop
    integer, intent(in)                :: numcrit
    integer, intent(in)                :: nu
    integer, intent(in)                :: lengthseq
    integer, intent(in)                :: realdata(:)
    real, intent(in)                   :: crit(:)
    logical, intent(out)               :: success(6)
    ! Local variables.
    integer              :: activepop
    integer              :: allocateStatus
    integer              :: bin(numcrit)
    integer              :: event
    integer              :: jpop
    integer              :: ntotalpop
    integer              :: numanctot
    integer, allocatable :: numstrain(:)
    real                 :: time
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
    integer, intent(in)                :: npop
    integer, intent(in)                :: nu
    integer, intent(out)               :: activepop
    integer, intent(out)               :: numstrain(:)
    integer, intent(out)               :: numanctot
    real, intent(out)                  :: time
    type(darrayReal), intent(inout)    :: div
    ! Local variables.
    integer :: i
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
    integer, intent(in)                :: bin(:)
    integer, intent(in)                :: numcrit
    integer, intent(in)                :: realdata(:)
    logical, intent(out)               :: success(6)
    ! Local variables.
    integer :: i
    integer :: jcrit
    real    :: xbin
    real    :: xreal
    ! Initialize the success array.
    success = (/ (.false., i = 1, 6) /)
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal / xbin .gt. 5.0 .or. xbin / xreal .gt. 5.0) return
    end do
    ! Meets the 500% criterion.
    success(1) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal / xbin .gt. 2.0 .or. xbin / xreal .gt. 2.0) return
    end do
    ! Meets the 200% criterion.
    success(2) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal / xbin .gt. 1.5 .or. xbin / xreal .gt. 1.5) return
    end do
    ! Meets the 150% criterion.
    success(3) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal / xbin .gt. 1.25 .or. xbin / xreal .gt. 1.25) return
    end do
    ! Meets the 125% criterion.
    success(4) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal / xbin .gt. 1.1 .or. xbin / xreal .gt. 1.1) return
    end do
    ! Meets the 110% criterion.
    success(5) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
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
    double precision, intent(in)       :: omega
    double precision, intent(in)       :: sigma
    integer, intent(in)                :: numstrain(:)
    integer, intent(in)                :: activepop
    integer, intent(out)               :: event
    real, intent(inout)                :: time
    ! Local variables.
    integer          :: eligibleNI
    integer          :: eligiblePS
    integer          :: nsubs
    double precision :: effectiveOmega
    double precision :: effectiveSigma
    double precision :: rateKey
    double precision :: timeWait
    real             :: x
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
