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
!> The correctpcr program corrects for PCR error by assuming that each error
!> would be a unique nucleotide at one site.  Of all the sites with a unique
!> nucleotide, we randomly choose some of them to be corrected.  Each strain
!> chosen to be corrected has the offending nucleotide modified to match that
!> of the strain most closely related.
!>
!> @pre  Requires that the population.dat and pcrerror.dat files exist and
!>         that they contain all of the variables necessary to perform the
!>         correction.
!>
!> @post Generates the correctpcr.out file that contains the corrected
!>         sequences.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program correctpcr
  use methods
  use darray
  use tmatrix
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), parameter   :: code(4) = (/'a', 't', 'g', 'c'/)
  character(len = 256)            :: populationFile
  character(len = 256)            :: pcrerrorFile
  character(len = 256)            :: outputFile
  integer                         :: allocateStatus
  integer                         :: numnuc(4)
  integer                         :: i
  integer                         :: iii
  integer                         :: num
  integer                         :: numtochange
  integer                         :: closest
  integer                         :: seqchange
  integer                         :: nucchange
  integer                         :: numsingletons
  integer                         :: iother
  integer                         :: jstrain
  integer                         :: knuc
  integer                         :: nucnuc
  integer                         :: jchange
  integer                         :: numstrain
  integer                         :: length
  type(darrayInteger)             :: singletonStrain
  type(darrayInteger)             :: singletonNuc
  integer, allocatable            :: tochange(:,:)
  integer, parameter              :: populationUnit = 1
  integer, parameter              :: outputUnit = 2
  integer, parameter              :: pcrerrorUnit = 4
  logical                         :: fileExists
  real                            :: x
  real                            :: diff
  real                            :: xmindiff
  real                            :: pcrerror
  type(tmatrixReal)               :: divergematrix
  ! Provide default file names to use.
  populationFile = 'population.dat'
  pcrerrorFile = 'pcrerror.dat'
  outputFile = 'correctpcr.out'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, populationFile)
      case (2)
        call getArgument (2, pcrerrorFile)
      case (3)
        call getArgument (3, outputFile)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: population.dat pcrerror.dat correctpcr.out"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input files exist.
  inquire (file = trim (populationFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The population.dat file was not found at: ", &
      trim (populationFile)
    ! Error, exit the program.
    stop
  end if
  inquire (file = trim (pcrerrorFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The pcrerror.dat file was not found at: ", &
      trim (pcrerrorFile)
    ! Error, exit the program.
    stop
  end if
  ! Open input files
  open (unit = populationUnit, file = trim (populationFile), &
    access = 'sequential', form = 'formatted')
  open (unit = pcrerrorUnit, file = trim (pcrerrorFile), &
    access = 'sequential', form = 'formatted')
  ! Open output file
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! Read in population file
  read (unit = populationUnit, fmt = *) numstrain, length
  ! numstrain is the number of strains in the sample
  ! length is the number of nucleotide sites
  ! Allocate variables.
  call darrayInitialize(singletonStrain, 4 * numstrain)
  call darrayInitialize(singletonNuc, 4 * numstrain)
  call tmatrixInitialize(divergematrix, numstrain)
  allocate (sequence(numstrain, length), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  allocate (tochange(numstrain, 2), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: tochange!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  ! read in the sequences
  do jstrain = 1, numstrain
    read (unit = populationUnit, fmt = trim(sequence_format)) &
      (sequence(jstrain, knuc), knuc = 1, length)
  end do
  ! Read in pcrerrorFile
  ! pcrerror is the per nucleotide site error in pcr; according to Akmaev
  ! and Wang, it is 1.37e-4, or 1/7300
  read (unit = pcrerrorUnit, fmt = *) pcrerror
  ! iii is the odd random number seed
  read (unit = pcrerrorUnit, fmt = *) iii
  ! Initialize the random number generator.
  call randomInitialize (iii)
  call diverge (numstrain, length, sequence, divergematrix)
  numsingletons = 0
  do knuc = 1, length
    ! numnuc(1) = number of a's at site knuc
    ! numnuc(2) = number of t's at site knuc
    ! numnuc(3) = number of g's at site knuc
    ! numnuc(4) = number of c's at site knuc
    numnuc = (/ (0, i = 1, 4) /)
    do nucnuc = 1, 4
      do jstrain = 1, numstrain
        if (sequence(jstrain, knuc) .eq. code(nucnuc)) then
          numnuc(nucnuc) = numnuc(nucnuc) + 1
        end if
      end do
      if (numnuc(nucnuc) .eq. 1) then
        do jstrain = 1, numstrain
          if (sequence(jstrain, knuc) .eq. code(nucnuc)) then
            numsingletons = numsingletons + 1
            call darraySet (singletonStrain, numsingletons, jstrain)
            call darraySet (singletonNuc, numsingletons, knuc)
            exit
          endif
        end do
      endif
    end do
  end do
  ! now, how many singletons to correct for pcr error?
  numtochange = nint (pcrerror * length * numstrain)
  if (numtochange .gt. numsingletons) then
    numtochange = numsingletons
  end if
  ! change the singletons
  if (numtochange .gt. 0) then
    do jchange = 1, numtochange
      call randomNumber (x)
      num = int (x * numsingletons) + 1
      tochange(jchange, 1) = darrayGet(singletonStrain, num)
      tochange(jchange, 2) = darrayGet(singletonNuc, num)
      call darraySet(singletonStrain, num, darrayGet(singletonStrain, numsingletons))
      call darraySet(singletonNuc, num, darrayGet(singletonNuc, numsingletons))
      numsingletons = numsingletons - 1
    end do
    ! now, change each of the chosen singletons to the value of the most
    ! similar strain
    closest = 1
    do jchange = 1, numtochange
      xmindiff = 1.0
      seqchange = tochange(jchange, 1)
      nucchange = tochange(jchange, 2)
      do iother = 1, numstrain
        if (seqchange .eq. iother) cycle
        diff = tmatrixGet(divergematrix, seqchange, iother)
        if (diff .lt. xmindiff) then
          closest = iother
          xmindiff = diff
        endif
      end do
      sequence(seqchange, nucchange) = sequence(closest, nucchange)
    end do
  end if
  ! Output to output file.
  write (unit = outputUnit, fmt = *) numstrain, length
  do jstrain = 1, numstrain
    write (unit = outputUnit, fmt = trim(sequence_format)) &
      (sequence(jstrain, knuc), knuc = 1, length)
  end do
  ! Deallocate memory.
  call darrayClose(singletonStrain)
  call darrayClose(singletonNuc)
  call tmatrixClose(divergematrix)
  deallocate (sequence, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  deallocate (tochange, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: tochange!"
    ! Error, exit the program.
    stop
  end if
  ! Close the random number generator.
  call randomClose ()
  ! Close data files.
  close (unit = populationUnit)
  close (unit = pcrerrorUnit)
  close (unit = outputUnit)
  ! Successful termination of program.
  stop
end program correctpcr
