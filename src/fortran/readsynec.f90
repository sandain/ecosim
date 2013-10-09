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
!> The readsynec program extracts sequence data and the name of each strain
!> into separate files, making sure that the sequence data is in lower case.
!>
!> @pre  Requires that the fasta.dat file exists and that it contains
!>         sequence data.
!>
!> @post Generates the population.dat and nameofstrains.dat files that contain
!>         just the sequence data, and the name of each strain respectively.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program readsynec
  use methods
  implicit none
  character(len = 20)             :: sequenceFormat
  character(len = 20)             :: strainnameFormat
  character(len = 256)            :: fastaFile
  character(len = 256)            :: populationFile
  character(len = 256)            :: namesofstrainsFile
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), allocatable :: strainname(:,:)
  integer                         :: allocateStatus
  integer                         :: i
  integer                         :: jallele
  integer                         :: idigit
  integer                         :: knuc
  integer                         :: numstrain
  integer                         :: length
  integer, parameter              :: strainnameLength = 20
  integer, parameter              :: fastaUnit = 1
  integer, parameter              :: populationUnit = 2
  integer, parameter              :: nameofstrainsUnit = 3
  logical                         :: fileExists
  ! Provide default file names to use.
  fastaFile = 'removegaps.dat'
  populationFile = 'population.dat'
  namesofstrainsFile = 'namesofstrains.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, fastaFile)
      case (2)
        call getArgument (2, populationFile)
      case (3)
        call getArgument (3, namesofstrainsFile)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: fasta.dat population.dat namesofstrains.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input file exist.
  inquire (file = trim (fastaFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The fasta.dat file was not found at: ", &
      trim (fastaFile)
    ! Error, exit the program.
    stop
  end if
  ! Open Input file
  open (unit = fastaUnit, file = trim (fastaFile), access = 'sequential', &
    form = 'formatted')
  ! Open Output files
  open (unit = populationUnit, file = trim (populationFile), &
    access = 'sequential', form = 'formatted')
  open (unit = nameofstrainsUnit, file = trim (namesofstrainsFile), &
    access = 'sequential', form = 'formatted')
  ! Read in the fasta file and output to the population and nameofstrains files
  read (unit = fastaUnit, fmt = *) numstrain, length
  write (unit = populationUnit, fmt = *) numstrain, length
  write (unit = nameofstrainsUnit, fmt = *) numstrain, length
  ! Allocate Variables.
  allocate (sequence(numstrain, length), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  allocate (strainname(numstrain, strainnameLength), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: strainname!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = strainnameFormat, fmt = "(a1,i0,a3)") "(", strainnameLength, &
    "a1)"
  write (unit = sequenceFormat, fmt = "(a1,i0,a3)") "(", length, "a1)"
  do jallele = 1, numstrain
    read (unit = fastaUnit, fmt = trim (strainnameFormat)) &
      (strainname(jallele, idigit), idigit = 1, strainnameLength)
    write (unit = nameofstrainsUnit, fmt = trim (strainnameFormat)) &
      (strainname(jallele, idigit), idigit = 2, strainnameLength)
    read (unit = fastaUnit, fmt = trim (sequenceFormat)) &
      (sequence(jallele, knuc), knuc = 1, length)
    do knuc = 1, length
      if (sequence(jallele, knuc) .eq. 'A') sequence(jallele, knuc) = 'a'
      if (sequence(jallele, knuc) .eq. 'G') sequence(jallele, knuc) = 'g'
      if (sequence(jallele, knuc) .eq. 'T') sequence(jallele, knuc) = 't'
      if (sequence(jallele, knuc) .eq. 'C') sequence(jallele, knuc) = 'c'
    end do
    write (unit = populationUnit, fmt = trim (sequenceFormat)) &
      (sequence(jallele, knuc), knuc = 1, length)
  end do
  ! Deallocate memory.
  deallocate (sequence, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  deallocate (strainname, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: strainname!"
    ! Error, exit the program.
    stop
  end if
  ! Close data files.
  close (unit = fastaUnit)
  close (unit = populationUnit)
  close (unit = nameofstrainsUnit)
  ! Successful termination of program.
  stop
end program readsynec
