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
!> The divergencematrix program generates a divergence matrix of the provided
!> sequence data.
!>
!> @pre  Requires that the correctpcr.out file exists and that it contains the
!>         sequence data necessary to generate the divergence matrix.
!> @post Generates the divergencematrix.dat file that contains the divergence
!>         matrix.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program divergencematrix
  use methods
  use tmatrix
  implicit none
  character(len = 20)             :: sequenceFormat
  character(len = 20)             :: outputFormat
  character(len = 256)            :: inputFile
  character(len = 256)            :: outputFile
  character(len = 1), allocatable :: sequence(:,:)
  integer                         :: allocateStatus
  integer                         :: numstrain
  integer                         :: length
  integer                         :: i
  integer                         :: jstrain
  integer                         :: knuc
  integer                         :: kstrain
  integer, parameter              :: inputUnit = 1
  integer, parameter              :: outputUnit = 2
  logical                         :: file_exists
  type(tmatrixReal)               :: divergematrix
  ! Provide default file names to use.
  inputFile = 'correctpcr.out'
  outputFile = 'divergencematrix.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, inputFile)
      case (2)
        call getArgument (2, outputFile)
      case (3)
        call getArgument (3, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: correctpcr.out divergencematrix.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input file exist.
  inquire (file = trim (inputFile), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The correctpcr.out file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Open Input file
  open (unit = inputUnit, file = trim (inputFile), &
    access = 'sequential', form = 'formatted')
  ! Open Output file
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! Read in input file.
  read (unit = inputUnit, fmt = *) numstrain, length
  ! numstrain is the number of strains in the sample;
  ! length is the number of nucleotide sites
  ! Allocate variables.
  call tmatrixInitialize(divergematrix, numstrain)
  allocate (sequence(numstrain, length), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = sequenceFormat, fmt = "(a1,i0,a3)") "(", length, "a1)"
  write (unit = outputFormat, fmt = "(a1,i0,a10)") "(", numstrain, &
    "(1x,f7.4))"
  ! Read in the sequence data
  do jstrain = 1, numstrain
    read (unit = inputUnit, fmt = trim (sequenceFormat)) &
      (sequence(jstrain, knuc), knuc = 1, length)
  end do
  call diverge (numstrain, length, sequence, divergematrix)
  ! Output to output file.
  write (unit = outputUnit, fmt = *) numstrain, length
  do jstrain = 1, numstrain
    write (unit = outputUnit, fmt = trim (outputFormat)) &
      (1.0 - tmatrixGet(divergematrix, jstrain, kstrain), kstrain = 1, numstrain)
  end do
  ! Deallocate memory.
  call tmatrixClose(divergematrix)
  deallocate (sequence, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  ! Close data files.
  close (unit = inputUnit)
  close (unit = outputUnit)
  ! Successful termination of program.
  stop
end program divergencematrix
