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


program removegaps
  use methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 20)             :: strainname_format
  character(len = 256)            :: sequencesfasta_txt
  character(len = 256)            :: numbers_dat
  character(len = 256)            :: fasta_dat
  character(len = 1), allocatable :: strainname(:,:)
  character(len = 1), allocatable :: sequence(:,:)
  logical                         :: badnuc
  integer                         :: allocate_status
  integer                         :: numstrain
  integer                         :: length
  integer                         :: jnuc
  integer                         :: knuc
  integer                         :: mnuc
  integer                         :: jseq
  integer                         :: idigit
  integer, parameter              :: max_strainname = 20
  integer, parameter              :: sequencesfasta_unit = 1
  integer, parameter              :: numbers_unit = 2
  integer, parameter              :: fasta_unit = 3
  ! Read in command line arguments
  ! Expect path and filename for: sequencesfasta.txt numbers.dat fasta.dat
  if (iargc() .ge. 3) then
    call getArgument (1, sequencesfasta_txt)
    call getArgument (2, numbers_dat)
    call getArgument (3, fasta_dat)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: sequencesfasta.txt numbers.dat fasta.dat"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 4) then
    call getArgument (4, debug)
  else
    debug = .false.
  end if
  ! Open Input files
  open (unit = sequencesfasta_unit, file = trim (sequencesfasta_txt), access = 'sequential', form = 'formatted')
  open (unit = numbers_unit, file = trim (numbers_dat), access = 'sequential', form = 'formatted')
  ! Open Output file
  open (unit = fasta_unit, file = trim (fasta_dat), access = 'sequential', form = 'formatted')
  ! Read in numbers_dat
  read (unit = numbers_unit, fmt = *) numstrain, length
  ! length is the number of nucleotides
  ! numstrain is the number of sequences
  ! Allocate Variables.
  allocate (sequence(numstrain, length), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (strainname(numstrain, max_strainname), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  ! Create the format for printing sequence data.
  write (unit = strainname_format, fmt = "(a1,i0,a3)") "(", max_strainname, "a1)"
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  ! Read in sequencesfasta_txt
  do jseq = 1, numstrain
    read (unit = sequencesfasta_unit, fmt = trim(strainname_format)) (strainname(jseq, idigit), idigit = 1, max_strainname)
    read (unit = sequencesfasta_unit, fmt = trim(sequence_format)) (sequence(jseq, jnuc), jnuc = 1, length)
    do mnuc = 1, length
      if (sequence(jseq, mnuc) .eq. 'A') sequence(jseq, mnuc) = 'a'
      if (sequence(jseq, mnuc) .eq. 'G') sequence(jseq, mnuc) = 'g'
      if (sequence(jseq, mnuc) .eq. 'T') sequence(jseq, mnuc) = 't'
      if (sequence(jseq, mnuc) .eq. 'C') sequence(jseq, mnuc) = 'c'
    end do
  end do
  ! ok, so we've read in all the sequences.  Now we get rid of
  ! all nucleotide sites in which even a single sequence has a '-'
  knuc = 0
  ! knuc is the number of nucleotides that have been found to be gap free
  do jnuc = 1, length
    badnuc = .false.
    do jseq = 1, numstrain
      if (sequence(jseq, jnuc) .ne. 'a' .and. sequence(jseq, jnuc) .ne. 'g' .and. &
        sequence(jseq, jnuc) .ne. 't' .and. sequence(jseq, jnuc) .ne. 'c') then
        badnuc = .true.
        exit
      end if
    end do
    if (badnuc) cycle
    knuc = knuc + 1
    do jseq = 1, numstrain
      sequence(jseq, knuc) = sequence(jseq, jnuc)
    end do
  end do
  ! Output to fasta_dat
  write (unit = fasta_unit, fmt = *) numstrain, knuc
  do jseq = 1, numstrain
    write (unit = fasta_unit, fmt = trim(strainname_format)) (strainname(jseq, idigit), idigit = 1, max_strainname)
    write (unit = fasta_unit, fmt = trim(sequence_format)) (sequence(jseq, jnuc), jnuc = 1, knuc)
  end do
  ! Close data files.
  close (unit = sequencesfasta_unit)
  close (unit = numbers_unit)
  close (unit = fasta_unit)
  stop
end program removegaps
