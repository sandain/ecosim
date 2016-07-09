!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation and periodic
!    selection, yielding a certain number of ecotypes.
!
!    Copyright (C) 2009-2015  Fred Cohan, Wesleyan University
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
!> The fredmethod program uses the Nelder-Mead simplex method for function
!> minimization to optimize the solution provided by the brute force search
!> method.
!>
!> @pre  Requires that the fredmethodIn.dat file exists and that it contains
!>         all of the variables necessary to run the fredmethod program.
!>
!> @post Generates the fredmethodOut.dat file that contains the output from
!>         the fredmethod.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program fredmethod
  ! Load intrinsic modules.
  use, intrinsic :: iso_fortran_env
  ! Load our modules.
  use :: methods

  implicit none

  ! Local variables
  character(len = 256)             :: inputFile
  character(len = 256)             :: outputFile
  logical                          :: fileExists
  integer(kind = int32)            :: npop
  integer(kind = int32)            :: i
  integer(kind = int32), parameter :: outputUnit = 4
  real(kind = real64)              :: omega
  real(kind = real64)              :: sigma
  real(kind = real64)              :: avgsuccess(6)
  ! bldanny common block
  integer(kind = int32) :: numcrit
  integer(kind = int32) :: nu
  integer(kind = int32) :: nrep
  integer(kind = int32) :: lengthseq
  integer(kind = int32) :: realdata(1000)
  integer(kind = int32) :: jwhichxavg
  real(kind = real32)   :: crit(1000)
  common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
  ! Provide default file names to use.
  inputFile = 'fredmethodIn.dat'
  outputFile = 'fredmethodOut.dat'
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
          "Expected: fredmethodIn.dat fredmethodOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) &
      "The fredmethod input file was not found at: ", trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read the input file.
  call readinput (trim (inputFile), npop, omega, sigma)
  ! Open the output file.
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! Make sure omega and sigma are greater than zero.
  if (omega .gt. 0.0 .and. sigma .gt. 0.0) then
    if (debug) then
      write (unit = *, fmt = *) 'npop= ', npop
      write (unit = *, fmt = *) 'omega= ', omega
      write (unit = *, fmt = *) 'sigma= ', sigma
    end if
    call runFredProgram (omega, sigma, npop, numcrit, nu, nrep, &
      lengthseq, realdata, crit, avgsuccess)
    if (debug) then
      write (unit = *, fmt = *) 'yvalue= ', avgsuccess(jwhichxavg)
    end if
  end if
  ! Output the answer.
  write (unit = outputUnit, fmt = *) &
    omega, sigma, npop, avgsuccess(jwhichxavg)
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
  !> @param[out]    npop          The npop value to be tested.
  !> @param[out]    omega         The omega value to be tested.
  !> @param[out]    sigma         The sigma value to be tested.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readinput (fname, npop, omega, sigma)
    character(len = *), intent(in)       :: fname
    integer(kind = int32), intent(out)   :: npop
    real(kind = real64), intent(out)     :: omega
    real(kind = real64), intent(out)     :: sigma
    ! Local variables
    integer(kind = int32)            :: iii
    integer(kind = int32)            :: jcrit
    integer(kind = int32), parameter :: input_unit = 1
    ! bldanny common block
    integer(kind = int32) :: numcrit
    integer(kind = int32) :: nu
    integer(kind = int32) :: nrep
    integer(kind = int32) :: lengthseq
    integer(kind = int32) :: realdata(1000)
    integer(kind = int32) :: jwhichxavg
    real(kind = real32)   :: crit(1000)
    common/bldanny/numcrit,nu,nrep,lengthseq,realdata,crit,jwhichxavg
    ! Open the input file.
    open (unit = input_unit, file = fname, action = 'read', &
      access = 'sequential', form = 'formatted')
    ! numcrit is the number of criteria for making cluster bins
    read (unit = input_unit, fmt = *) numcrit
    do jcrit = 1, numcrit
      ! crit() is the criterion value
      ! realdata() is the number of bins in the real data
      read (unit = input_unit, fmt = *) crit(jcrit), realdata(jcrit)
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

end program fredmethod
