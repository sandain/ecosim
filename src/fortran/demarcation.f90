!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation and periodic
!    selection, yielding a certain number of ecotypes.
!
!    Copyright (C) 2009-2018  Fred Cohan, Wesleyan University
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
!> The demarcation program calculates the most likely value of npop (number of
!> ecotypes) using the Nelder-Mead simplex method for function minimization.
!> This npop value is then used to demarcate ecotypes.
!>
!> @pre  Requires that the demarcationIn.dat file exist and that it contains
!>         all of the variables necessary to calculate the most likely value
!>         of npop.
!>
!> @post Generates the demarcationOut.dat file that contains the most likely
!>         value of npop and its likelihood value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program demarcation
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
  integer(kind = int32)            :: bestnpop
  integer(kind = int32)            :: testednpop
  integer(kind = int32)            :: i
  integer(kind = int32)            :: istep
  integer(kind = int32), parameter :: outputUnit = 4
  real(kind = real64)              :: ratio
  real(kind = real64)              :: likelihoodone
  real(kind = real64)              :: bestlikelihood
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
        write (unit = *, fmt = *) &
          "Expected: demarcationIn.dat demarcationOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) &
      "The demarcationIn.dat file was not found at: ", trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read the input file.
  call readinput (trim (inputFile), omega, sigma, npop, istep)
  ! Open the output file.
  open (unit = outputUnit, file = trim (outputFile), &
    access = 'sequential', form = 'formatted')
  ! Start off with the best npop value equal to the predicted value.
  bestnpop = npop
  bestlikelihood = 0.0d0
  likelihoodone = 0.0d0
  ! Make sure omega and sigma are greater than zero.
  if (omega .gt. 1.0d-6 .and. sigma .gt. 1.0d-6) then
    ! Test npop value = 1.
    if (debug) then
      write (unit = *, fmt = *) 'omega= ', omega
      write (unit = *, fmt = *) 'sigma= ', sigma
      write (unit = *, fmt = *) 'npop= ', 1
    end if
    call runFredProgram (omega, sigma, 1, numcrit, nu, nrep, lengthseq, &
      realdata, crit, avgsuccess)
    likelihoodone = avgsuccess(jwhichxavg)
    if (debug) then
      write (unit = *, fmt = *) 'yvalue= ', likelihoodone
    end if
    if (likelihoodone .gt. 1.0d-6) then
      bestnpop = 1
      bestlikelihood = likelihoodone
      ! Test npop values from istep + 1 to the npop estimate.
      do testednpop = istep + 1, npop, istep
        if (debug) then
          write (unit = *, fmt = *) 'omega= ', omega
          write (unit = *, fmt = *) 'sigma= ', sigma
          write (unit = *, fmt = *) 'npop= ', testednpop
        end if
        call runFredProgram (omega, sigma, testednpop, numcrit, nu, nrep, &
          lengthseq, realdata, crit, avgsuccess)
        if (debug) then
          write (unit = *, fmt = *) 'yvalue= ', avgsuccess(jwhichxavg)
        end if
        if (avgsuccess(jwhichxavg) .lt. 1.0d-6) cycle
        ! Do the likelihood ratio test.
        ratio = -2.0 * log (bestlikelihood / avgsuccess(jwhichxavg))
        if (avgsuccess(jwhichxavg) .gt. bestlikelihood .and. ratio .gt. 3.84) then
!        if (avgsuccess(jwhichxavg) .gt. bestlikelihood .and. ratio .gt. 6.83) then
          bestnpop = testednpop
          bestlikelihood = avgsuccess(jwhichxavg)
        end if
      end do
    end if
  end if
  ! Output the answer.
  write (unit = outputUnit, fmt = *) &
    'npop ', 1, ' likelihood ', likelihoodone
  write (unit = outputUnit, fmt = *) &
    'npop ', bestnpop, ' likelihood ', bestlikelihood
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
  !> @param[out]    omega         The omega value to be tested.
  !> @param[out]    sigma         The sigma value to be tested.
  !> @param[out]    npop          The npop value to be tested.
  !> @param[out]    istep         The factor by which we tweak npop.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readinput (fname, omega, sigma, npop, istep)
    character(len = *), intent(in)       :: fname
    integer(kind = int32), intent(out)   :: npop
    integer(kind = int32), intent(out)   :: istep
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
    read (unit = input_unit, fmt = *) istep
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

end program demarcation
