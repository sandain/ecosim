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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! On September 16, 2006, I now invent a new strain from the
! parent population when a new ecotype is created.  This new
! strain is the parent founder of the new ecotype.
! These are the programs that were changed between Sept. 21, 2005,
! and Sept. 16, 2006:
! donicheinv, choosepopnicheinv, coalescence, startpops
! Also, some of the CALL statements in the main programs were 
! changed.
! This is a non-ultrametric divergence design, with drift but no
! recombination. September 21, 2005
! note, this will soon be replaced by a quicker method that Danny
! is writing, which will not require the n-square binning algorithm,
! but rather an nlogn sorting routine
! note, the programs I've changed for nonultra are:
! main, whicheventandwhen, 
! divergenonultra (and deleted binning), binningdanny, coalescence,
! and startpops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program bruteforce
  use methods
  implicit none
  ! Local variables
  character(len = 256)                :: inputFile
  character(len = 256)                :: outputFile
  character(len = 256)                :: timeFile
  character(len = *), parameter       :: outputFormat = &
    "(f15.7,',',f15.7,',',i5,',',1e9.2,6(',',1e18.10))"
  logical                             :: fileExists
  integer                             :: i
  integer                             :: indexxn
  integer                             :: indexsigma
  integer                             :: indexomega
  integer                             :: indexnpop
  integer                             :: lengthseq
  integer                             :: npop
  integer                             :: npoprange(2)
  integer                             :: nrep
  integer                             :: nu
  integer                             :: numcrit
  integer                             :: numincsomega
  integer                             :: numincssigma
  integer                             :: numincsnpop
  integer                             :: numincsxn
  integer                             :: whichavg
  integer                             :: realdata(1000)
  integer, parameter                  :: outputUnit = 2
  integer, parameter                  :: timeUnit = 4
  double precision                    :: omega
  double precision                    :: omegarange(2)
  double precision                    :: sigma  
  double precision                    :: sigmarange(2)
  double precision                    :: xn
  double precision                    :: xnrange(2)
  double precision                    :: diffomega
  double precision                    :: diffsigma
  double precision                    :: diffnpop
  double precision                    :: diffxn
  real                                :: time1
  real                                :: time2
  real                                :: delta
  real                                :: probthreshold
  real                                :: avgsuccess(6)
  real                                :: crit(1000)
  ! Provide default file names to use.
  inputFile = 'bruteforceIn.dat'
  outputFile = 'bruteforceOut.dat'
  timeFile = 'time.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, inputFile)
      case (2)
        call getArgument (2, outputFile)
      case (3)
        call getArgument (3, timeFile)
      case (4)
        call getArgument (4, numberThreads)
      case (5)
        call getArgument (5, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) "Expected: bruteforceIn.dat bruteforceOut.dat time.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that the input file exists.
  inquire (file = trim (inputFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The bruteforceIn.dat file was not found at: ", &
      trim (inputFile)
    ! Error, exit the program.
    stop
  end if
  ! Read in the input file.
  call readinput (trim (inputFile), omegarange, sigmarange, npoprange, &
    xnrange, numincsomega, numincssigma, numincsnpop, numincsxn, &
    numcrit, nu, nrep, lengthseq, realdata, crit, whichavg, probthreshold)
  ! Open the output files.
  open (unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  open (unit = timeUnit, file = trim (timeFile), access = 'sequential', form = 'formatted')
  call cpu_time (time1)
  diffomega = log10 (omegarange(2)) - log10 (omegarange(1))
  diffsigma = log10 (sigmarange(2)) - log10 (sigmarange(1))
  diffnpop = log10 (real (npoprange(2))) - log10 (real (npoprange(1)))
  diffxn = log10 (xnrange(2)) - log10 (xnrange(1))
  if (numincsomega .eq. 0) then
    omegarange(2) = omegarange(1) + 1.0
    diffomega = 2.0 * (log10 (omegarange(2)) - log10 (omegarange(1)))
    numincsomega = 1
  endif
  if (numincssigma .eq. 0) then
    sigmarange(2) = sigmarange(1) + 1.0
    diffsigma = 2.0 * (log10 (sigmarange(2)) - log10 (sigmarange(1)))
    numincssigma = 1
  endif
  if (numincsnpop .eq. 0) then
    npoprange(2) = npoprange(1) + 1
    diffnpop = 2.0 * (log10 (real (npoprange(2))) - log10 (real (npoprange(1))))
    numincsnpop = 1
  endif
  if (numincsxn .eq. 0) then
    xnrange(2) = xnrange(1) + 1.0
    diffxn = 2.0 * (log10 (xnrange(2)) - log10 (xnrange(1)))
    numincsxn = 1
  endif
  indexomega = 0
  do while (indexomega .lt. numincsomega)
    omega = 10.0 ** (log10 (omegarange(1)) + indexomega * (diffomega / numincsomega) * 0.999)
    indexsigma = 0
    do while (indexsigma .lt. numincssigma)
      sigma = 10.0 ** (log10 (sigmarange(1)) + indexsigma * (diffsigma / numincssigma) * 0.999)
      indexnpop = 0
      do while (indexnpop .lt. numincsnpop)
        npop = nint (10.0 ** (log10 (real (npoprange(1))) + indexnpop * (diffnpop / numincsnpop) * 0.999))
        if (npop .gt. nu) npop = nu
        indexxn = 0
        do while (indexxn .lt. numincsxn)
          xn = 10.0 ** (log10 (xnrange(1)) + indexxn * (diffxn / numincsxn) * 0.999)
          if (debug) then
            write (unit = *, fmt = *) 'omega= ', omega
            write (unit = *, fmt = *) 'sigma= ', sigma
            write (unit = *, fmt = *) 'npop= ', npop
            write (unit = *, fmt = *) 'xn= ', xn
          end if
          ! avgsuccess is a count of the number of results that are within
          ! X% tolerance for a particular set of parameter values.
          call runProgram (omega, sigma, npop, xn, numcrit, nu, nrep, &
            lengthseq, realdata, crit, avgsuccess)
          if (avgsuccess(1) .gt. 0.0) then
            write (unit = outputUnit, fmt = outputFormat) &
              omega, sigma, npop, xn, (avgsuccess(i), i = 1, 6)
          end if
          if (debug) then
            write (unit = *, fmt = *) 'avgsuccess= ', (avgsuccess(i), i = 1, 6)
            write (unit = *, fmt = *)
          end if
          indexxn = indexxn + 1
        end do
        indexnpop = indexnpop + 1
      end do
      indexsigma = indexsigma + 1
    end do
    indexomega = indexomega + 1
  end do
  ! delta is the time elapsed
  call cpu_time (time2)
  delta = time2 - time1
  write (unit = timeUnit, fmt = *) delta, '    time'
  ! Close the random number generator.
  call randomClose ()
  ! Close the output files.
  close (unit = outputUnit)
  close (unit = timeUnit)
  stop
end program bruteforce
