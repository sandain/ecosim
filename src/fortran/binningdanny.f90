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
!> The binningdanny program performs a complete linkage hierarchical
!> clustering based on sequence identity.
!>
!> @pre  Requires that the divergencematrix.dat and binlevels.dat files exist
!>         and that they contain all of the variables necessary to perform
!>         the binning.
!>
!> @post Generates the binningdannyOut.dat file that contains the number of
!>         clusters at each identity level.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program binningdanny
  use methods
  implicit none
  character(len = 256)          :: divergencematrixFile
  character(len = 256)          :: binlevelsFile
  character(len = 256)          :: outputFile
  character(len = *), parameter :: outputFormat = "(F12.8, I16)"
  integer                       :: allocateStatus
  integer                       :: i
  integer                       :: j
  integer                       :: k
  integer                       :: l
  integer                       :: xi
  integer                       :: xj
  integer                       :: n_seq
  integer                       :: n_clusters
  integer                       :: length_seq
  integer                       :: n_level
  integer, allocatable          :: cluster(:,:)
  integer, allocatable          :: cluster_size(:)
  integer, parameter            :: identitymatrixUnit = 1
  integer, parameter            :: outputUnit = 2
  integer, parameter            :: binlevelsUnit = 3
  logical                       :: fileExists
  real                          :: identity_level
  real                          :: x
  real, allocatable             :: seq_identity(:,:)
  real, allocatable             :: cluster_dist(:,:)
  real, allocatable             :: bin_id_level(:)
  ! Provide default file names to use.
  divergencematrixFile = 'divergencematrix.dat'
  binlevelsFile = 'binlevels.dat'
  outputFile = 'binningdannyOut.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, divergencematrixFile)
      case (2)
        call getArgument (2, binlevelsFile)
      case (3)
        call getArgument (3, outputFile)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: divergencematrix.dat binlevels.dat binningdannyOut.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input files exist.
  inquire (file = trim (divergencematrixFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) &
      "The divergencematrix.dat file was not found at: ", &
      trim (divergencematrixFile)
    ! Error, exit the program.
    stop
  end if
  inquire (file = trim (binlevelsFile), exist = fileExists)
  if (fileExists .neqv. .true.) then
    write (unit = *, fmt = *) "The binlevels.dat file was not found at: ", &
      trim (binlevelsFile)
    ! Error, exit the program.
    stop
  end if
  ! Open Input files
  open (unit = identitymatrixUnit, file = trim (divergencematrixFile), &
    access = "sequential", status = "unknown")
  open (unit = binlevelsUnit, file = trim (binlevelsFile), &
    access = "sequential", status = "unknown")
  ! Open Output file
  open (unit = outputUnit, file = trim (outputFile), &
    access = "sequential", status = "unknown")
  ! Read in divergencematrixFile
  read (unit = identitymatrixUnit, fmt = *) n_seq, length_seq
  ! Allocate variables.
  allocate (seq_identity(n_seq, n_seq), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: seq_identity!"
    ! Error, exit the program.
    stop
  end if
  allocate (cluster(n_seq, n_seq), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: cluster!"
    ! Error, exit the program.
    stop
  end if
  allocate (cluster_dist(n_seq, n_seq), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: cluster_dist!"
    ! Error, exit the program.
    stop
  end if
  allocate (cluster_size(n_seq), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: cluster_size!"
    ! Error, exit the program.
    stop
  end if
  do i = 1, n_seq
    read (unit = identitymatrixUnit, fmt = *) &
      (seq_identity(i, j), j = 1, n_seq)
  end do
  ! Read in binlevels file.
  read (unit = binlevelsUnit, fmt = *) n_level
  allocate (bin_id_level(n_level), stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: bin_id_level!"
    ! Error, exit the program.
    stop
  end if
  do i = 1, n_level
    read (unit = binlevelsUnit, fmt = *) bin_id_level(i)
  end do
  ! Loop: for each identity level find number of bins
  do l = 1, n_level
    ! Initialize variables
    ! identity_level = max identity difference in cluster
    ! cluster(i,j) = jth element of ith cluster
    ! cluster_size(i) = size of ith cluster
    ! cluster_dist(i,j) = min percent identity seq in cluster i vs j
    ! n_clusters = number of clusters
    identity_level = bin_id_level(l)
    do i = 1, n_seq
      do j = 1, n_seq
        cluster(i, j) = 0
        cluster_dist(i, j) = seq_identity(i, j)
      end do
      cluster(i, 1) = i
      cluster_size(i) = 1
    end do
    n_clusters = n_seq
    ! Loop: find closest pair of clusters, combine, recompute distances
    do
      ! Find closest pair, x = distance between pair xi,xj; xi < xj
      x = -1.0
      xi = 1
      xj = 2
      do i = 1, n_seq
        do j = i + 1, n_seq
          if (cluster_dist(i, j) .lt. x) cycle
          x = cluster_dist(i, j)
          xi = i
          xj = j
        end do
      end do
      ! If closest pair further apart than identity_level then done
      if (x .lt. identity_level) exit
      ! Otherwise, combine the closest pair into xi and eliminate xj
      do i = cluster_size(xi) + 1, cluster_size(xi) + cluster_size(xj)
        cluster(xi, i) = cluster(xj, i - cluster_size(xi))
      end do
      cluster_size(xi) = cluster_size(xi) + cluster_size(xj)
      n_clusters = n_clusters - 1
      do i = 1, n_seq
        cluster_dist(xj, i) = -1.0
        cluster_dist(i, xj) = -1.0
        cluster(xj, i) = 0
      end do
      cluster_size(xj) = 0
      ! Recalculate distance from xi to all other active clusters and repeat
      do i = 1, n_seq
        if (cluster_size(i) .eq. 0) cycle
        x = 1.0
        do j = 1, cluster_size(xi)
          do k = 1, cluster_size(i)
            if (seq_identity(cluster(xi, j), cluster(i, k)) .lt. x) then
              x = seq_identity(cluster(xi, j), cluster(i, k))
            end if
          end do
        end do
        cluster_dist(xi, i) = x
        cluster_dist(i, xi) = x
      end do
    end do
    ! Output to the output file
    write (unit = outputUnit, fmt = outputFormat) identity_level, n_clusters
  end do
  ! Deallocate memory.
  deallocate (seq_identity, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: seq_identity!"
    ! Error, exit the program.
    stop
  end if
  deallocate (cluster, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: cluster!"
    ! Error, exit the program.
    stop
  end if
  deallocate (cluster_dist, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: cluster_dist!"
    ! Error, exit the program.
    stop
  end if
  deallocate (cluster_size, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: cluster_size!"
    ! Error, exit the program.
    stop
  end if
  deallocate (bin_id_level, stat = allocateStatus)
  if (allocateStatus .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: bin_id_level!"
    ! Error, exit the program.
    stop
  end if
  ! Close data files.
  close (unit = identitymatrixUnit)
  close (unit = binlevelsUnit)
  close (unit = outputUnit)
  ! Successful termination of program.
  stop
end program binningdanny
