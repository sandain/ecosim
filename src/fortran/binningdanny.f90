!    Ecotype Simulation models the sequence diversity within a bacterial clade as
!    the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
! 
!    Copyright (C) 2009  Fred Cohan, Wesleyan University
!                        Danny Krizanc, Wesleyan University
!                        Jason M. Wood, Montana State University
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complete linkage hierarchical clustering based sequence identity
!   Input:
!     n_seq - number of sequences (upto 2000)
!     identity_level - percent identity of bin size sought (between 0.0 and 1.0)
!     seq_identity - matrix of sequence percent identity (between 0.0 and 1.0)
!     n_level - number of bin levels (upto 100)
!     bin_id_level - list of bin identity levels (betweem 0.0 and 1.0)
!   Output:
!     n_clusters - number of clusters (bins) at identity level sought
!     clusters - matrix of clusters
!     cluster_size - vector of cluster sizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program binningdanny
  use methods
  implicit none
  character(len = 256) :: identitymatrix_dat
  character(len = 256) :: binlevels_dat
  character(len = 256) :: output_dat
  integer              :: allocate_status
  integer              :: i
  integer              :: j
  integer              :: k
  integer              :: l
  integer              :: xi
  integer              :: xj
  integer              :: n_seq
  integer              :: n_clusters
  integer              :: length_seq
  integer              :: n_level
  integer, allocatable :: cluster(:,:)
  integer, allocatable :: cluster_size(:)
  integer, parameter   :: identitymatrix_unit = 1
  integer, parameter   :: output_unit = 2
  integer, parameter   :: binlevels_unit = 3
  real                 :: identity_level
  real                 :: x
  real, allocatable    :: seq_identity(:,:)
  real, allocatable    :: cluster_dist(:,:)
  real, allocatable    :: bin_id_level(:)
  ! Read in command line arguments
  ! Expect path and filename for: identitymatrix.dat binlevels.dat output.dat
  if (iargc() .ge. 3) then
    call getArgument (1, identitymatrix_dat)
    call getArgument (2, binlevels_dat)
    call getArgument (3, output_dat)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: identitymatrix.dat binlevels.dat output.dat"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 4) then
    call getArgument (4, debug)
  else
    debug = .false.
  end if
  ! Open Input files
  open (unit = identitymatrix_unit, file = trim (identitymatrix_dat), access = "sequential", status = "unknown")
  open (unit = binlevels_unit, file = trim (binlevels_dat), access = "sequential", status = "unknown")
  ! Open Output file
  open (unit = output_unit, file = trim (output_dat), access = "sequential", status = "unknown")
  ! Read in identitymatrix_dat
  read (unit = identitymatrix_unit, fmt = *) n_seq, length_seq
  ! Allocate variables.
  allocate (seq_identity(n_seq, n_seq), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (cluster(n_seq, n_seq), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (cluster_dist(n_seq, n_seq), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (cluster_size(n_seq), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  do i = 1, n_seq
    read (unit = identitymatrix_unit, fmt = *) (seq_identity(i, j), j = 1, n_seq)
  end do
  ! Read in binlevels_dat
  read (unit = binlevels_unit, fmt = *) n_level
  allocate (bin_id_level(n_level), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  do i = 1, n_level
    read (unit = binlevels_unit, fmt = *) bin_id_level(i)
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
        cluster(xi,i) = cluster(xj, i - cluster_size(xi))
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
            if (seq_identity(cluster(xi, j),cluster(i, k)) .lt. x) then
              x = seq_identity(cluster(xi, j),cluster(i, k))
            end if
          end do
        end do
        cluster_dist(xi,i) = x
        cluster_dist(i,xi) = x
      end do
    end do
    ! Output to output_dat
    write (unit = output_unit, fmt = *) identity_level, n_clusters
  end do
  ! Close data files.
  close (unit = identitymatrix_unit)
  close (unit = binlevels_unit)
  close (unit = output_unit)
  stop
end program binningdanny
