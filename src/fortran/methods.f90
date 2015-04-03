!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
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

  ! Acinas data type.
  type acinas_data
    integer              :: numcrit
    integer              :: nu
    integer              :: nrep
    integer              :: lengthseq
    integer              :: whichavg
    integer, allocatable :: realdata(:)
    real                 :: probthreshold
    real, allocatable    :: crit(:)
  end type acinas_data

  ! Number of increments data type.
  type number_increments
    integer :: npop
    integer :: omega
    integer :: sigma
    integer :: xn
  end type number_increments

  ! Parameters data type.
  type parameters_data
    integer          :: npop
    double precision :: omega
    double precision :: sigma
    double precision :: xn
  end type parameters_data

  ! Declare public methods.
  public :: diverge
  public :: getArgument
  public :: randomClose
  public :: randomInitialize
  public :: randomNumber
  public :: readinput
  public :: simulation

  ! Declare public types.
  public :: acinas_data
  public :: number_increments
  public :: parameters_data

  ! Declare public variables.
  logical, public :: debug = .false.       !< Display debug information.
  integer, public :: numberThreads = 1     !< The number of threads to start.

  ! Declare private global parameters.
  integer, parameter :: &
    EVENT_NICHE_INVASION     = 1001        !< Niche invasion event.
  integer, parameter :: &
    EVENT_PERIODIC_SELECTION = 1002        !< Periodic selection event.
  integer, parameter :: &
    EVENT_POPULATION_DRIFT   = 1003        !< Population drift event.

  type(ziggurat_t), allocatable :: rng(:)  !< The state variables for the RNG.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets an argument from the command line.
  !
  !  Params:
  !    arg : Argument number to get.
  !    out : Value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface getArgument
    module procedure getLogicalArgument
    module procedure getIntegerArgument
    module procedure getRealArgument
    module procedure getStringArgument
  end interface getArgument

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Complete linkage hierarchical clustering based sequence identity.
  !
  !  Params:
  !    numcrit            : Number of bin levels, input.
  !    crit               : List of bin identity levels (between 0.0 and 1.0), input.
  !    nu                 : Number of sequences, input.
  !    nclustersarray     : Number of clusters (bins) at identity level sought, output.
  !    identitymatrix     : Matrix of clusters, output.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine binningdanny (numcrit, crit, nu, nclustersarray, identitymatrix)
    integer, intent(in)               :: numcrit
    integer, intent(in)               :: nu
    integer, intent(inout)            :: nclustersarray(:)
    real, intent(in)                  :: crit(:)
    real, intent(in)                  :: identitymatrix(:,:)
    ! Local variables.
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: xi
    integer :: xj
    integer :: n_clusters
    integer :: cluster_size(nu)
    integer :: cluster(nu, nu)
    real    :: x
    real    :: identity_level
    real    :: cluster_dist(nu, nu)
    ! Loop: for each identity level find number of bins
    do l = 1, numcrit
      ! Initialize variables
      ! identity_level = max identity difference in cluster
      ! cluster(i,j) = jth element of ith cluster
      ! cluster_size(i) = size of ith cluster
      ! cluster_dist(i,j) = min percent identity seq in cluster i vs j
      ! n_clusters = number of clusters
      identity_level = crit(l)
      do i = 1, nu
        do j = 1, nu
          cluster(i, j) = 0
          cluster_dist(i, j) = identitymatrix(i, j)
        end do
        cluster(i, 1) = i
        cluster_size(i) = 1
      end do
      n_clusters = nu
      ! Loop: find closest pair of clusters, combine, recompute distances
      do
        ! Find closest pair, x = distance between pair xi,xj; xi < xj
        x = -1.0
        xi = 1
        xj = 2
        do i = 1, nu
          do j = i + 1, nu
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
        do i = 1, nu
          cluster_dist(xj, i) = -1.0
          cluster_dist(i, xj) = -1.0
          cluster(xj, i) = 0
        end do
        cluster_size(xj) = 0
        ! Recalculate distance from xi to all other active clusters and repeat
        do i = 1, nu
          if (cluster_size(i) .eq. 0) cycle
          x = 1.0
          do j = 1, cluster_size(xi)
            do k = 1, cluster_size(i)
              if (identitymatrix(cluster(xi, j), cluster(i, k)) .lt. x) then
                x = identitymatrix(cluster(xi, j), cluster(i, k))
              end if
            end do
          end do
          cluster_dist(xi, i) = x
          cluster_dist(i, xi) = x
        end do
      end do
      nclustersarray(l) = n_clusters
    end do
    return
  end subroutine binningdanny

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Canonical
  !
  !  Params:
  !    npop      :
  !    nu        :
  !    numstrain :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine canonical (npop, nu, numstrain)
    integer, intent(in)               :: npop
    integer, intent(in)               :: nu
    integer, intent(inout)            :: numstrain(:)
    ! Local variables.
    integer :: upto
    integer :: numberofstrains
    integer :: istrain
    integer :: ipop
    real    :: x
    numstrain(1) = nu
    if (npop .eq. 1) return
    ! for the upto-1 population, choose a random population (ipop)
    do upto = 2, npop
      do
        call randomNumber (x)
        ipop = int ((upto - 1) * x) + 1
        numberofstrains = numstrain(ipop)
        if (numberofstrains .ne. 1) exit
      end do
      ! next, pick a random breakpoint
      call randomNumber (x)
      istrain = int ((numberofstrains - 1) * x) + 1
      numstrain(ipop) = istrain
      numstrain(upto) = numberofstrains - istrain
    end do
    return
  end subroutine canonical

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Choosepopdrift
  !
  !  Params:
  !    popfordrift :
  !    activepop   :
  !    realname    :
  !    xn          :
  !    numstrain   :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine choosepopdrift (popfordrift, activepop, realname, xn, numstrain)
    integer, intent(out)              :: popfordrift
    integer, intent(in)               :: activepop
    integer, intent(in)               :: realname(:)
    integer, intent(in)               :: numstrain(:)
    double precision, intent(in)      :: xn
    ! Local variables.
    integer :: jpop
    real    :: x
    real    :: totaldrift
    real    :: probdrift(4000)
    totaldrift = 0.0
    do jpop = 1, activepop
      probdrift(realname(jpop)) = totaldrift + real (1.0 / xn) * numstrain(realname(jpop)) * (numstrain(realname(jpop)) - 1.0) * 0.5
      totaldrift = probdrift(realname(jpop))
    end do
    if (totaldrift .gt. 0.0) then
      do jpop = 1, activepop
        probdrift(realname(jpop)) = probdrift(realname(jpop)) / totaldrift
      end do
      call randomNumber (x)
      do jpop = 1, activepop
        if (x .lt. probdrift(realname(jpop))) then
          popfordrift = realname(jpop)
          return
        endif
      end do
    end if
    return
  end subroutine choosepopdrift

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Choosepopnicheinv
  !
  !  Params:
  !    nu                    :
  !    numstrain             :
  !    realname              :
  !    popforparent          :
  !    popfornascent         :
  !    strainparent          :
  !    activepop             :
  !    ntotalpopwithinvented :
  !    numanc                :
  !    numanctot             :
  !    ancdist               :
  !    numdesc               :
  !    div                   :
  !    nameanc               :
  !    namedesc              :
  !    namestrain            :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine choosepopnicheinv (numstrain, realname, popforparent, popfornascent, strainparent, activepop, &
    ntotalpopwithinvented, numanc, numanctot, ancdist, numdesc, div, nameanc, namedesc, namestrain)
    integer, intent(inout)            :: numstrain(:)
    integer, intent(in)               :: realname(:)
    integer, intent(out)              :: popforparent
    integer, intent(out)              :: popfornascent
    integer, intent(out)              :: strainparent
    integer, intent(in)               :: activepop
    integer, intent(out)              :: ntotalpopwithinvented
    integer, intent(inout)            :: numanc(:)
    integer, intent(out)              :: numanctot
    integer, intent(inout)            :: numdesc(:)
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    integer, intent(inout)            :: namestrain(:,:)
    real, intent(inout)               :: ancdist(:,:)
    real, intent(inout)               :: div(:)
    ! Local variables.
    integer :: chosen
    integer :: numeligible
    integer :: jpop
    integer :: eligibleparentpop(4000)
    integer :: eligiblenascentpop(4000)
    real    :: x
    ! first, choose the population for the nascent population
    numeligible = 0
    do jpop = 1, activepop
      numeligible = numeligible + 1
      eligiblenascentpop(numeligible) = realname(jpop)
    end do
    call randomNumber (x)
    chosen = int (x * numeligible) + 1
    popfornascent = eligiblenascentpop(chosen)
    ! above, we use the first strain of the chosen population as the strain
    ! next, choose the population for the parent population
    numeligible = 0
    do jpop = 1, activepop
      if (realname(jpop) .eq. popfornascent) cycle
      numeligible = numeligible + 1
      eligibleparentpop(numeligible) = realname(jpop)
    end do
    call randomNumber (x)
    chosen = int (x * numeligible) + 1
    popforparent = eligibleparentpop(chosen)
    ! now we create a parent strain from population POPFORPARENT
    numstrain(popforparent) = numstrain(popforparent) + 1
    chosen = numstrain(popforparent)
    ntotalpopwithinvented = ntotalpopwithinvented + 1
    strainparent = ntotalpopwithinvented
    namestrain(popforparent, chosen) = ntotalpopwithinvented
    numanctot = numanctot + 1
    numanc(strainparent) = 1
    nameanc(strainparent, 1) = numanctot
    ancdist(strainparent, 1) = 0
    numdesc(numanctot) = 1
    namedesc(numanctot, 1) = strainparent
    div(numanctot) = 0
    return
  end subroutine choosepopnicheinv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Choosepopps
  !
  !  Params:
  !    numstrain :
  !    realname  :
  !    popforps  :
  !    activepop :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine choosepopps (numstrain, realname, popforps, activepop)
    integer, intent(in)               :: numstrain(:)
    integer, intent(in)               :: realname(:)
    integer, intent(out)              :: popforps
    integer, intent(in)               :: activepop
    ! Local variables.
    integer :: chosen
    integer :: eligiblepspop(activepop)
    integer :: jpop
    integer :: numeligible
    real    :: x
    numeligible = 0
    do jpop = 1, activepop
      if (numstrain(realname(jpop)) .le. 1) cycle
      numeligible = numeligible + 1
      eligiblepspop(numeligible) = realname(jpop)
    end do
    call randomNumber (x)
    chosen = int (x * numeligible) + 1
    popforps = eligiblepspop(chosen)
    return
  end subroutine choosepopps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Coalescence
  !
  !  Params:
  !    numcoalesce  :
  !    namecoalesce :
  !    numanc       :
  !    numanctot    :
  !    numdesc      :
  !    time         :
  !    lengthseq    :
  !    div          :
  !    nameanc      :
  !    namedesc     :
  !    ancdist      :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine coalescence (numcoalesce, namecoalesce, numanc, numanctot, numdesc, time, lengthseq, div, &
    nameanc, namedesc, ancdist)
    integer, intent(in)               :: numcoalesce
    integer, intent(in)               :: namecoalesce(:)
    integer, intent(inout)            :: numanc(:)
    integer, intent(out)              :: numanctot
    integer, intent(inout)            :: numdesc(:)
    integer, intent(in)               :: lengthseq
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    double precision, intent(in)      :: time
    real, intent(inout)               :: div(:)
    real, intent(inout)               :: ancdist(:,:)
    ! Local variables.
    integer :: strainname
    integer :: idivaddactual
    integer :: jcoal
    integer :: iancnow
    integer :: nameofdesc
    integer :: namelastancest
    integer :: js
    integer :: numofanc
    integer :: numancest
    integer :: jdesc
    real    :: divaddexpect
    real    :: divaddactual
    numanctot = numanctot + 1
    ! we add one to the total number of ancestors NUMANCTOT
    ! next, for the new ancestor being created (number numanctot),
    ! we must list all its descendants
    numdesc(numanctot) = 0
    ! numdesc(numanctot) is the number of descendants of the ancestor
    ! currently being described
    do jcoal = 1, numcoalesce
      strainname = namecoalesce(jcoal)
      numancest = numanc(strainname)
      namelastancest = nameanc(strainname, numancest)
      ! namelastancest is the name of strainname's most distant ancestor,
      ! until now
      do js = 1, numdesc(namelastancest)
        numdesc(numanctot) = numdesc(numanctot) + 1
        nameofdesc = namedesc(namelastancest, js)
        namedesc(numanctot, numdesc(numanctot)) = nameofdesc
        ! next, give nameofdesc its oldest ancestor
        numanc(nameofdesc) = numanc(nameofdesc) + 1
        nameanc(nameofdesc, numanc(nameofdesc)) = numanctot
      end do
    end do
    div(numanctot) = real (time / lengthseq)
    ! note do-loop below is an addition for nonultra
    do jcoal = 1, numcoalesce
      strainname = namecoalesce(jcoal)
       ! note, this penultimate ancestor is:
      iancnow = nameanc(strainname, numanc(strainname) - 1)
      divaddexpect = (div(numanctot) - div(iancnow))
      ! now get poisson--note we have to deal with integer number
      ! of mutations, so I multiply the per nucleotide divergence
      ! by lengthseq, and then divide later
      call poisson(divaddexpect * lengthseq, idivaddactual)
      divaddactual = real (idivaddactual) / lengthseq
      ! divaddactual is the actual number of substitutions
      ! in the node going to the last ancestor from a lineage,
      ! compared to the previous ancestor of strain strainname
      ! now we deal with all the descendants of ancestor iancnow
      do jdesc = 1, numdesc(iancnow)
        nameofdesc = namedesc(iancnow, jdesc)
        numofanc = numanc(nameofdesc)
        ancdist(nameofdesc, numofanc) = divaddactual + ancdist(nameofdesc, numofanc - 1)
      end do
    end do
    return
  end subroutine coalescence

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
    character(len = 1), intent(in)    :: seq(:,:)
    integer, intent(in)               :: numstrains
    integer, intent(in)               :: length
    type(tmatrixReal), intent(inout)  :: divergematrix
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
  !  Divergenonultra
  !
  !  Params:
  !    nu             :
  !    numanc         :
  !    nameanc        :
  !    ancdist        :
  !    identitymatrix :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine divergenonultra (nu, numanc, nameanc, ancdist, identitymatrix)
    integer, intent(in)               :: nu
    integer, intent(in)               :: numanc(:)
    integer, intent(in)               :: nameanc(:,:)
    real, intent(in)                  :: ancdist(:,:)
    real, intent(inout)               :: identitymatrix(:,:)
    ! Local variables.
    integer :: jstrain
    integer :: kstrain
    integer :: janc
    integer :: kanc
    real    :: xjukesdivergence
    real    :: divergence
    logical :: equal
    do jstrain = 1, nu
      do kstrain = 1, nu
        do janc = 1, numanc(jstrain)
          equal = .false.
          do kanc = 1, numanc(kstrain)
            if (nameanc(jstrain, janc) .eq. nameanc(kstrain, kanc)) then
              equal = .true.
              exit
            end if 
          end do
          if (equal) exit
        end do
        divergence = ancdist(jstrain, janc) + ancdist(kstrain, kanc)
        xjukesdivergence = -0.75 * (exp ((-4.0 / 3.0) * divergence) - 1.0)
        identitymatrix(jstrain, kstrain) = 1.0 - xjukesdivergence
      end do
    end do
    return
  end subroutine divergenonultra

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  doDrift
  !
  !  Params:
  !    popfordrift :
  !    ntotalpop   :
  !    numanc      :
  !    numanctot   :
  !    numdesc     :
  !    time        :
  !    lengthseq   :
  !    numstrain   :
  !    div         :
  !    nameanc     :
  !    namedesc    :
  !    namestrain  :
  !    ancdist     :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doDrift (popfordrift, ntotalpop, numanc, numanctot, numdesc, time, lengthseq, &
    numstrain, div, nameanc, namedesc, namestrain, ancdist)
    integer, intent(in)               :: popfordrift
    integer, intent(out)              :: ntotalpop
    integer, intent(inout)            :: numanc(:)
    integer, intent(out)              :: numanctot
    integer, intent(inout)            :: numdesc(:)
    integer, intent(in)               :: lengthseq
    integer, intent(inout)            :: numstrain(:)
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    integer, intent(inout)            :: namestrain(:,:)
    double precision, intent(in)      :: time
    real, intent(inout)               :: div(:)
    real, intent(inout)               :: ancdist(:,:)
    ! Local variables.
    integer :: driftee
    integer :: drifter
    integer :: numcoalesce
    integer :: namecoalesce(4000)
    real    :: x
   ! Pick two strains from population popfordrift
    call randomNumber (x)
    driftee = int (x * numstrain(popfordrift)) + 1
    do
      call randomNumber (x)
      drifter = int (x * numstrain(popfordrift)) + 1
      if (drifter .ne. driftee) exit
    end do
    numcoalesce = 2
    namecoalesce(1) = namestrain(popfordrift, drifter)
    namecoalesce(2) = namestrain(popfordrift, driftee)
    call coalescence (numcoalesce, namecoalesce, numanc, numanctot, numdesc, time, lengthseq, div, nameanc, namedesc, ancdist)
    namestrain(popfordrift, driftee) = namestrain(popfordrift, numstrain(popfordrift))
    numstrain(popfordrift) = numstrain(popfordrift) - 1
    ntotalpop = ntotalpop - 1
    return
  end subroutine doDrift

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  doNicheInvasion
  !
  !  Params:
  !    numanc        :
  !    numanctot     :
  !    numdesc       :
  !    time          :
  !    lengthseq     :
  !    activename    :
  !    activepop     :
  !    realname      :
  !    ntotalpop     :
  !    strainparent  :
  !    popfornascent :
  !    numstrain     :
  !    div           :
  !    nameanc       :
  !    namedesc      :
  !    namestrain    :
  !    ancdist       :
  !    numcoalesce   :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doNicheInvasion (numanc, numanctot, numdesc, time, lengthseq, activename, activepop, realname, &
    ntotalpop, strainparent, popfornascent, numstrain, div, nameanc, namedesc, namestrain, ancdist)
    integer, intent(inout)            :: numanc(:)
    integer, intent(out)              :: numanctot
    integer, intent(inout)            :: numdesc(:)
    integer, intent(in)               :: lengthseq
    integer, intent(inout)            :: activename(:)
    integer, intent(inout)            :: activepop
    integer, intent(inout)            :: realname(:)
    integer, intent(out)              :: ntotalpop
    integer, intent(in)               :: strainparent
    integer, intent(in)               :: popfornascent
    integer, intent(in)               :: numstrain(:)
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    integer, intent(in)               :: namestrain(:,:)
    double precision, intent(in)      :: time
    real, intent(inout)               :: div(:)
    real, intent(inout)               :: ancdist(:,:)
    ! Local variables.
    integer :: numcoalesce
    integer :: namecoalesce(4000)
    integer :: jpop
    integer :: jstrain
    ! Now do the coalescence of the whole nascent population (popfornascent) and strainparent
    numcoalesce = 0
    do jstrain = 1, numstrain(popfornascent)
      numcoalesce = numcoalesce + 1
      namecoalesce(numcoalesce)=namestrain(popfornascent,jstrain)
    end do
    numcoalesce = numcoalesce + 1
    namecoalesce(numcoalesce) = strainparent
    call coalescence (numcoalesce, namecoalesce, numanc, numanctot, numdesc, time, lengthseq, div, nameanc, namedesc, ancdist)
    ! Next, delete the nascent population
    realname(activename(popfornascent)) = realname(activepop)
    activename(realname(activepop)) = activename(popfornascent)
    activepop = activepop - 1
    ntotalpop = 0
    do jpop = 1, activepop
      ntotalpop = ntotalpop + numstrain(realname(jpop))
    end do
    return
  end subroutine doNicheInvasion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Dops
  !
  !  Params:
  !    popforps   :
  !    numanc     :
  !    numstrain  :
  !    numanctot  :
  !    numdesc    :
  !    time       :
  !    lengthseq  :
  !    activepop  :
  !    realname   :
  !    ntotalpop  :
  !    div        :
  !    nameanc    :
  !    namedesc   :
  !    namestrain :
  !    ancdist    :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doPeriodicSelection (popforps, numanc, numstrain, numanctot, numdesc, time, lengthseq, activepop, &
    realname, ntotalpop, div, nameanc, namedesc, namestrain, ancdist)
    integer, intent(in)               :: popforps
    integer, intent(inout)            :: numanc(:)
    integer, intent(inout)            :: numstrain(:)
    integer, intent(out)              :: numanctot
    integer, intent(inout)            :: numdesc(:)
    integer, intent(in)               :: lengthseq
    integer, intent(inout)            :: activepop
    integer, intent(in)               :: realname(:)
    integer, intent(inout)            :: ntotalpop
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    integer, intent(in)               :: namestrain(:,:)
    double precision, intent(in)      :: time
    real, intent(inout)               :: div(:)
    real, intent(inout)               :: ancdist(:,:)
    ! Local variables.
    integer :: numcoalesce
    integer :: namecoalesce(4000)
    integer :: jpop
    integer :: jstrain
    ! First, make an array of the strains that will be coalesced
    numcoalesce = 0
    do jstrain = 1, numstrain(popforps)
      numcoalesce = numcoalesce + 1
      namecoalesce(numcoalesce) = namestrain(popforps,jstrain)
    end do
    call coalescence (numcoalesce, namecoalesce, numanc, numanctot, numdesc, time, lengthseq, div, nameanc, namedesc, ancdist)
    numstrain(popforps) = 1
    ntotalpop = 0
    do jpop = 1, activepop
      ntotalpop = ntotalpop + numstrain(realname(jpop))
    end do
    return
  end subroutine doPeriodicSelection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Returns the number of populations eligible to be the nascent and parent ecotypes, as
  !  well as eligible for periodic selection for a population to be eligible to be a nascent
  !  ecotype, it can have any number of strain members. For a population to be eligible for
  !  periodic selection, it must have more than one strain member.
  !
  !  Params:
  !    numstrain       :
  !    activepop       :
  !    eligiblenascent :
  !    eligibleparent  :
  !    eligibleps      :
  !    realname        :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eligible (numstrain, activepop, eligiblenascent, eligibleparent, eligibleps, realname)
    integer, intent(in)               :: numstrain(:)
    integer, intent(in)               :: activepop
    integer, intent(out)              :: eligiblenascent
    integer, intent(out)              :: eligibleparent
    integer, intent(out)              :: eligibleps
    integer, intent(in)               :: realname(:)
    ! Local variables.
    integer :: jpop
    eligiblenascent = activepop
    eligibleparent = activepop
    if (activepop .eq. 1) then
      eligiblenascent = 0
      eligibleparent = 0
    endif
    eligibleps = 0
    do jpop = 1, activepop
      if (numstrain(realname(jpop)) .ne. 1) then
        eligibleps = eligibleps + 1
      end if
    end do
    return
  end subroutine eligible

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets an integer argument from the command line.
  !
  !  Params:
  !    arg : Argument number to get.
  !    out : Integer to return.
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
      write (unit = *, fmt = *) "ERROR: Integer number not supplied with argument ", arg
      stop
    end if
    return
  end subroutine getIntegerArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a logical argument from the command line.
  !
  !  Params:
  !    arg : Argument number to get.
  !    out : logical value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getLogicalArgument (arg, out)
    integer, intent(in)               :: arg
    logical, intent(out)              :: out
    ! Local variables.
    character(len = 100) :: buffer
    call getarg (arg, buffer)
    ! Unless we receive an appropriate true string, return false.
    out = .false.
    if (trim(buffer) .eq. 'true' .or. trim(buffer) .eq. 't') then
      out = .true.
    else if (trim(buffer) .eq. 'True' .or. trim(buffer) .eq. 'T') then
      out = .true.
    end if
    return
  end subroutine getLogicalArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a real argument from the command line.
  !
  !  Params:
  !    arg : Argument number to get.
  !    out : Real to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getRealArgument (arg, out)
    integer, intent(in)               :: arg
    real, intent(out)                 :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer              :: error
    call getarg (arg, buffer)
    read (unit = buffer, fmt = *, iostat = error) out
    if (error .ne. 0) then
      write (unit = *, fmt = *) "ERROR: Real number not supplied with argument ", arg
      stop
    end if
    return
  end subroutine getRealArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a string argument from the command line.
  !
  !  Params:
  !    arg : Argument number to get.
  !    out : String to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getStringArgument (arg, out)
    integer, intent(in)               :: arg
    character(len = *), intent(out)   :: out
    call getarg (arg, out)
    return
  end subroutine getStringArgument

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
    integer, intent(in)               :: iii
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
      call randomNumber (x)
      call ziggurat_seed (rng(i), int (iii * x))
    end do
  end subroutine randomInitialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generate a random number in the range [0,1).
  !>
  !> @param[out]    x             The random number generated
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine randomNumber (x)
    real, intent(out)                 :: x
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
  !  Takes the expected number of substitutions, and gives back the actual number of
  !  substitutions using the Poisson.
  !
  !  Params:
  !    xmean :
  !    nsubs :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson (xmean, nsubs)
    integer, intent(out)              :: nsubs
    real, intent(in)                  :: xmean
    ! Local variables.
    integer :: jmut
    real    :: x
    real    :: accumprob
    real    :: prob
    real    :: expect
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
  !  Read in the data for this simulation.
  !
  !  Params:
  !    fname   : File name to load.
  !    bottom  : Bottom parameter values.
  !    top     : Top parameter values.
  !    numincs : Number of increments for the parameters.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readinput (fname, acinas, bottom, top, numincs)
    character(len = *), intent(in)       :: fname
    type(acinas_data), intent(out)       :: acinas
    type(parameters_data), intent(out)   :: bottom
    type(parameters_data), intent(out)   :: top
    type(number_increments), intent(out) :: numincs
    ! Local variables.
    integer            :: allocate_status
    integer            :: iii
    integer            :: jcrit
    integer            :: io_status
    integer, parameter :: io_unit = 101
    ! Open the file
    open (unit = io_unit, file = fname, action = 'read', access = 'sequential', form = 'formatted')
    ! numcrit is the number of criteria for making cluster bins
    read (unit = io_unit, fmt = *) acinas%numcrit
    ! Allocate variables.
    allocate (acinas%crit(acinas%numcrit), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
    end if
    allocate (acinas%realdata(acinas%numcrit), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
    end if
    do jcrit = 1, acinas%numcrit
      ! realdata(jcrit) is the number of bins in the real data at criterion jcrit / 1000.0
      read (unit = io_unit, fmt = *) acinas%realdata(jcrit)
    end do
    ! realdata(jcrit) is the actual number of bins for the jcrit th criterion
    ! crit(jcrit) is the jcrit th criterion value
    do jcrit = 1, acinas%numcrit
      read (unit = io_unit, fmt = *) acinas%crit(jcrit)
    end do
    ! omega is the rate of niche invasion per eligible parental population
    ! measured as niche invasions per nucleotide substitution in a given gene
    ! omegabot and omegatop are the range
    read (unit = io_unit, fmt = *) bottom%omega, top%omega
    ! sigma is the rate of periodic selection per eligible population,
    ! measured as periodic selection events per population per nucleotide substitution in
    ! a given gene
    read (unit = io_unit, fmt = *) bottom%sigma, top%sigma
    ! npop is the number of ecotypes assumed to be in the environmental DNA sample
    read (unit = io_unit, fmt = *) bottom%npop, top%npop
    ! note that omega, sigma, and npop are the three parameters we will estimate by
    ! maximum likelihood
    ! xn is the actual population size in nature
    read (unit = io_unit, fmt = *) bottom%xn, top%xn
    ! numincs values are the numbers of increments to be investigated
    read (unit = io_unit, fmt = *) numincs%omega, numincs%sigma, numincs%npop, numincs%xn
    ! nu is the number of homologous gene sequences in the environmental sample
    ! following Acinas et al., this should be in the thousands.  The program is
    ! currently written to allow 10000 samples
    read (unit = io_unit, fmt = *) acinas%nu
    ! nrep is the number of replicate simulations for a given set of sigma, omega, and npop
    read (unit = io_unit, fmt = *) acinas%nrep
    ! iii is the odd random number seed (up to nine digits)
    read (unit = io_unit, fmt = *) iii
    call randomInitialize (iii)
    ! lengthseq is the length in nucleotides of the sequence analyzed
    read (unit = io_unit, fmt = *) acinas%lengthseq
    ! jwhichxavg indicates success level considered:
    !   1 represents 5.00 times
    !   2 represents 2.00 times
    !   3 represents 1.50 times
    !   4 represents 1.25 times
    !   5 represents 1.10 times
    !   6 represents 1.05 times
    read (unit = io_unit, fmt = *, iostat = io_status) acinas%whichavg
    if (io_status .ne. 0) then
      acinas%whichavg = 1
    end if
    ! get the number of individuals in each ecotype
    read (unit = io_unit, fmt = *, iostat = io_status) acinas%probthreshold
    if (io_status .ne. 0) then
      acinas%probthreshold = 1.0
    end if
    ! Close the file
    close (unit = io_unit)
    return
  end subroutine readinput

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Perform nrep repetitions of the fred method, and calculate the average number of results
  !  within (500%, 200%, 150%, 125%, 110%, and 105%) tolerance for a particular set of parameter values.
  !
  !  Params:
  !    acinas     : The acinas data for this simulation.
  !    parameters : The parameters for this simulation (omega, sigma, npop, xn).
  !    avgsuccess : The average success for this run.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation (acinas, parameters, avgsuccess)
    type(acinas_data), intent(in)     :: acinas
    type(parameters_data), intent(in) :: parameters
    real, intent(out)                 :: avgsuccess(6)
    ! Local variables.
    integer             :: i
    integer             :: irep
    logical             :: success(6, acinas%nrep)
    ! Initialize the avgsuccess array with values of 0.
    avgsuccess = (/ (0.0, i = 1, 6) /)
    ! Make sure that omega has valid values.
    if (isnan (parameters%omega)) return ! XXX Should use ieee_is_nan when supported.
    if (parameters%omega .lt. epsilon (0.0)) return
    if (parameters%omega .gt. huge (0.0) - 1.0) return
    ! Make sure that sigma has valid values.
    if (isnan (parameters%sigma)) return ! XXX Should use ieee_is_nan when supported.
    if (parameters%sigma .lt. epsilon (0.0)) return
    if (parameters%sigma .gt. huge (0.0) - 1.0) return
    ! Make sure that drift has valid values.
    if (isnan (parameters%xn)) return ! XXX Should use ieee_is_nan when supported.
    if (parameters%xn .lt. epsilon (0.0)) return
    if (parameters%xn .gt. huge (0.0) - 1.0) return
    ! Make sure that npop does not exceed the number of homologous gene
    ! sequences in the environmental sample, and is greater than zero.
    if (parameters%npop .gt. acinas%nu .or. parameters%npop .le. 0) return
#ifdef _OPENMP
    call omp_set_num_threads (numberThreads)
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(irep)
    !$OMP DO
#endif
    ! Run the simulation nrep number of times.
    do irep = 1, acinas%nrep
      call fredMethod (parameters%omega, parameters%sigma, parameters%xn, &
        parameters%npop, acinas%numcrit, acinas%nu, acinas%lengthseq, &
        acinas%realdata, acinas%crit, success(:, irep))
    end do
#ifdef _OPENMP
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
#endif
    ! Calculate the average success level for each tolerance range.
    do irep = 1, acinas%nrep
      do i = 1, 6
        if (success(i, irep)) avgsuccess(i) = avgsuccess(i) + 1.0
      end do
    end do
    avgsuccess = (/ (avgsuccess(i) / real (acinas%nrep), i = 1, 6) /)
    return
  end subroutine simulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Perform the fred method for the given values.
  !
  !  Params:
  !    omega      :
  !    sigma      :
  !    npop       :
  !    xn         :
  !    numcrit    :
  !    lengthseq  :
  !    realdata   :
  !    crit       :
  !    success    :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fredMethod (omega, sigma, xn, npop, numcrit, nu, lengthseq, &
    realdata, crit, success)
    double precision, intent(in)      :: omega
    double precision, intent(in)      :: sigma
    double precision, intent(in)      :: xn
    integer, intent(in)               :: npop
    integer, intent(in)               :: numcrit
    integer, intent(in)               :: nu
    integer, intent(in)               :: lengthseq
    integer, intent(in)               :: realdata(1000)
    real, intent(in)                  :: crit(1000)
    logical, intent(out)              :: success(6)
    ! Local variables.
    integer          :: event
    integer          :: activepop
    integer          :: ntotalpop
    integer          :: ntotalpopwithinvented
    integer          :: numanctot
    integer          :: popfordrift
    integer          :: popforparent
    integer          :: popfornascent
    integer          :: popforps
    integer          :: strainparent
    integer          :: eligiblenascent
    integer          :: eligibleparent
    integer          :: eligibleps
    integer          :: nclustersarray(numcrit)
    integer          :: numstrain(npop)
    integer          :: realname(npop)
    integer          :: nameanc(4000, 4000)        ! bl1 common block
    integer          :: namedesc(4000, 4000)       ! bl1 common block
    integer          :: namestrain(4000, 4000)     ! bl1 common block
    integer          :: numanc(4000)
    integer          :: numdesc(4000)
    integer          :: activename(4000)
    double precision :: time
    real             :: ancdist(2000, 2000)        ! bl10 common block
    real             :: div(4000)
    real             :: identitymatrix(nu, nu)     ! bl10 common block
    ! Start the simulation.
    call startpops(npop, nu, activepop, ntotalpop, numstrain, realname, activename, numanc, numdesc, &
      numanctot, time, ntotalpopwithinvented, nameanc, namedesc, namestrain, ancdist)
    do
      ! subroutine eligible returns the number of populations eligible to be the nascent and
      ! parent ecotypes, as well as eligible for periodic selection
      call eligible (numstrain, activepop, eligiblenascent, eligibleparent, eligibleps, realname)
      call whichEventAndWhen (eligibleparent, eligibleps, omega, sigma, xn, time, event, activepop, realname, numstrain)
      ! subroutine WhichEventAndWhen returns the name of the key event
      ! 'nicheinv' or 'periodic' as the value of EVENT, and gives the
      ! waiting time until the event (since the previous key event), TIME
      if (event .eq. EVENT_PERIODIC_SELECTION) then
        call choosepopps (numstrain, realname, popforps, activepop)
        ! subroutine choosepopps chooses the population for periodic selection;
        ! it is the population numbered popforps
        call doPeriodicSelection (popforps, numanc, numstrain, numanctot, numdesc, time, lengthseq, activepop, &
          realname, ntotalpop, div, nameanc, namedesc, namestrain, ancdist)
      endif
      if (event .eq. EVENT_NICHE_INVASION) then
        call choosepopnicheinv (numstrain, realname, popforparent, popfornascent, strainparent, activepop, &
          ntotalpopwithinvented, numanc, numanctot, ancdist, numdesc, div, nameanc, namedesc, namestrain)
        call doNicheInvasion (numanc, numanctot, numdesc, time, lengthseq, activename, activepop, realname, &
          ntotalpop, strainparent, popfornascent, numstrain, div, nameanc, namedesc, namestrain, ancdist)
      endif
      if (event .eq. EVENT_POPULATION_DRIFT) then
        call choosepopdrift (popfordrift, activepop, realname, xn, numstrain)
        call doDrift (popfordrift, ntotalpop, numanc, numanctot, numdesc, time, lengthseq, &
          numstrain, div, nameanc, namedesc, namestrain, ancdist)
      endif
      if (ntotalpop .le. 1) exit
    end do
    call divergenonultra (nu, numanc, nameanc, ancdist, identitymatrix)
    ! next, do binning
    call binningdanny (numcrit, crit, nu, nclustersarray, identitymatrix)
    ! Check whether a replicate is within a factor of (5.00, 2.00, 1.50, 1.25,
    ! 1.10, 1.05) for all bins
    call testForSuccessFit (nclustersarray, numcrit, realdata, success)
    return
  end subroutine fredMethod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Start the population simulation.
  !>
  !> @param[in]     npop           The number of ecotypes assumed to be in
  !>                                the environmental DNA sample.
  !> @param[in]     nu             The number of homologous gene sequences
  !>                                in the environmental sample.
  !> @param[out]    activepop      The number of active populations.
  !> @param[out]    ntotalpop      
  !> @param[out]    numstrain      An array of the number of strains in each
  !>                                ecotype.
  !> @param[out]    realname       The original number of a population.
  !> @param[out]    activename     The active name from the list of
  !>                                still-active populations.
  !> @param[out]    numanc         An array of the number of ancestors for
  !>                                each strain.
  !> @param[out]    numdesc        An array of the the number of descendants
  !>                                of each ancestor.
  !> @param[out]    numanctot      The total number of ancestors, over all
  !>                                strains.
  !> @param[out]    time           The time of the last ancestor described.
  !> @param[out]    ntotalpopwithinvented
  !> @param[out]    nameanc        A 2D array holding the names each strains
  !>                                ancestors.  nameanc(strain, ancestor)
  !> @param[out]    namedesc
  !> @param[out]    namestrain
  !> @param[out]    ancdist        A 2d array holding the actual divergence
  !>                                between each contemporary strain and its
  !>                                ancestors.  Note, the first ancestor of
  !>                                each strain is itself.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine startpops (npop, nu, activepop, ntotalpop, numstrain, realname, activename, numanc, numdesc, &
    numanctot, time, ntotalpopwithinvented, nameanc, namedesc, namestrain, ancdist)
    integer, intent(in)               :: npop
    integer, intent(in)               :: nu
    integer, intent(out)              :: activepop
    integer, intent(out)              :: ntotalpop
    integer, intent(inout)            :: numstrain(:)
    integer, intent(inout)            :: realname(:)
    integer, intent(inout)            :: activename(:)
    integer, intent(inout)            :: numanc(:)
    integer, intent(inout)            :: numdesc(:)
    integer, intent(out)              :: numanctot
    integer, intent(out)              :: ntotalpopwithinvented
    integer, intent(inout)            :: nameanc(:,:)
    integer, intent(inout)            :: namedesc(:,:)
    integer, intent(inout)            :: namestrain(:,:)
    double precision, intent(out)     :: time
    real, intent(inout)               :: ancdist(:,:)
    ! Local variables.
    integer :: strainname
    integer :: jpop
    integer :: jstrain
    ! Time is the time of the last ancestor described; it starts at 0
    activepop = npop
    call canonical (npop, nu, numstrain)
    ! activepop is the number of active populations (note some will disappear
    ! in the backwards simulation as they become invented in reverse)
    time = 0.0
    ntotalpop = nu
    ntotalpopwithinvented = nu
    ! here we assign strains to the various ecotypes
    strainname = 0
    do jpop = 1, npop
      realname(jpop) = jpop
      activename(jpop) = jpop
      ! the realname of a population is its original number;
      ! the activename is the name from the list of still-active populations
      do jstrain = 1, numstrain(jpop)
        strainname = strainname + 1
        namestrain(jpop, jstrain) = strainname
        ! here we set up each strain as its own ancestor;
        ! that is, the first ancestor of each strain is itself
        numanc(strainname) = 1
        nameanc(strainname, 1) = strainname
        ! numanc(i) is the number of ancestors for the ith strain
        ! nameanc(i,j) is the name of the jth ancestor of the ith strain
        ! each ancestor has all its descendants listed, as below
        ancdist(strainname, 1) = 0.0
        ! ancdist(i,k) is the actual divergence between contemporary strain i
        ! and its kth ancestor.  Note, the first ancestor of each strain is
        ! itself.
        numdesc(strainname) = 1
        ! numdesc (i) is the number of descendants of the ith ancestor
        ! note that the 1 through nu terminal nodes have their first
        ! ancestor as 1 through nu, respectively
        namedesc(strainname, 1) = strainname
      end do
    end do
    numanctot = nu
    ! numanctot is the total number of ancestors, over all strains
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
    integer, intent(in)               :: bin(:)
    integer, intent(in)               :: numcrit
    integer, intent(in)               :: realdata(:)
    logical, intent(out)              :: success(6)
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
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 5.0 .or. xbin / xreal .gt. 5.0) return
    end do
    ! Meets the 500% criterion.
    success(1) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 2.0 .or. xbin / xreal .gt. 2.0) return
    end do
    ! Meets the 200% criterion.
    success(2) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.5 .or. xbin / xreal .gt. 1.5) return
    end do
    ! Meets the 150% criterion.
    success(3) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.25 .or. xbin / xreal .gt. 1.25) return
    end do
    ! Meets the 125% criterion.
    success(4) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
      if (xreal / xbin .gt. 1.1 .or. xbin / xreal .gt. 1.1) return
    end do
    ! Meets the 110% criterion.
    success(5) = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = bin(jcrit)
      if (xreal .lt. 1.0d-6 .or. xbin .lt. 1.0d-6) return
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
  !> @param[in]     eligibleparent
  !> @param[in]     eligibleps    
  !> @param[in]     omega         The rate of niche invasion per eligible
  !>                                parental population, measured as niche
  !>                                invasions per nucleotide substitution in a
  !>                                given gene.
  !> @param[in]     sigma         The rate of periodic selection per eligible
  !>                                population, measured as periodic selection
  !>                                events per population per nucleotide
  !>                                substitution in a given gene.
  !> @param[in]     xn            The rate of population drift.
  !> @param[in,out] time          The total time (in substitutions) before the
  !>                                present for a coalescence event.
  !> @param[out]    event         The type of key event.
  !> @param[in]     activepop     The number of active populations.
  !> @param[in]     realname      
  !> @param[in]     numstrain     An array of the number of strains in each
  !>                                ecotype.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine whichEventAndWhen (eligibleparent, eligibleps, omega, sigma, xn, time, event, activepop, realname, numstrain)
    integer, intent(out)              :: event
    integer, intent(in)               :: eligibleparent
    integer, intent(in)               :: eligibleps
    integer, intent(in)               :: activepop
    integer, intent(in)               :: numstrain(:)
    integer, intent(in)               :: realname(:)
    double precision, intent(in)      :: omega
    double precision, intent(in)      :: sigma
    double precision, intent(in)      :: xn
    double precision, intent(inout)   :: time
    ! Local variables.
    integer          :: jpop
    double precision :: effectiveOmega
    double precision :: effectiveSigma
    double precision :: effectiveDrift
    double precision :: rateKey
    double precision :: timeWait
    real             :: x
    ! note that omega is a per population rate of niche invasion,
    ! so we multiply by the number of eligible parent populations
    effectiveOmega = eligibleparent * omega
    effectiveSigma = eligibleps * sigma
    effectiveDrift = 0.0
    do jpop = 1, activepop
      effectiveDrift = effectiveDrift + (1.0 / xn) * numstrain(realname(jpop)) * (numstrain(realname(jpop)) - 1.0) * 0.5
    end do
    ! rateKey is the total rate of all key events.
    rateKey = effectiveOmega + effectiveSigma + effectiveDrift
    ! The following uses 0.01 / rateKey as the length of time over which the probability of getting a key event is 0.01.
    ! The full formula for timeWait is:  timeWait = (-1 * log (x) / 0.01) * (0.01 / rateKey)
    call randomNumber (x)
    timeWait = -1.0 * log (x) / rateKey
    ! time is incremented by the timeWait
    ! time is the total time (in substitutions) before the present for a coalescence event
    time = time + timeWait
    ! now, what kind of key event
    call randomNumber (x)
    if (x .lt. effectiveOmega / rateKey) then
      event = EVENT_NICHE_INVASION
    else if (x .lt. (effectiveOmega + effectiveSigma) / rateKey) then
      event = EVENT_PERIODIC_SELECTION
    else if (x .lt. (effectiveOmega + effectiveSigma + effectiveDrift) / rateKey) then
      event = EVENT_POPULATION_DRIFT
    endif
    return
  end subroutine whichEventAndWhen

end module methods
