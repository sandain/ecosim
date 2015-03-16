!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A program for function minimization using the simplex method.
! The minimum found will often be a local, not a global, minimum.
!
! For details, see Nelder & Mead, The Computer Journal, January 1965 [1].
!
! Programmed by DE Shaw,
! CSIRO, Division of Mathematics & Statistics
! P.O. Box 218, Lindfield, N.S.W. 2070
!
! With amendments by RWM Wedderburn
! Rothamsted Experimental Station
! Harpenden, Hertfordshire, England
!
! Further amended by Alan Miller
! CSIRO Division of Mathematics & Statistics
! Private Bag 10, Clayton, VIC. 3168
!
! Further amended by Fred Cohan and Jason M Wood
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> The simplexmethod module contains the Nelder-Mead simplex method for
!> function minimization.
!>
!> #### Reference
!>  [1] JA Nelder, and R Mead. 1965. A simplex Method for Function
!>        Minimization. The Computer Journal. 7(4):308-313.
!>        http://dx.doi.org/10.1093/comjnl/7.4.308
!>
!> @author DE Shaw
!> @author RWM Wedderburn
!> @author Alan Miller
!> @author Fred Cohan
!> @author Jason M Wood
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module simplexmethod
  use ISO_FORTRAN_ENV

  implicit none

  private

  ! Declare public methods.
  public :: nelmead
  public :: nelmeadFunction

  ! Declare private global parameters.
  real(kind = real64), parameter :: eta = epsilon(1.0d0) !< Epsilon.
  real(kind = real64), parameter :: neta = -1.0d0 * eta  !< Negative epsilon.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> An interface for functions to be optimized using the simplex method.
  !>
  !> @param[in]     nop           The number of parameters.
  !> @param[inout]  p             The parameters to be optimized.
  !> @param[out]    func          The function value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine nelmeadFunction (nop, p, func)
      use ISO_FORTRAN_ENV
      integer(kind = int32), intent(in)  :: nop
      real(kind = real64), intent(inout) :: p(nop)
      real(kind = real64), intent(out)   :: func
    end subroutine nelmeadFunction
  end interface

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A program for function minimization using the simplex method.
  !>
  !> @param[in,out] p             The values of the parameters.
  !> @param[in,out] step          The initial step sizes.
  !> @param[in]     nop           The number of parameters, including any to
  !>                                be held fixed.
  !> @param[out]    func          The function value corresponding to the
  !>                                final parameter values.
  !> @param[in]     max           The maximum number of function evaluations
  !>                                allowed. Say, 20 times the number of
  !>                                parameters, nop.
  !> @param[in]     iprint        The print control parameter.
  !>                                < 0 No printing.
  !>                                = 0 Printing of parameter values and the
  !>                                    function value after initial evidence
  !>                                    of convergence.
  !>                                > 0 As for iprint = 0 plus progress
  !>                                    report after every iprint
  !>                                    evaluations, plus printing for the
  !>                                    initial simplex.
  !> @param[in]     stopcr        Stopping criterion. The criterion is
  !>                                applied to the standard deviation of the
  !>                                values of func at the points of the
  !>                                simplex.
  !> @param[in]     nloop         The stopping rule is applied after every
  !>                                nloop function evaluations. Normally
  !>                                nloop should be slightly greater than
  !>                                nop, say nloop = 2 * nop.
  !> @param[in]     iquad         = 1 if fitting of a quadratic surface is
  !>                                required.
  !>                              = 0 if not.
  !>                              The fitting of a quadratic surface is
  !>                                strongly recommended, provided that the
  !>                                fitted function is continuous in the
  !>                                vicinity of the minimum. It is often a
  !>                                good indicator of whether a premature
  !>                                termination of the search has occurred.
  !> @param[in]     simp          Criterion for expanding the simplex to
  !>                                overcome rounding errors before fitting
  !>                                the quadratic surface. The simplex is
  !>                                expanded so that the function values at
  !>                                the points of the simplex exceed those
  !>                                at the supposed minimum by at least an
  !>                                amount simp.
  !> @param[out]    var           Contains the diagonal elements of the
  !>                                inverse of the information matrix.
  !> @param[in]     functn        Name of the user's subroutine - Arguments
  !>                                (nop, p, func) which returns the function
  !>                                value for a given set of parameter values
  !>                                in array p.
  !> @param[out]    ifault        = 0 for successful termination.
  !>                              = 1 if maximum number of function
  !>                                evaluations exceeded.
  !>                              = 2 if information matrix is not +ve
  !>                                semi-definite.
  !>                              = 3 if nop < 1.
  !>                              = 4 if nloop < 1.
  !> @param[in]     lout          The file handle to output to.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nelmead (p, step, nop, func, max, iprint, stopcr, nloop, &
    iquad, simp, var, functn, ifault, lout)
    real(kind = real64), intent(inout)   :: p(nop)
    real(kind = real64), intent(inout)   :: step(nop)
    real(kind = real64), intent(out)     :: func
    real(kind = real64), intent(in)      :: stopcr
    real(kind = real64), intent(in)      :: simp
    real(kind = real64), intent(out)     :: var(nop)
    integer(kind = int32), intent(in)    :: nop
    integer(kind = int32), intent(in)    :: max
    integer(kind = int32), intent(in)    :: iprint
    integer(kind = int32), intent(in)    :: nloop
    integer(kind = int32), intent(in)    :: iquad
    integer(kind = int32), intent(out)   :: ifault
    integer(kind = int32), intent(in)    :: lout
    procedure(nelmeadFunction), pointer  :: functn
    ! Local variables.
    real(kind = real64)            :: g(nop + 1, nop)
    real(kind = real64)            :: h(nop + 1)
    real(kind = real64)            :: pbar(nop)
    real(kind = real64)            :: pstar(nop)
    real(kind = real64)            :: pstst(nop)
    real(kind = real64)            :: aval(nop)
    real(kind = real64)            :: bmat(nop * (nop + 1) / 2)
    real(kind = real64)            :: pmin(nop)
    real(kind = real64)            :: vc(nop * (nop + 1) / 2)
    real(kind = real64)            :: temp(nop)
    real(kind = real64)            :: a0
    real(kind = real64)            :: hmin
    real(kind = real64)            :: hmax
    real(kind = real64)            :: hmean
    real(kind = real64)            :: hstar
    real(kind = real64)            :: hstd
    real(kind = real64)            :: hstst
    real(kind = real64)            :: rmax
    real(kind = real64)            :: savemn
    real(kind = real64)            :: test
    real(kind = real64)            :: ymin
    real(kind = real64), parameter :: a = 1.0d0 !< Reflection Coefficient.
    real(kind = real64), parameter :: b = 0.5d0 !< Contraction Coefficient.
    real(kind = real64), parameter :: c = 2.0d0 !< Expansion Coefficient.
    integer(kind = int32)          :: i
    integer(kind = int32)          :: i1
    integer(kind = int32)          :: i2
    integer(kind = int32)          :: iflag
    integer(kind = int32)          :: ii
    integer(kind = int32)          :: ij
    integer(kind = int32)          :: ijk
    integer(kind = int32)          :: imin
    integer(kind = int32)          :: imax
    integer(kind = int32)          :: irank
    integer(kind = int32)          :: irow
    integer(kind = int32)          :: j
    integer(kind = int32)          :: j1
    integer(kind = int32)          :: jj
    integer(kind = int32)          :: k
    integer(kind = int32)          :: l
    integer(kind = int32)          :: loop
    integer(kind = int32)          :: nap
    integer(kind = int32)          :: neval
    integer(kind = int32)          :: nmore
    integer(kind = int32)          :: np1
    integer(kind = int32)          :: nullty
    !
    ! if progress reports have been requested, print heading.
    !
    if (iprint .gt. 0) write (unit = lout, fmt = 1000) iprint
1000 format (' Progress Report every', I4, ' function evaluations'/, &
      ' EVAL.   FUNC.VALUE.', 10X, 'PARAMETER VALUES')
    !
    ! check input arguments.
    !
    ifault = 0
    if (nop .le. 0) ifault = 3
    if (nloop .le. 0) ifault = 4
    if (ifault .ne. 0) return
    !
    ! set nap = no. of parameters to be varied, i.e. with step .ne. 0.
    !
    nap = 0
    neval = 0
    loop = 0
    iflag = 0
    do i = 1, nop
      ! check to see if step(i) is not equal to zero.
      if (step(i) .lt. neta .or. step(i) .gt. eta) then
        nap = nap + 1
      end if
    end do
    !
    ! if nap = 0 evaluate function at the starting point and return.
    !
    if (nap .gt. 0) goto 30
    call functn (nop, p, func)
    return
    !
    ! set up the initial simplex.
    !
30  continue
    do i = 1, nop
      g(1, i) = p(i)
    end do
    irow = 2
    do 60 i = 1, nop
    ! check to see if step(i) is equal to zero.
    if (step(i) .gt. neta .and. step(i) .lt. eta) then
      goto 60
    end if
    do j = 1, nop
      g(irow, j) = p(j)
    end do
    g(irow, i) = p(i) + step(i)
    irow = irow + 1
60  CONTINUE
    !
    np1 = nap + 1
    do 90 i = 1, np1
    do j = 1, nop
      p(j) = g(i, j)
    end do
    call functn (nop, p, h(i))
    neval = neval + 1
    if (iprint .le. 0) goto 90
    write (unit = lout, fmt = 1010) neval, h(i), (p(j), j = 1, nop)
1010 format (/1X, I4, 2X, G12.5, 2X, 5G11.4, 3(/21X, 5G11.4))
90  CONTINUE
    !
    ! start of main cycle.
    !
    ! find max. & min. values for current simplex (hmax & hmin).
    !
100 loop = loop + 1
    imax = 1
    imin = 1
    hmax = h(1)
    hmin = h(1)
    do 120 i = 2, np1
    if (h(i) .le. hmax) goto 110
    imax = i
    hmax = h(i)
    goto 120
110 if (h(i) .ge. hmin) goto 120
    imin = i
    hmin = h(i)
120 CONTINUE
    !
    ! find the centroid of the vertices other than p(imax).
    !
    do i = 1, nop
      pbar(i) = 0.0d0
    end do
    do 150 i = 1, np1
    if (i .eq. imax) goto 150
    do J = 1, nop
      pbar(j) = pbar(j) + g(i, j)
    end do
150 CONTINUE
    do j = 1, nop
      pbar(j) = pbar(j) / float (nap)
    end do
    !
    ! reflect maximum through pbar to pstar,
    ! hstar = function value at pstar.
    !
    do i = 1, nop
      pstar(i) = a * (pbar(i) - g(imax, i)) + pbar(i)
    end do
    call functn (nop, pstar, hstar)
    neval = neval + 1
    if (iprint .le. 0) goto 180
    if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) &
      neval, hstar, (pstar(j), j = 1, nop)
    !
    ! if hstar < hmin, reflect pbar through pstar,
    ! hstst = function value at pstst.
    !
180 if (hstar .ge. hmin) goto 220
    do i = 1, nop
      pstst(i) = c * (pstar(i) - pbar(i)) + pbar(i)
    end do
    call functn (nop, pstst, hstst)
    neval = neval + 1
    if (iprint .le. 0) goto 200
    if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) &
      neval, hstst, (pstst(j), j = 1, nop)
    !
    ! if hstst < hmin replace current maximum point by pstst and
    ! hmax by hstst, then test for convergence.
    !
200 if (hstst .ge. hmin) goto 320
    do i = 1, nop
      ! check to see if step(i) is not equal to zero.
      if (step(i) .lt. neta .or. step(i) .gt. eta) then
        g(imax, i) = pstst(i)
      end if
    end do
    h(imax) = hstst
    goto 340
    !
    ! hstar is not < hmin.
    ! test whether it is < function value at some point other than
    ! p(imax). if it is replace p(imax) by pstar & hmax by hstar.
    !
220 do 230 i = 1, np1
    if (i .eq. imax) goto 230
    if (hstar .lt. h(i)) goto 320
230 CONTINUE
    !
    ! hstar > all function values except possibly hmax.
    ! if hstar < = hmax, replace p(imax) by pstar & hmax by hstar.
    !
    if (hstar .gt. hmax) goto 260
    do i = 1, nop
      ! check to see if step(i) is not equal to zero.
      if (step(i) .lt. neta .or. step(i) .gt. eta) then
        g(imax, i) = pstar(i)
      end if
    end do
    hmax = hstar
    h(imax) = hstar
    !
    ! contracted step to the point pstst,
    ! hstst = function value at pstst.
    !
260 do i = 1, nop
      pstst(i) = b * g(imax, i) + (1.0 - b) * pbar(i)
    end do
    call functn (nop, pstst, hstst)
    neval = neval + 1
    if (iprint .le. 0) goto 280
    if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) &
      neval, hstst, (pstst(j), j = 1, nop)
    !
    ! if hstst < hmax replace p(imax) by pstst & hmax by hstst.
    !
280 if (hstst .gt. hmax) goto 300
    do i = 1, nop
      ! check to see if step(i) is not equal to zero.
      if (step(i) .lt. neta .or. step(i) .gt. eta) then
        g(imax, i) = pstst(i)
      end if
    end do
    h(imax) = hstst
    goto 340
    !
    ! hstst > hmax.
    ! shrink the simplex by replacing each point, other than the current
    ! minimum, by a point mid-way between its current position and the
    ! minimum.
    !
300 do 315 i = 1, np1
    if (i .eq. imin) goto 315
    do j = 1, nop
      ! check to see if step(j) is not equal to zero.
      if (step(j) .lt. neta .or. step(j) .gt. eta) then
        g(i, j) = (g(i, j) + g(imin, j)) * 0.5
      end if
      p(j) = g(i, j)
    end do
    call functn (nop, p, h(i))
    neval = neval + 1
    if (iprint .le. 0) goto 315
    if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) &
      neval, h(i), (p(j), j = 1, nop)
315 CONTINUE
    goto 340
    !
    ! replace maximum point by pstar & h(imax) by hstar.
    !
320 do i = 1, nop
      ! check to see if step(i) is not equal to zero.
      if (step(i) .lt. neta .or. step(i) .gt. eta) then
        g(imax, i) = pstst(i)
      end if
    end do
    h(imax) = hstar
    !
    ! if loop = nloop test for convergence, otherwise repeat main cycle.
    !
340 if (loop .lt. nloop) goto 100
    !
    ! calculate mean & standard deviation of function values for the
    ! current simplex.
    !
    hstd = 0.0d0
    hmean = 0.0d0
    do i = 1, np1
      hmean = hmean + h(i)
    end do
    hmean = hmean / float (np1)
    do i = 1, np1
      hstd = hstd + (h(i) - hmean) ** 2
    end do
    hstd = sqrt (hstd / float (np1))
    !
    ! if the rms > stopcr, set iflag & loop to zero and goto the
    ! start of the main cycle again.
    !
    if (hstd .le. stopcr .or. neval .gt. max) goto 410
    iflag = 0
    loop = 0
    goto 100
    !
    ! find the centroid of the current simplex and the function value there.
    !
410 do 380 i = 1, nop
    ! check to see if step(i) is equal to zero.
    if (step(i) .gt. neta .and. step(i) .lt. eta) then
      goto 380
    end if
    p(i) = 0.0d0
    do j = 1, np1
      p(i) = p(i) + g(j, i)
    end do
    p(i) = p(i) / float (np1)
380 CONTINUE
    call functn (nop, p, func)
    neval = neval + 1
    if (iprint .le. 0) goto 390
    if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) &
      neval, func, (p(j), j = 1, nop)
    !
    ! test whether the no. of function values allowed, max, has been
    ! overrun; if so, exit with ifault = 1.
    !
390 if (neval .le. max) goto 420
    ifault = 1
    if (iprint .lt. 0) return
    write (unit = lout, fmt = 1020) max
1020 format (' No. of function evaluations > ', I5)
    write (unit = lout, fmt = 1030) hstd
1030 format (' RMS of function values of last simplex =', G14.6)
    write (unit = lout, fmt = 1040) (p(i), i = 1, nop)
1040 format (' Centroid of last simplex =', 4(/1X, 6G13.5))
    write (unit = lout, fmt = 1050) func
1050 format (' Function value at centroid =', G14.6)
    return
    !
    ! convergence criterion satisfied.
    ! if iflag = 0, set iflag & save hmean.
    ! if iflag = 1 & change in hmean < = stopcr then search is complete.
    !
420 if (iprint .lt. 0) goto 430
    write (unit = lout, fmt = 1060)
1060 format (/' EVIDENCE OF CONVERGENCE')
    write (unit = lout, fmt = 1040) (p(i), i = 1, nop)
    write (unit = lout, fmt = 1050) func
430 if (iflag .gt. 0) goto 450
    iflag = 1
440 savemn = hmean
    loop = 0
    goto 100
450 if (abs (savemn - hmean) .ge. stopcr) goto 440
    if (iprint .lt. 0) goto 460
    write (unit = lout, fmt = 1070) neval
1070 format (//' Minimum found after', I5, ' function evaluations')
    write (unit = lout, fmt = 1080) (p(i), i = 1, nop)
1080 format (' Minimum at', 4(/1X, 6G13.6))
    write (unit = lout, fmt = 1090) func
1090 format (' Function value at minimum =', G14.6)
460 if (iquad .le. 0) return
    !
    !------------------------------------------------------------------
    !
    ! quadratic surface fitting.
    !
    if (iprint .ge. 0) write (unit = lout, fmt = 1110)
1110 format (/' Fitting quadratic surface about supposed minimum'/)
    !
    ! expand the final simplex, if necessary, to overcome rounding
    ! errors.
    !
    hmin = func
    nmore = 0
    do 490 i = 1, np1
470 test = abs (h(i) - func)
    if (test .ge. simp) goto 490
    do j = 1, nop
      ! check to see if step(j) is not equal to zero.
      if (step(j) .lt. neta .or. step(j) .gt. eta) then
        g(i, j) = (g(i, j) - p(j)) + g(i, j)
      end if
      pstst(j) = g(i, j)
    end do
    call functn (nop, pstst, h(i))
    nmore = nmore + 1
    neval = neval + 1
    if (h(i) .ge. hmin) goto 470
    hmin = h(i)
    if (iprint .ge. 0) write (unit = lout, fmt = 1010) &
      neval, hmin, (pstst(j), j = 1, nop)
    goto 470
490 CONTINUE
    !
    ! function values are calculated at an additional nap points.
    !
    do i = 1, nap
      i1 = i + 1
      do j = 1, nop
        pstar(j) = (g(1, j) + g(i1, j)) * 0.5
      end do
      call functn (nop, pstar, aval(i))
      nmore = nmore + 1
      neval = neval + 1
    end do
    !
    ! the matrix of estimated second derivatives is calculated and its
    ! lower triangle is stored in bmat.
    !
    a0 = h(1)
    do 540 i = 1, nap
    i1 = i - 1
    i2 = i + 1
    if (i1 .lt. 1) goto 540
    do j = 1, i1
      j1 = j + 1
      do k = 1, nop
        pstst(k) = (g(i2, k) + g(j1, k)) * 0.5
      end do
      call functn (nop, pstst, hstst)
      nmore = nmore + 1
      neval = neval + 1
      l = i * (i - 1) / 2 + j
      bmat(l) = 2.0 * (hstst + a0 - aval(i) - aval(j))
    end do
540 CONTINUE
    l = 0
    do i = 1, nap
      i1 = i + 1
      l = l + i
      bmat(l) = 2.0 * (h(i1) + a0 - 2.0 * aval(i))
    end do
    !
    ! the vector of estimated first derivatives is calculated and
    ! stored in aval.
    !
    do i = 1, nap
      i1 = i + 1
      aval(i) = 2.0 * aval(i) - (h(i1) + 3.0 * a0) * 0.5
    end do
    !
    ! the matrix q of Nelder & Mead is calculated and stored in g.
    !
    do i = 1, nop
      pmin(i) = g(1, i)
    end do
    do 580 i = 1, nap
    i1 = i + 1
    do 580 j = 1, nop
    g(i1, j) = g(i1, j) - g(1, j)
580 CONTINUE
    do 590 i = 1, nap
    i1 = i + 1
    do 590 j = 1, nop
    g(i, j) = g(i1, j)
590 CONTINUE
    !
    ! invert bmat.
    !
    call syminv (bmat, nap, bmat, temp, nullty, ifault, rmax)
    if (ifault .ne. 0) goto 600
    irank = nap - nullty
    goto 610
600 if (iprint .ge. 0) write (unit = lout, fmt = 1120)
1120 format (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/ &
      ' MINIMUM PROBABLY NOT FOUND'/)
    ifault = 2
    if (neval .gt. max) return
    write (unit = lout, fmt = 1130)
1130 format (/10X, 'Search restarting'/)
    do 605 i = 1, nop
605 step(i) = 0.5 * step(i)
    goto 30
    !
    ! bmat * a / 2 is calculated and stored in h.
    !
610 do 650 i = 1, nap
    h(i) = 0.0d0
    do 640 j = 1, nap
    if (j .gt. i) goto 620
    l = i * (i - 1) / 2 + j
    goto 630
620 l = j * (j - 1) / 2 + i
630 h(i) = h(i) + bmat(l) * aval(j)
640 CONTINUE
650 CONTINUE
    !
    ! find the position, pmin, & value, ymin, of the minimum of the
    ! quadratic.
    !
    ymin = 0.0d0
    do i = 1, nap
      ymin = ymin + h(i) * aval(i)
    end do
    ymin = a0 - ymin
    do 670 i = 1, nop
    pstst(i) = 0.0d0
    do 670 j = 1, nap
670 pstst(i) = pstst(i) + h(j) * g(j, i)
    do 680 i = 1, nop
680 pmin(i) = pmin(i) - pstst(i)
    if (iprint .lt. 0) goto 690
    write (unit = lout, fmt = 1140) ymin, (pmin(i), i = 1, nop)
1140 format (' Minimum of quadratic surface =', G14.6, ' at', &
      4(/1X, 6G13.5))
    write (unit = lout, fmt = 1150)
1150 format (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED', &
      1X, 'FROM THE MINIMIZATION, '/ &
      ' THE MINIMUM MAY BE FALSE &/OR THE INFORMATION MATRIX MAY BE', &
      1X, 'INACCURATE'/)
    !
    ! q * bmat * q' / 2 is calculated & its lower triangle stored in vc.
    !
690 do 760 i = 1, nop
    do 730 j = 1, nap
    h(j) = 0.0d0
    do 720 k = 1, nap
    if (k .gt. j) goto 700
    l = j * (j - 1) / 2 + k
    goto 710
700 l = k * (k - 1) / 2 + j
710 h(j) = h(j) + bmat(l) * g(k, i) * 0.5
720 CONTINUE
730 CONTINUE
    do 750 j = i, nop
    l = j * (j - 1) / 2 + i
    vc(l) = 0.0d0
    do k = 1, nap
      vc(l) = vc(l) + h(k) * g(k, j)
    end do
750 CONTINUE
760 CONTINUE
    !
    ! the diagonal elements of vc are copied into var.
    !
    j = 0
    do i = 1, nop
      j = j + i
      var(i) = vc(j)
    end do
    if (iprint .lt. 0) return
    write (unit = lout, fmt = 1160) irank
1160 format (' Rank of information matrix =', I3/ &
      ' Inverse of information matrix: - ')
    ijk = 1
    goto 880
    !
790 continue
    write (unit = lout, fmt = 1170)
1170 format (/' if the function minimized was - Log(likelihood), '/ &
      ' this is the covariance matrix of the parameters.'/ &
      ' if the function was a sum of squares of residuals, '/ &
      ' this matrix must be multiplied by twice the estimated', &
      1X, 'residual variance'/' to obtain the covariance matrix.'/)
    call syminv (vc, nap, bmat, temp, nullty, ifault, rmax)
    !
    ! bmat now contains the information matrix.
    !
    write (unit = lout, fmt = 1190)
1190 format (' INFORMATION MATRIX: - '/)
    ijk = 3
    goto 880
    !
800 ijk = 2
    ii = 0
    ij = 0
    do i = 1, nop
      ii = ii + i
      if (vc(ii) .gt. 0.0d0) then
        vc(ii) = 1.0 / sqrt (vc(ii))
      else
        vc(ii) = 0.0d0
      end if
      jj = 0
      do j = 1, i - 1
        jj = jj + j
        ij = ij + 1
        vc(ij) = vc(ij) * vc(ii) * vc(jj)
      end do
      ij = ij + 1
    end do
    !
    write (unit = lout, fmt = 1200)
1200 format (//' CORRELATION MATRIX: - ')
    ii = 0
    do i = 1, nop
      ii = ii + i
      ! check to see if vc(ii) is not equal to zero.
      if (vc(ii) .lt. neta .or. vc(ii) .gt. eta) then
        vc(ii) = 1.0d0
      end if
    end do
    goto 880
    !
    ! exit, on successful termination.
    !
860 write (unit = lout, fmt = 1210) nmore
1210 format (/' A further', I4, ' function evaluations have been used'/)
    return
    !
880 l = 1
890 if (l .gt. nop) goto (790, 860, 800), ijk
    ii = l * (l - 1) / 2
    do 910 i = l, nop
    i1 = ii + l
    ii = ii + i
    i2 = min0 (ii, i1 + 5)
    if (ijk .eq. 3) goto 900
    write (unit = lout, fmt = 1230) (vc(j), j = i1, i2)
    goto 910
900 write (unit = lout, fmt = 1230) (bmat(j), j = i1, i2)
910 CONTINUE
1230 format (1X, 6G13.5)
    write (unit = lout, fmt = 1240)
1240 format (/)
    l = l + 6
    goto 890
  end subroutine nelmead

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Algorithm AS7, Applied Statistics, Volume 17, 1968.
  !>
  !> @param[in,out] a             The symmetric matrix to be inverted, stored
  !>                                in lower-triangular form.
  !> @param[in]     n             Order of the matrix.
  !> @param[in,out] c             The inverse of a (a generalized inverse if
  !>                                c is singular), also stored in lower
  !>                                triangular. c and a may occupy the same
  !>                                locations.
  !> @param[in,out] w             Workspace.
  !> @param[out]    nullty        The rank deficiency of a.
  !> @param[out]    ifault        Error indicator.
  !>                                = 1 if n < 1.
  !>                                = 2 if a is not +ve semi-definite.
  !>                                = 0 otherwise.
  !> @param[out]    rmax          Approximate bound on the accuracy of the
  !>                                diagonal elements of c. if
  !>                                rmax = 1.e-04 then the diagonal
  !>                                elements of c will be accurate to about
  !>                                4 decimal digits.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine syminv (a, n, c, w, nullty, ifault, rmax)
    real(kind = real64), intent(inout)   :: a(:)
    real(kind = real64), intent(inout)   :: c(:)
    real(kind = real64), intent(inout)   :: w(n)
    real(kind = real64), intent(out)     :: rmax
    integer(kind = int32), intent(in)    :: n
    integer(kind = int32), intent(out)   :: nullty
    integer(kind = int32), intent(out)   :: ifault
    ! Local Variables
    real(kind = real64)   :: x
    integer(kind = int32) :: i
    integer(kind = int32) :: icol
    integer(kind = int32) :: irow
    integer(kind = int32) :: j
    integer(kind = int32) :: jcol
    integer(kind = int32) :: k
    integer(kind = int32) :: l
    integer(kind = int32) :: ndiag
    integer(kind = int32) :: nn
    integer(kind = int32) :: nrow
    integer(kind = int32) :: mdiag
    nrow = n
    ifault = 1
    if (nrow .le. 0) goto 100
    ifault = 0
    !
    ! Cholesky factorization of a, result in c
    !
    call chola (a, nrow, c, nullty, ifault, rmax, w)
    if (ifault .ne. 0) goto 100
    !
    ! invert c & form the product (cinv)' * cinv, where cinv is the inverse
    ! of c, row by row starting with the last row.
    ! irow = the row number, ndiag = location of last element in the row.
    !
    nn = nrow * (nrow + 1) / 2
    irow = nrow
    ndiag = nn
    ! check to see if c(ndiag) is equal to zero.
10  if (c(ndiag) .gt. neta .and. c(ndiag) .lt. eta) goto 60
    l = ndiag
    do i = irow, nrow
      w(i) = c(l)
      l = l + i
    end do
    icol = nrow
    jcol = nn
    mdiag = nn
30  l = jcol
    x = 0.0d0
    if (icol .eq. irow) x = 1.0d0 / w(irow)
    k = nrow
40  if (k .eq. irow) goto 50
    x = x - w(k) * c(l)
    k = k - 1
    l = l - 1
    if (l .gt. mdiag) l = l - k + 1
    goto 40
50  c(l) = x / w(irow)
    if (icol .eq. irow) goto 80
    mdiag = mdiag - icol
    icol = icol - 1
    jcol = jcol - 1
    goto 30
60  l = ndiag
    do j = irow, nrow
      c(l) = 0.0d0
      l = l + j
    end do
80  ndiag = ndiag - irow
    irow = irow - 1
    if (irow .ne. 0) goto 10
100 return
  end subroutine syminv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Algorithm AS6, Applied Statistics, Volume 17, 1968
  !>
  !> Modifications by AJ Miller
  !>
  !> @param[in]  a                A +ve definite matrix stored in
  !>                                lower-triangular form.
  !> @param[in]  n                The order of a.
  !> @param[out] u                A lower-triangular matrix such that
  !>                                u * u' = a.
  !>                                a & u may occupy the same locations.
  !> @param[out] nullty           The rank deficiency of a.
  !> @param[out] ifault           Error indicator.
  !>                                = 1 if n < 1.
  !>                                = 2 if a is not +ve semi-definite.
  !>                                = 0 otherwise.
  !> @param[out] rmax             An estimate of the relative accuracy of the
  !>                                diagonal elements of u.
  !> @param[out] r                Array containing bounds on the relative
  !>                                 accuracy of each diagonal element of u.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine chola (a, n, u, nullty, ifault, rmax, r)
    real(kind = real64), intent(in)      :: a(:)
    real(kind = real64), intent(out)     :: u(:)
    real(kind = real64), intent(out)     :: rmax
    real(kind = real64), intent(out)     :: r(:)
    integer(kind = int32), intent(in)    :: n
    integer(kind = int32), intent(out)   :: nullty
    integer(kind = int32), intent(out)   :: ifault
    ! Local Variables
    real(kind = real64)            :: rsq
    real(kind = real64)            :: w
    integer(kind = int32)          :: i
    integer(kind = int32)          :: icol
    integer(kind = int32)          :: irow
    integer(kind = int32)          :: j
    integer(kind = int32)          :: k
    integer(kind = int32)          :: l
    integer(kind = int32)          :: m
    ifault = 1
    if (n .le. 0) return
    ifault = 2
    nullty = 0
    rmax = eta
    r(1) = eta
    j = 1
    k = 0
    rsq = 0.0d0
    !
    ! factorize column by column, icol = column no.
    !
    do icol = 1, n
      l = 0
      !
      ! irow = row number within column icol.
      !
      do irow = 1, icol
        k = k + 1
        w = a(k)
        if (irow .eq. icol) rsq = (w * eta) ** 2
        m = j
        do i = 1, irow
          l = l + 1
          if (i .eq. irow) exit
          w = w - u(l) * u(m)
          if (irow .eq. icol) rsq = rsq + (u(l) ** 2 * r(i)) ** 2
          m = m + 1
        end do
        if (irow .ne. icol) exit
        if (u(l) .gt. 0.0d0) then
          u(k) = w / u(l)
        else
          u(k) = 0.0d0
          if (abs (w) .gt. abs (rmax * a(k))) return
        end if
      end do
      !
      ! end of row, estimate relative accuracy of diagonal element.
      !
      rsq = sqrt (rsq)
      if (abs (w) .gt. 5.0d0 * rsq) then
        if (w .lt. 0.0d0) return
        u(k) = sqrt (w)
        r(i) = rsq / w
        if (r(i) .gt. rmax) rmax = r(i)
      else
        u(k) = 0.0d0
        nullty = nullty + 1
      end if
      j = j + icol
    end do
    ifault = 0
    return
  end subroutine chola

end module simplexmethod
