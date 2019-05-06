!*******************************************************************************
!  Read and interpolate Hades Acceptance Matrix for Theorists
!
!  Code version 2.0 of February 15, 2011
!
!  Created :  07/04/05 R. Holzmann
!  Modified : 25/04/05 R. Holzmann  added new interpolation modes
!  Modified : 04/05/05 R. Holzmann  added resolution function
!  Modified : 24/06/05 R. Holzmann  added pair acceptance
!  Modified : 20/07/05 R. Holzmann  fixed error messages on input
!  Modified : 06/02/06 R. Holzmann  support for diff. matrix sizes (v1.0)
!  Modified : 15/04/08 R. Holzmann  added function to set resolution (v1.1)
!  Modified : 07/02/09 R. Holzmann  added p-dependant Eloss in smear4momentum()
!  Modified : 16/04/09 R. Holzmann  theta and phi resolution from Ar+KCl embedding
!  Modified:  15/02/11 R. Holzmann  added support for non-gaussian momentum smearing
!  Mofified:  20/04/11 J. Weil      use stream I/O, generic intrinsics, etc
!*******************************************************************************
! This file contains three modules:
! 1) HAFT_aux    (auxiliary routines)
! 2) HAFT_single (single-particle acceptance filtering)
! 3) HAFT_pair   (pair acceptance filtering)
!*******************************************************************************


!*******************************************************************************
!*******************************************************************************
! This module contains auxiliary functions and general constants,
! which are used in the modules HAFT_single and HAFT_pair.
!*******************************************************************************
!*******************************************************************************
module HAFT_aux

  implicit none
  private

  public :: kernel, sampleGauss, sampleMP, interpol

  real, parameter, public :: pi = 3.141592654

  real, parameter :: lg10 = log(10.)

contains

  !*****************************************************************************
  !  Compute interpolation kernel
  !
  !  mode = 0: nearest neighbour
  !       = 1: piece-wise linear
  !       = 2: piece-wise quadratic
  !       =-2: quadratic B-spline
  !       = 3: cubic spline (B=0,C=1)
  !       =-3: cubic B-spline (B=1, C=0)
  !       = 4: cubic Catmull-Rom spline (B=0, C=1/2)
  !       =-4: cubic "optimal" cardinal spline (B=1/3, C=1/3)
  !
  !  mode >=0: interpolating (i.e. exact at grid points)
  !  mode < 0: approximating (i.e. not exact, but smoother)
  !*****************************************************************************
  real function kernel(u,mode)

      real, intent(in) :: u
      integer, intent(in) :: mode

      real :: ua

      select case (mode)
      case (0)  ! nearest neighbour
        if (u>-0.5 .and. u<=0.5) then
          kernel = 1.
        else
          kernel = 0.
        end if

      case (1)  ! linear
        ua = abs(u)
        if (ua<1.) then
          kernel = 1.-ua
        else
          kernel = 0.
        end if

      case (2)  ! quadratic
        ua = abs(u)
        if (ua<=0.5) then
          kernel = 1.0-2.*ua*ua
        else if (ua<=1.5) then
          kernel = 1.5+(ua-2.5)*ua
        else
          kernel = 0.
        end if

      case (-2)  ! quadratic B-spline
        ua = abs(u)
        if (ua<=0.5) then
          kernel = 0.75-ua*ua
        else if (ua<=1.5) then
          kernel = 1.125+0.5*(ua-3.)*ua
        else
          kernel = 0.
        end if

      case (3)  ! cubic
        ua = abs(u)
        if (ua<=1.) then
          kernel = 1.+(ua-2.)*ua*ua
        else if (ua<=2.) then
          kernel = 4.+((5.-ua)*ua-8.)*ua
        else
          kernel = 0.
        end if

      case (-3)  ! cubic B-spline
        ua = abs(u)
        if (ua<=1.) then
          kernel = 2./3.+(0.5*ua-1.)*ua*ua
        else if (ua<=2.) then
          kernel = 4./3.+((1.-ua/6.)*ua-2.)*ua
        else
          kernel = 0.
        end if

      case (4)  ! cubic Catmull-Rom
        ua = abs(u)
        if (ua<=1.) then
          kernel = 1.+(1.5*ua-2.5)*ua*ua
        else if (ua<=2.) then
          kernel = 2.+((2.5-0.5*ua)*ua-4.)*ua
        else
          kernel = 0.
        end if

      case (-4)  ! optimal cubic cardinal spline
                 ! (=compromise between blurring and ringing)
        ua = abs(u)
        if (ua<=1.) then
          kernel = 8./9.+(7./6.*ua-2.)*ua*ua
        else if (ua<=2.) then
          kernel = 16./9.+((2.-7./18.*ua)*ua-10./3.)*ua
        else
          kernel = 0.
        end if

      case default  ! undefined mode
        kernel = 0.

      end select

  end function kernel


  !*****************************************************************************
  !  Return random number according to a normal distribution.
  !*****************************************************************************
  real function sampleGauss(mean,sigma)
      real, intent(in) :: mean, sigma

      real :: theta, r(2)

      sampleGauss = mean
      if (sigma<=0.) return
      call random_number(r)
      theta = 2.*pi*r(1)
      sampleGauss = mean + sigma*cos(theta)*sqrt(-2.*log(r(2)))
  end function sampleGauss


  !*****************************************************************************
  !     HADES momentum spread (pRec-pSim)/pSim
  !*****************************************************************************
  real function momSpread(x,respar,ns)

      real, intent(in) :: x, respar(10), ns

      real :: pos, sig, left, right, farleft, argn, argp, argn2, e2, amp

      e2 = exp(-0.5*ns*ns)

      pos = respar(1)     ! Mean
      sig = respar(2)     ! Sigma
      left = respar(3)    ! Par3 (>0)
      right = respar(4)   ! Par4 (<0)
      farleft = respar(5) ! Par5 (>0)

      if (x>=(pos-ns*sig) .or. x<=(-lg10/left+pos-ns*sig)) then
         argn = 0.
      else
         argn = 1.
      end if

      if (x>=(pos-ns*sig)) then
        argp = 1.
      else
        argp = 0.
      end if

      if (x>(-lg10/left+pos-ns*sig)) then
        argn2 = 0.
      else
        argn2 = 1.
      end if

      amp = e2  ! Gauss amplitude at +/-2 sigma

      momSpread = exp( -0.5*((x-pos)/sig)*((x-pos)/sig) )               &  ! Gauss
             + amp*exp(  left*(x-(pos-ns*sig)) )*argn                   &  ! left tail (connects to Gauss at pos-ns*sig
             + amp*exp( right*(x-(pos-ns*sig)) )*argp                   &  ! right tail (Gauss sits on top of it)
             + 0.1*amp*exp( farleft*(x-(-lg10/left+pos-ns*sig)) )*argn2    ! far left tail
                                                                           ! (joins left tail where decayed to 1/10)
  end function momSpread


  !*****************************************************************************
  !  Return random number according to the normalized HADES momentum distribution.
  !*****************************************************************************
  real function sampleMP(respar,ns)

      real, intent(in) :: respar(10), ns

      real :: pos, sig, left, right, farleft, A(0:3), F(0:3)
      real :: dx, ftest, r1, r2, r3, e2
      integer :: cnt, cnt1, cnt2, cnt3

      e2 = exp(-0.5*ns*ns)

      pos = respar(1)      ! centroid
      sig = respar(2)      ! width
      left = respar(3)     ! left slope
      right = respar(4)    ! right slope
      farleft = respar(5)  ! far left slope

      ! compute function amplitudes
      A(0) = 0.1*(1. + e2)
      A(1) = 1. + e2
      A(2) = 1. + exp(right*ns*sig)
      A(3) = (e2 + exp(right*2.*ns*sig))/exp(right*2.*ns*sig)

      ! compute function areas
      F(0) = A(0)/farleft * (1. - exp(farleft*(lg10/left+ns*sig-pos-1.)))  ! [-1, 1/10left]
      F(1) = A(1)/left * 9./10.                                            ! [1/10, pos-ns*sig]
      F(2) = A(2) * 2.*ns*sig                                              ! [pos-ns*sig, pos+ns*sig]
      F(3) = A(3)/right * (exp(right*(1.-pos+ns*sig))-exp(right*2.*ns*sig))! [pos + ns*sig, 1]

      F(0:3) = F(0:3)/sum(F)   ! normalize areas

      ! sample dx by comparing with piece-wise function

      do cnt=1,1000    ! allow max 1000 trials
         cnt1 = 0
         cnt2 = 0
         cnt3 = 0
         call random_number(r1)
         ! select region and sample test function
         if  (r1 < F(0)) then             ! far left tail

            do
               cnt1 = cnt1 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/farleft - ns*sig + pos - lg10/left
               if (cnt1==1000) write(6,*) 'cnt1=1000 ', pos, sig, farleft
               if (dx>=-1. .or. cnt1>=1000) exit  ! limit range to >=-1
            end do
            ftest = A(0) * exp(farleft*(dx+lg10/left-pos+ns*sig))

         else if (r1 < sum(F(0:1))) then      ! left tail

            do
               cnt2 = cnt2 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/left - ns*sig + pos
               if (cnt2==1000) write(6,*) 'cnt2=1000', pos, sig, left
               if (dx>=-lg10/left-ns*sig+pos .or. cnt2>=1000) exit
            end do
            ftest = A(1) * exp(left*(dx-pos+ns*sig))

         else if (r1 < sum(F(0:2))) then   ! peak region

            call random_number(r2)
            r2 = r2 - 0.5
            dx = 2.*ns*sig*r2 + pos
            ftest = A(2)

         else                            ! right tail

            do
               cnt3 = cnt3 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/right + ns*sig + pos
               if (cnt3==1000) write(6,*) 'cnt3=1000', pos, sig, right
               if (dx<=1. .or. cnt3>=1000) exit  ! limit range to <=1
            end do
            ftest = A(3) * exp(right*(dx-pos+ns*sig))

         end if

         ! do rejection test
         sampleMP = dx

         call random_number(r3)
         if ( r3<momSpread(dx,respar,ns)/ftest ) return

      end do

      write(6,*) 'cnt=1000'
      sampleMP = 0.

  end function sampleMP


  !*****************************************************************************
  !     linear interpolation in table (xtab,ytab)
  !*****************************************************************************
  real function interpol(x,xtab,ytab,n)
      real, intent(in) :: x, xtab(*), ytab(*)
      integer, intent(in) :: n

      integer :: i
      real :: a, b

      if (x<=xtab(1)) then ! below table range
        interpol = ytab(1)
        return
      else if (x>=xtab(n)) then ! above table range
        interpol = ytab(n)
        return
      end if

      do i=2,n
        interpol = ytab(i)
        if (x==xtab(i)) return
        if (x<xtab(i)) exit
      end do

      a = ytab(i-1)
      b = (ytab(i)-ytab(i-1))/(xtab(i)-xtab(i-1))
      interpol = a + (x-xtab(i-1))*b
  end function interpol


end module HAFT_aux



!*******************************************************************************
!*******************************************************************************
! This module contains routines for single-particle acceptance filtering.
!
!  Usage : 1) set acceptance file name with
!                 call setFileName(fname)
!          2) sample single acceptance values with calls to
!                 acc = getHadesAcceptance(id,mom,theta,phi,mode)
!          3) apply detector resolution with calls to
!                 call smearhadesmomentum(...)
!*******************************************************************************
!*******************************************************************************
module HAFT_single

  implicit none
  private

  public :: setFileName
  public :: getHadesAcceptance
  public :: smearhadesmomentum, smearHades3Momentum

  character(len=200), save :: fname  = 'HadesAcceptanceFilter.acc'

  !  HAFT declaration of acceptance matrix arrays and resolution tables
  integer, parameter :: nids  = 14         ! <<== change if < max id

  type :: AcceptanceMatrix
    real(4), dimension(:), allocatable :: matrix     ! the actual matrix
    integer(4) :: xdim = 0, ydim = 0, zdim = 0       ! dimensions
    real(4) :: pmin = 0., thmin = 0., phmin = 0., &  ! limits
               pmax = 0., thmax = 0., phmax = 0.
    real :: dp = 0., dth = 0., dph = 0.              ! step sizes
  end type

  ! acceptance matrices for e+, e-, pi+, pi-, K+, K- and p
  type(AcceptanceMatrix), dimension(nids), save :: acc

  type :: ResParameterTable
    real(4), dimension(:,:), allocatable :: tab  ! the actual table
    integer(4) :: xdim, ydim                     ! dimensions
    real(4) :: pmin, pmax, thmin, thmax          ! limits
    real :: dp, dth                              ! step sizes
    logical :: resflag = .false.                 ! parameters loaded?
  end type

  ! resolution parameter tables for e+ and e-
  type(ResParameterTable), dimension(nids), save :: par

  logical :: readflag = .false.
  character(len=80) :: comment
  integer(4) :: ntab
  real(4) sigpA(3), sigpB(3), sigth, sigph, XX0  ! resolution parameters

contains

  !*****************************************************************************
  !  Sets name of input file containing the filter
  !*****************************************************************************
  subroutine setFileName(name)
      character*(*), intent(in) :: name

      integer :: dummy

      fname = name
      dummy = readHAFTmatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname |'',a80,''|'')') fname
  end subroutine setFileName


  !*****************************************************************************
  !  Opens file in unformatted direct access mode
  !  and reads HADES acceptance matrices (as linearized arrays)
  !*****************************************************************************
  integer function readHAFTmatrix()

      integer, parameter :: runit = 77  ! change if input unit is already busy

      integer(4) :: pid
      integer :: i, bins, bytes

      readHAFTmatrix = 0

      if (readflag) return

      readflag = .false.

      open(unit=runit,file=fname,access='stream',status='old',err=99)
      bytes=1
      !write(6,*) ' '
      read(runit,pos=bytes,err=100) comment
      bytes = bytes + 80
      !write(6,'(a80)') comment
      !write(6,*) '--------------------------------------'
      read(runit,pos=bytes,err=100) sigpA(1:3), sigpB(1:3), sigth, sigph, XX0
      bytes = bytes + 9*4

      do

        read(runit,pos=bytes,end=50,err=100) pid  ! break out if EOF reached
        bytes = bytes + 4

        if (pid>=0) then

          ! read acceptance matrix
          if (pid<1 .or. pid>nids) then
            write(6,*) 'acceptance not yet supported for PID ', pid, ' File = ',trim(fname)
            stop
          end if

          read(runit,pos=bytes,err=100) acc(pid)%xdim, acc(pid)%ydim, acc(pid)%zdim
          bytes = bytes + 3*4

          bins = acc(pid)%xdim * acc(pid)%ydim * acc(pid)%zdim

          read(runit,pos=bytes,err=100) acc(pid)%pmin, acc(pid)%pmax,   &
                                        acc(pid)%thmin, acc(pid)%thmax, &
                                        acc(pid)%phmin, acc(pid)%phmax
          bytes = bytes + 6*4

          allocate(acc(pid)%matrix(bins), source=0._4)
          read(runit,pos=bytes,err=100) acc(pid)%matrix(1:bins)
          bytes = bytes + bins*4

          !write(6,*) 'Acceptance matrix read for ID= ', pid
          !write(6,*) 'dims= ',acc(pid)%xdim, ' ', acc(pid)%ydim, ' ', acc(pid)%zdim
          !write(6,*) 'lims= ',acc(pid)%pmin, ' ', acc(pid)%pmax, ' ', acc(pid)%thmin, &
          !          ' ', acc(pid)%thmax, ' ', acc(pid)%phmin, ' ', acc(pid)%phmax
          !write(6,*) 'size of matrix :', bins
          !write(6,*) '--------------------------------------'
          acc(pid)%dp = (acc(pid)%pmax-acc(pid)%pmin)/real(acc(pid)%xdim)
          acc(pid)%dth = (acc(pid)%thmax-acc(pid)%thmin)/real(acc(pid)%ydim)
          acc(pid)%dph = (acc(pid)%phmax-acc(pid)%phmin)/real(acc(pid)%zdim)

        else

          ! read resolution parameters
          pid = -pid
          if (pid<1 .or. pid>nids) then
            write(6,*) 'resolution not yet supported for PID ', pid, ' File = ',trim(fname)
            stop
          end if

          read(runit,pos=bytes,err=100) par(pid)%xdim, par(pid)%ydim
          bytes = bytes + 2*4

          bins = par(pid)%xdim*par(pid)%ydim

          read(runit,pos=bytes,err=100) par(pid)%pmin, par(pid)%pmax, &
                                        par(pid)%thmin, par(pid)%thmax
          bytes = bytes + 4*4

          read(runit,pos=bytes,err=100) ntab ! nb. of parameter tables
          bytes = bytes + 4

          allocate(par(pid)%tab(ntab,bins), source=0._4)
          do i=1,ntab
            read(runit,pos=bytes,err=100) par(pid)%tab(i,1:bins)
            bytes = bytes + bins*4
          end do

          par(pid)%resflag = .true.

          !write(6,*) 'Resolution tables read for ID= ', pid
          !write(6,*) 'dims= ',par(pid)%xdim, ' ', par(pid)%ydim
          !write(6,*) 'lims= ',par(pid)%pmin, ' ', par(pid)%pmax, ' ', &
          !                    par(pid)%thmin, ' ', par(pid)%thmax
          !write(6,*) 'size of parameter tables :', ntab, ' x', bins
          !write(6,*) '--------------------------------------'
          par(pid)%dp  = (par(pid)%pmax-par(pid)%pmin)   / real(par(pid)%xdim)
          par(pid)%dth = (par(pid)%thmax-par(pid)%thmin) / real(par(pid)%ydim)

        end if

      end do

 50   close(runit)
      readHAFTmatrix = bytes-1 ! return number of bytes read
      readflag = .true.
      return

      ! Error opening or reading
 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',trim(fname)
      readHAFTMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',trim(fname)
      readHAFTmatrix = -1
      return
  end function readHAFTmatrix


  !*****************************************************************************
  !  Returns acceptance value at cell (i,j,k) of linearized matrix
  !  for particle ID
  !*****************************************************************************
  real function getMatrixVal(acc,i,j,k)
      type(AcceptanceMatrix), intent(in) :: acc
      integer, intent(in) :: i, j, k

      integer :: i1, j1, k1, ilin

      i1 = min(max(1,i),acc%xdim)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),acc%ydim)  ! This effectively extrapolates matrix
      k1 = min(max(1,k),acc%zdim)  ! beyond table boundaries.

      ilin = i1+acc%xdim*(j1-1)+acc%xdim*acc%ydim*(k1-1)  ! linearized index
      getMatrixVal = acc%matrix(ilin)
  end function getMatrixVal


  !*****************************************************************************
  !  Returns HADES acceptance for particle id of given momentum (in GeV/c),
  !  polar angle theta (in deg.) and azimuthal angle phi (in deg.)
  !  by table interpolation.
  !
  !  Lab frame used: z = beam axis, y = vertical axis
  !
  !                 ^ y
  !            x \  |
  !               \ |
  !                \|    z
  !                 O---->
  !
  !  id = 2 : positron
  !     = 3 : electron
  !     = 8 : pi+
  !     = 9 : pi-
  !     = 10: K+
  !     = 12: K-
  !     = 14: proton
  !
  !  mode = 0 : nearest-neighbour interpolation
  !       = 1 : tri-linear interpolation
  !       = 2 : tri-quadratic interpolation
  !       =-2 : tri-quadratic B-spline interpolation
  !       = 3 : tri-cubic interpolation
  !       =-3 : tri-cubic B-spline interpolation
  !       = 4 : tri-cubic Catmull-Rom spline
  !       =-4 : tri-cubic optimal cardinal spline
  !*****************************************************************************
  real function getHadesAcceptance(pid,p0,theta0,phi0,mode)
      use HAFT_aux, only: kernel

      integer, intent(in) :: pid
      real, intent(in) :: p0, theta0, phi0
      integer, intent(in), optional :: mode

      real :: p, theta, phi, u, v, w, sum_, Kx, Ky, Kz
      integer :: ix, iy, iz, i, j, k, ilo, ihi, jlo, jhi, klo, khi, mod_

      if (present(mode)) then
        mod_ = mode
      else
        mod_ = -2   ! use B-spline interpolation
      end if

      getHadesAcceptance = 0.

      if (readHAFTmatrix()==-1) return

      if (p0<acc(pid)%pmin .or. theta0<acc(pid)%thmin .or. theta0>acc(pid)%thmax .or. &
          phi0<acc(pid)%phmin) return

      p = min(p0,acc(pid)%pmax-2.01*acc(pid)%dp)  ! level off acceptance at high p
      theta = theta0
      phi = phi0
      if (phi > 60.) phi = mod(phi,60.)   ! modulo HADES sector

      ix = int(acc(pid)%xdim*((p-0.5*acc(pid)%dp-acc(pid)%pmin)/(acc(pid)%pmax-acc(pid)%pmin))) + 1  ! floor indices
      iy = int(acc(pid)%ydim*((theta-0.5*acc(pid)%dth-acc(pid)%thmin)/(acc(pid)%thmax-acc(pid)%thmin))) + 1
      iz = int(acc(pid)%zdim*((phi-0.5*acc(pid)%dph-acc(pid)%phmin)/(acc(pid)%phmax-acc(pid)%phmin))) + 1

      select case (mod_)
      case (0,1)  ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      case (2,3,4,-2,-3,-4)
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      case default  ! mode not defined
        return
      end select

      if (ilo<0 .or. jlo<0 .or. klo<0) return
      if (ihi>acc(pid)%xdim+1 .or. jhi>acc(pid)%ydim+1 .or. khi>acc(pid)%zdim+1) return

      sum_ = 0.
      do i=ilo,ihi                      ! triple interpolation loop
        u = (p - (real(i)-0.5)*acc(pid)%dp-acc(pid)%pmin)/acc(pid)%dp
        Kx = kernel(u,mod_)
        do j=jlo,jhi
          v = (theta - (real(j)-0.5)*acc(pid)%dth-acc(pid)%thmin)/acc(pid)%dth
          Ky = kernel(v,mod_)
          do k=klo,khi
            w = (phi - (real(k)-0.5)*acc(pid)%dph-acc(pid)%phmin)/acc(pid)%dph
            Kz = kernel(w,mod_)
            sum_ = sum_ + getMatrixVal(acc(pid),i,j,k)*Kx*Ky*Kz
          end do
        end do
      end do
      sum_ = max(min(sum_, 1.01),0.)  ! clip over/undershoots

      getHadesAcceptance = sum_

  end function getHadesAcceptance


  !*****************************************************************************
  !  Returns acceptance value at cell (i,j) of linearized
  !  parameter table for particle ID
  !*****************************************************************************
  real function getTableVal(par,i,j,itab)
      type(ResParameterTable), intent(in) :: par
      integer, intent(in) :: i, j, itab

      integer :: i1, j1, ilin

      i1 = min(max(1,i),par%xdim)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),par%ydim)  ! This effectively extrapolates matrix

      ilin = i1+par%xdim*(j1-1)    ! linearized index

      getTableVal = par%tab(itab,ilin)
  end function getTableVal


  !*****************************************************************************
  !     Interpolate resolution parameter table as function
  !     of momentum and theta (pin in GeV/c and theta in degree)
  !*****************************************************************************
  real function param(par,pin,thin,itab)
      use HAFT_aux, only: kernel

      type(ResParameterTable), intent(in) :: par
      real, intent(in) :: pin, thin
      integer, intent(in) :: itab

      integer :: i, j, ix, iy, ilo, ihi, jlo, jhi, mod_
      real :: p, th, sum_, u, v, Kx, Ky

      mod_ = 1
      param = 0.

      p = min(max(pin,par%pmin),par%pmax)  ! safety fence
      th = min(max(thin,par%thmin),par%thmax)

      ix = int(par%xdim*((p-0.5*par%dp-par%pmin)/(par%pmax-par%pmin))) + 1      ! floor indices
      iy = int(par%ydim*((th-0.5*par%dth-par%thmin)/(par%thmax-par%thmin))) + 1

      select case (mod_)
      case (0,1)  ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
      case (2,3,4,-2,-3,-4)
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
      case default  ! mode not defined
        return
      end select

      if (ilo<0 .or. jlo<0) return
      if (ihi>par%xdim+1 .or. jhi>par%ydim+1) return

      sum_ = 0.
      do i=ilo,ihi                      ! double interpolation loop
        u = (p - (real(i)-0.5)*par%dp-par%pmin)/par%dp
        Kx = kernel(u,mod_)
        do j=jlo,jhi
          v = (th - (real(j)-0.5)*par%dth-par%thmin)/par%dth
          Ky = kernel(v,mod_)
          sum_ = sum_ + getTableVal(par,i,j,itab)*Kx*Ky
        end do
      end do

      param = sum_
  end function param


  !*****************************************************************************
  !  Apply the Hades momentum resolution to a 4-momentum vector (in GeV/c)
  !
  !  Lab frame used: z = beam axis, y = vertical axis
  !
  !                 ^ y
  !            x \  |
  !               \ |
  !                \|    z
  !                 O---->
  !
  !  All components of the 4-momentum vector are changed by this call.
  !
  !  If parameter tables are loaded, the particle momentum is smeared according
  !  to an asymmetric response function, if not, default gaussian sampling is used.
  !
  !  The default resolution mode is determined by:
  !
  !   mode = 1 : low-resolution     (MDC 1+2)
  !        = 2 : medium-resolution  (MDC 1+2+3)
  !        = 3 : high-resolution    (MDC 1+2+3+4)
  !*****************************************************************************
  subroutine smearHades4Momentum(mom,mode,pid)
      use HAFT_aux, only: pi, sampleGauss, sampleMP

      real, intent(inout) :: mom(4)
      integer, intent(in) :: mode,pid

      integer :: i
      real :: mass, mass2, pt, pt2, ptot, ptot2, theta, phi, sinth
      real :: sigp, betainv, sigms, sigms2, sigthms, sigphms, ploss, respar(10)
      real, parameter :: r2d = 180./pi, twopi = 2.*pi

      if (.not.readflag) then
        if (readHAFTmatrix()==-1) return
      end if

      pt2 = mom(1)**2 + mom(2)**2
      pt = sqrt(pt2)
      ptot2 = pt2 + mom(3)**2
      ptot = sqrt(ptot2)                      ! total momentum in GeV/c
      mass2 = mom(4)**2 - ptot2               ! particle mass does not change
      if (mass2<0.) mass2 = 0.
      mass = sqrt(mass2)
      betainv = 1.
      if (mass2>0.) betainv = sqrt(1. + mass2/ptot2)  ! 1/beta
      sigms = 0.0136*betainv/ptot*sqrt(XX0)*(1+0.038*log(XX0))
                                              ! multiple scattering angle
      sigms2 = sigms*sigms
      sigthms = sqrt(sigth*sigth+sigms2)      ! add quadratically to resolution
      sigphms = sqrt(sigph*sigph+sigms2)

      if (pid==2 .or. pid==3) then  ! leptons, use Ar+KCl embedding values
         sigthms = 0.00246/ptot + 0.00045
         sigphms = 0.00294/ptot + 0.00059
      endif

      theta = acos(mom(3)/ptot)               ! polar angle in radian
      sinth = sin(theta)
      phi = 0.
      if (pt>0.) phi = acos(mom(1)/pt)
      if (mom(2)<0.) phi = twopi - phi     ! azimuthal angle in radian

      !  If resolution parameters are available, use dedicated smearing

      if (par(pid)%resflag) then   ! resolution parameters are loaded for pid?

!        write(6,*) mom(1), mom(2), mom(3), mass, ptot, r2d*theta
         do i=1,ntab    ! look up parameters from tables
            respar(i) = param(par(pid),ptot,r2d*theta,i)
!           write(6,*) i, respar(i)
         end do

         ptot = ptot*(1.+sampleMP(respar,2.))  ! randomize momentum

      else  ! use default gaussian smearing

         if (mode<1 .or. mode>3) return  ! unknown mode
         sigp = 0.01*ptot*(sigpA(mode)+sigpB(mode)*ptot)  ! momentum resolution
         ptot = max(0.,sampleGauss(ptot,sigp))  ! smear total momentum

         if (pid==2 .or. pid==3) then
            ploss = 0.0018
         else
            ploss = 0.0
         end if
         if (ptot>0.1) then
            if (pid==2) then
               ploss = ploss + 0.005*(1.-exp(-(ptot-0.1)/0.053)) ! positron
            else if (pid==3) then
               ploss = ploss + 0.006*(1.-exp(-(ptot-0.1)/0.078)) ! electron
            end if
         end if
         ptot = ptot - ploss                    ! Eloss of electrons in target

      end if

      theta = abs(sampleGauss(theta,sigthms)) ! smear polar angle
      if (theta>pi) theta = twopi - theta  ! if > pi, mirror angle
      if (sinth>0.01) phi = sampleGauss(phi,sigphms/sinth) ! smear azimuth
      if (phi<0.) phi = phi + twopi        ! and check if within range
      if (phi>twopi) phi = phi - twopi

      sinth = sin(theta)
      mom(1) = ptot*sinth*cos(phi)            ! new momentum components
      mom(2) = ptot*sinth*sin(phi)
      mom(3) = ptot*cos(theta)
      mom(4) = sqrt(sum(mom(1:3)**2) + mass2)  ! total energy

  end subroutine smearHades4Momentum


  subroutine smearHadesMomentum(p,mode,pid)
      real(8), intent(inout) :: p(0:3)
      integer, intent(in) :: mode, pid

      real :: mom4(4)

      mom4(1:3) = p(1:3)
      mom4(4) = p(0)

      call smearHades4Momentum(mom4,mode,pid)

      p(1:3) = mom4(1:3)
      p(0) = mom4(4)
  end subroutine smearHadesMomentum


  !*****************************************************************************
  !  Apply Hades momentum resolution to a 3-momentum (calculate multiple
  !  scattering assuming the particle is an electron)
  !*****************************************************************************
  subroutine smearHades3Momentum(mom3,mode,pid)
      real, intent(inout) :: mom3(3)
      integer, intent(in) :: mode,pid

      real :: mom4(4), mass

      if (mode<1 .or. mode>3) return    ! unknown mode

      if (pid==2 .or. pid==3) then        ! e+ or e-
         mass = 0.000510998918
      else if (pid==8 .or. pid==9) then   ! pi+ or pi-
         mass = 0.13957018
      else if (pid==11 .or. pid==12) then ! K+ or K-
         mass = 0.493677
      else if (pid==14) then                ! proton
         mass = 0.938272029
      else
         return  ! particle is not supported
      end if

      mom4(1:3) = mom3(1:3)
      mom4(4) = sqrt(sum(mom3(:)**2) + mass**2)
      call smearHades4Momentum(mom4,mode,pid)
      mom3(1:3) = mom4(1:3)
  end subroutine smearHades3Momentum


end module HAFT_single



!*******************************************************************************
!*******************************************************************************
! This modules contains routines for pair acceptance filtering.
!
!  Usage : 1) set acceptance file name with
!                 call setPairFileName(fname)
!          2) sample pair acceptance values with calls to
!                 acc = getHadesPairAcceptance(mass,pt,rapidity,mode)
!          3) set resolution parameters via
!                 call setResolutionParameters(...)
!          4) apply detector resolution via
!                 call smearHadesPair(...)
!*******************************************************************************
!*******************************************************************************
module HAFT_pair

  implicit none
  private

  public :: setPairFileName
  public :: getHadesPairAcceptance
  public :: setResolutionParameters
  public :: smearHadesPair

  character(len=200), save :: fname = 'HadesPairAcceptanceFilter.acc'

  !  HAFT declaration of acceptance matrix arrays
  real(4), dimension(:), allocatable :: matrix51

  character(len=80) :: comment
  logical :: readflag = .false.
  integer(4) :: xdim, ydim, zdim
  real(4) :: mmin, mmax, ptmin, ptmax, rapmin, rapmax
  real :: dm, dpt, drap
  real :: sigpA(3), sigpB(3)  ! resolution parameters

contains

  !*****************************************************************************
  !  Returns HADES pair acceptance for given mass (in GeV/c**2),
  !  transverse momentum (in GeV/c) and rapidity (in lab frame)
  !  by table interpolation.
  !
  !  Lab frame used: z = beam axis, y = vertical axis
  !
  !                 ^ y
  !            x \  |
  !               \ |
  !                \|    z
  !                 O---->
  !
  !
  !  mode = 0 : nearest-neighbour interpolation
  !       = 1 : tri-linear interpolation
  !       = 2 : tri-quadratic interpolation
  !       =-2 : tri-quadratic B-spline interpolation
  !       = 3 : tri-cubic interpolation
  !       =-3 : tri-cubic B-spline interpolation
  !       = 4 : tri-cubic Catmull-Rom spline
  !       =-4 : tri-cubic optimal cardinal spline
  !*****************************************************************************
  real function getHadesPairAcceptance(mass0,pt0,rap0,mode)
      use HAFT_aux, only: kernel

      real, intent(in) :: mass0, pt0, rap0
      integer, intent(in), optional :: mode

      real :: mass, pt, rap, u, v, w, sum_, Kx, Ky, Kz
      integer :: ix, iy, iz, i, j, k, ilo, ihi, jlo, jhi, klo, khi, mod_

      getHadesPairAcceptance = 0.

      if (readHAFTPairMatrix()==-1) return

      if (present(mode)) then
        mod_ = mode  ! (use mode = 0 or 1, otherwise problems at pt=0!)
      else
        mod_ = 1    ! use trilinear interpolation
      end if

      if (mass0<mmin .or. pt0<ptmin .or. pt0>ptmax .or. rap0<rapmin .or. rap0>rapmax) &
        return

      mass = min(mass0,mmax-2.01*dm)  ! level off acceptance at high mass
      pt = pt0
      rap = rap0

      ix = int(xdim*((mass-0.5*dm-mmin)/(mmax-mmin))) + 1      ! floor indices
      iy = int(ydim*((pt-0.5*dpt-ptmin)/(ptmax-ptmin))) + 1
      iz = int(zdim*((rap-0.5*drap-rapmin)/(rapmax-rapmin))) + 1

      select case (mod_)
      case (0,1)  ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      case (2,3,4,-2,-3,-4)
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      case default  ! mode not defined
        return
      end select

      if (ilo<0 .or. jlo<0 .or. klo<0) return
      if (ihi>xdim+1 .or. jhi>ydim+1 .or. khi>zdim+1) return

      sum_ = 0.
      do i=ilo,ihi                      ! triple interpolation loop
        u = (mass - (real(i)-0.5)*dm-mmin)/dm
        Kx = kernel(u,mod_)
        do j=jlo,jhi
          v = (pt - (real(j)-0.5)*dpt-ptmin)/dpt
          Ky = kernel(v,mod_)
          do k=klo,khi
            w = (rap - (real(k)-0.5)*drap-rapmin)/drap
            Kz = kernel(w,mod_)
            sum_ = sum_ + getMatrixVal(i,j,k)*Kx*Ky*Kz
          end do
        end do
      end do
      sum_ = max(min(sum_, 1.01),0.)  ! clip over/undershoots

      getHadesPairAcceptance = sum_

  end function getHadesPairAcceptance


  !*****************************************************************************
  !  Opens file in unformatted direct access mode
  !  and reads HADES pair acceptance matrix (as linearized array)
  !*****************************************************************************
  integer function readHAFTPairMatrix()

      integer, parameter :: runit = 78  ! change if input unit is already busy

      integer :: bins, bytes

      readHAFTPairMatrix = 0

      if (readflag) return

      readflag = .false.

      open(unit=runit,file=fname,access='stream',status='old',err=99)
      bytes=1
      read(runit,pos=bytes,err=100) comment
      write(6,'(a80)') comment
      bytes = bytes + 80
      read(runit,pos=bytes,err=100) xdim, ydim, zdim

      bins = xdim*ydim*zdim
      allocate(matrix51(bins), source=0._4)

      bytes = bytes + 3*4
      read(runit,pos=bytes,err=100) mmin,mmax,ptmin,ptmax,rapmin,rapmax
      bytes = bytes + 6*4
      read(runit,pos=bytes,err=100) matrix51(1:bins)
      write(6,'(''Matrix read for e+e- pairs'')')
      bytes = bytes + bins*4
      close(runit)

      dm = (mmax-mmin)/real(xdim)
      dpt = (ptmax-ptmin)/real(ydim)
      drap = (rapmax-rapmin)/real(zdim)

      write(6,'('' coms= '',a80)') comment
      write(6,*) 'dims= ',xdim, ' ', ydim, ' ', zdim
      write(6,*) 'lims= ',mmin, ' ', mmax, ' ', ptmin, ' ', ptmax, &
                 ' ', rapmin, ' ', rapmax
      write(6,*) 'size of matrix :', bins

      readHAFTPairMatrix = bytes-1 ! return number of bytes read
      readflag = .true.
      return

      ! Error opening or reading
 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',trim(fname)
      readHAFTPairMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',trim(fname)
      readHAFTPairMatrix = -1
      return
  end function readHAFTPairMatrix


  !*****************************************************************************
  !  Sets name of input file containing the filter
  !*****************************************************************************
  subroutine setPairFileName(name)
      character*(*), intent(in) :: name

      integer :: dummy

      fname = name
      dummy = readHAFTPairMatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname |'',a80,''|'')') fname
  end subroutine setPairFileName


  !*****************************************************************************
  !  Returns acceptance value at cell (i,j,k) of linearized matrix
  !  for particle ID
  !*****************************************************************************
  real function getMatrixVal(i,j,k)
      integer, intent(in) :: i, j, k

      integer :: i1, j1, k1, ilin

      i1 = min(max(1,i),xdim)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),ydim)  ! This effectively extrapolates matrix
      k1 = min(max(1,k),zdim)  ! beyond table boundaries.

      ilin = i1+xdim*(j1-1)+xdim*ydim*(k1-1)  ! linearized index

      getMatrixVal = matrix51(ilin)

  end function getMatrixVal


  !*****************************************************************************
  !  Apply Hades momentum resolution to a pair (calculate multiple
  !  scattering assuming the particle is an electron)
  !*****************************************************************************
  subroutine smearHadesPair(pair,mode)
      use HAFT_aux, only: interpol, sampleGauss, sampleMP

      real, intent(inout) :: pair(3)
      integer, intent(in) :: mode

      real :: m, pt, rap, sigpt, sigm, sigrap, respar(10)

      real, parameter :: mtab(10) = (/ 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2 /)  ! mass grid
      real, parameter :: par1(10) = (/ 0.0077, -0.0082, -0.0125, -0.0120, -0.0114, &
                                        -0.0106, -0.0098, -0.0085, -0.0078, -0.0075 /)
      real, parameter :: par2(10) = (/ 0.0820, 0.0460, 0.0260, 0.0210, 0.0190, &
                                         0.0183, 0.0182, 0.0181, 0.0180, 0.0180 /)
      real, parameter :: par3(10) = (/ 13.5, 16.2, 19.9, 20.2, 19.2, &
                                         18.0, 16.9, 14.8, 12.8, 11.0 /)
      real, parameter :: par4(10) = (/-10.0, -15.4, -20.4, -23.1, -22.3, &
                                        -21.6, -20.6, -19.4, -18.2, -17.9 /)
      real, parameter :: par5(10) = (/ 18.1, 11.6, 10.4, 10.0, 9.4, 8.5, 7.8, 7.0, 6.2, 5.7 /)

      if (.not.readflag) then
        if (readHAFTPairMatrix()==-1) return
      end if

      m = pair(1)
      pt = pair(2)
      rap = pair(3)

      if (mode==4) then  ! use skewd mass function
        sigpt = 0.01*pt*(sigpA(3)+sigpB(3)*pt)/sqrt(2.)       ! pt resolution
        pt = max(0.,sampleGauss(pt,sigpt))                  ! smear pt

        respar(1) = interpol(m,mtab,par1,10)
        respar(2) = interpol(m,mtab,par2,10)
        respar(3) = interpol(m,mtab,par3,10)
        respar(4) = interpol(m,mtab,par4,10)
        respar(5) = interpol(m,mtab,par5,10)
!        do i=1,5    ! interpolate parameters
!          write(6,*) i, respar(i)
!        end do
        m = m*(1.+sampleMP(respar,sqrt(2.)))                  ! randomize mass
      else ! use gaussian smearing
        sigpt = 0.01*pt*(sigpA(mode)+sigpB(mode)*pt)/sqrt(2.) ! pt resolution
        pt = max(0.,sampleGauss(pt,sigpt))                  ! smear pt

        sigm = 0.01*m*(sigpA(mode)+sigpB(mode)*m)/sqrt(2.)
        m = max(0.,sampleGauss(m,sigm))                     ! smear mass
      end if

      sigrap = 0.1    !  just a guess!
      rap = sampleGauss(rap,sigrap)

      pair(1) = m
      pair(2) = pt
      pair(3) = rap
  end subroutine smearHadesPair


  !*****************************************************************************
  !     Set momentum resolution parameters
  !*****************************************************************************
  subroutine setResolutionParameters(mode,a,b)
      real, intent(in) :: a, b
      integer, intent(in) :: mode

      if (mode>=1 .and. mode<=3) then
         sigpA(mode) = a
         sigpB(mode) = b
      endif

      write(6,*) 'mode=1   ', ' ', sigpA(1), ' ', sigpB(1)
      write(6,*) 'mode=2   ', ' ', sigpA(2), ' ', sigpB(2)
      write(6,*) 'mode=3   ', ' ', sigpA(3), ' ', sigpB(3)
  end subroutine setResolutionParameters


end module HAFT_pair
