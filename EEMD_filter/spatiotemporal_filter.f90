module spatiotemporal_filter
      use omp_lib
!
      IMPLICIT NONE
      contains
      subroutine filter(gp,gpf,ni,nj,nij,ntime,nmodes)
!
      integer, parameter ::ikind=selected_real_kind(8)
      integer ij,nall,ntime,ni,nj,i,j,it,n,nvar,nk,n1,nij,npc,m
      integer nmodes,nyear,ndays,nn,ii
      real(kind=ikind) :: threshold1,threshold2,corr
      parameter (npc=12) ! number of pc's and eof's desired
      parameter (threshold1=1.0e-20_ikind) ! rmse threshold for pc iteration
      parameter (threshold2=0.999995_ikind)! correlation thresh for pc iter
!
      real gp(ni,nj,ntime),gpf(ni,nj,ntime,nmodes)
      real (kind=ikind) :: wt(nij)
      real sst(ni,nj),ssto(ni,nj),pcn(ntime,npc),pcf(ntime,nmodes,npc)
      real (kind=ikind) :: t1, t2, t3, t4, t5, t6, t7, t8,t9,t10,t11,t12
!
      real(kind=ikind) , allocatable :: reconst(:,:,:,:)
      real(kind=ikind) , allocatable :: eof(:,:)
      real(kind=ikind) , allocatable :: eofall(:,:)
      real(kind=ikind) , allocatable :: tmean(:)
      real(kind=ikind) , allocatable :: orig(:,:)
      real(kind=ikind) , allocatable :: origo(:,:)
      real(kind=ikind) , allocatable :: origo2(:,:)
      real(kind=ikind) , allocatable :: pc(:,:)
      real(kind=ikind) , allocatable :: pcall(:,:)
      real(kind=ikind) , allocatable :: pc2(:,:)
!
      allocate (reconst(ni,nj,ntime,npc)) ! reconstructed data from eof's
      allocate (eof(1,nij))         ! dummy's for eof's
      allocate (eofall(npc,nij))    ! contains all eof's!
      allocate (tmean(nij))    ! contains all eof's!
      allocate (orig(ntime,nij))   ! dummy's for reconstruction
      allocate (origo(ntime,nij))   ! dummy's for reconstruction
      allocate (origo2(ntime,nij))  ! dummy's for reconstruction
      allocate (pc(ntime,1))        ! dummy's for pc's
      allocate (pc2(ntime,1))       ! dummy's for pc's
      allocate (pcall(ntime,npc))   ! contains all pc
!      
!      
!
!
      do 300 n=1,ntime
       ij=0
       do i=1,ni
        do j=1,nj
         if(gp(i,j,n).ne.-9999.0)then
          ij=ij+1
          origo(n,ij)=gp(i,j,n)
         endif
        enddo
       enddo
300    continue

       call removetimemean(origo,nij,ntime,tmean)
       orig=origo
!
! compute first and successive eof & pc iteratively
!
       wt=1.0_ikind
       call eofpc(origo,pc,pc2,eof,nij,ntime,threshold1,&
       threshold2,wt,1)
      t7=omp_get_wtime()
!
! Original data being computed from first eof and pc
!
       origo=matmul(pc2,eof)
      t8=omp_get_wtime()
!
! first PC written to a consolidated array for all pc's to write
! to file finally
!
       do it=1,ntime
        pcn(it,1)=pc2(it,1)
       end do
!
!
! storing the first eof and pc in a consolidated array for all eof's and
! pc's
!
       pcall(:,1)=pc2(:,1)
       eofall(1,:)=eof(1,:)
!
! if only one EOF is being computed
!
       if(npc.eq.1)go to 500
!
! Now on to the EOF's and pc's beyond first
!
       do n=2,npc
!
! preparing input data for computing next EOF's and PC's
!
        origo2=origo
        orig=orig-origo2
        origo=orig
        call eofpc(origo,pc,pc2,eof,nij,ntime,threshold1,&
        threshold2,wt,n)
        pcall(:,n)=pc2(:,1)
        eofall(n,:)=eof(1,:)
        pcn(:,n)=pc2(:,1)
!
! Original data being computed from higer order eof's and pc's
!
        origo=matmul(pc2,eof)
       end do
500   continue
      print*, 'done with eof '
!
! filtering pc using eemd
!
      do n=1,npc
      m=3571
      call eemd (ntime,pcn(:,n),0.2,250,nmodes,m,pcf(:,:,n))
      end do
      print*, 'done with eemd'
!
      print*, 'nmodes , npc, ntime ', nmodes, npc, ntime
      do m=1,nmodes
       do n=1,npc
        pc2(:,1)=pcf(:,m,n)
        eof(1,:)=eofall(n,:)
        origo=matmul(pc2,eof)
        do it=1,ntime
         ij=0
         do i=1,ni
          do j=1,nj
           if(gp(i,j,it).ne.-9999.0_ikind)then
            ij=ij+1
            reconst(i,j,it,n)=origo(it,ij)
           else
            reconst(i,j,it,n)=-9999.0_ikind
           endif
          enddo
         enddo
        end do ! it loop
       end do  ! n loop
       do it=1,ntime
        ssto=0.0
        do n=1,npc
         ssto(:,:)=ssto(:,:)+reconst(:,:,it,n)
        end do
        gpf(:,:,it,m)=ssto(:,:)
       end do  ! it loop
      end do  ! m loop
      print*, 'out of filter'

      end subroutine filter 
!
     subroutine removetimemean(difft,nx,ntime,tmean)
     use omp_lib
     implicit none
     integer, parameter ::ikind=selected_real_kind(8)
     integer, parameter :: nthreads = 32
     integer,intent (in):: nx,ntime
     integer ix,it
     real(kind=ikind),intent(inout) :: difft(ntime,nx)
     real(kind=ikind),intent(out) :: tmean(nx)
     real(kind=ikind) :: sum,dev(nx)
     call omp_set_num_threads(nthreads)

       !$omp parallel shared(difft,tmean) private(ix,it,sum)
       !$omp do schedule(static)

       do ix=1,nx
        sum=0.0_ikind
        do it=1,ntime
         sum=sum+difft(it,ix)
        end do
        tmean(ix)=sum/float(ntime)
      end do
      !$omp end do
      !$omp end parallel
!
       !$omp parallel shared(difft,tmean) private(ix,it)
       !$omp do schedule(static)
       do it=1,ntime
        do ix=1,nx
         difft(it,ix)=difft(it,ix)-tmean(ix)
        end do
       end do
       !$omp end do
       !$omp end parallel
!
       !$omp parallel shared(difft,dev) private(ix,it)
       !$omp do schedule(static)
       do ix=1,nx
        dev(ix)=0.0_ikind
        do it=1,ntime
         dev(ix)=dev(ix)+difft(it,ix)**2
        end do
        dev(ix)=sqrt(dev(ix)/float(ntime))
       end do
       !$omp end do
       !$omp end parallel
       return
       end subroutine removetimemean
!
        subroutine eofpc(origo,pc,pc2,eof,nij,ntime,threshold1,&
                         threshold2,wt,n)
        use omp_lib
        implicit none
        integer, parameter ::ikind=selected_real_kind(8)
        integer,intent(in) :: nij,ntime,n
        integer ij
        real(kind=ikind),intent(in) :: wt(nij),threshold1,threshold2
        real(kind=ikind),intent(inout) :: eof(1,nij)
        real(kind=ikind),intent(inout) :: origo(ntime,nij)
        real(kind=ikind),intent(inout) :: pc(ntime,1),pc2(ntime,1)
        real(kind=ikind) :: corr,tol,dmean(ntime)
        real(kind=ikind) :: s1,s2,s3,s4,s5,s6,s7,s8,s9
!
! Doing the first iteration with pc acquiring the domain average value
!
        s1=omp_get_wtime()
        call domainavg(origo,nij,ntime,dmean,n)
        pc(:,1)=dmean(:)
!
! computing eof with first guess of pc
!
        call geteof(pc,origo,eof,nij,ntime,n)
!
! computing the next iteration of pc
!
        call getpc(origo,eof,pc2,wt,nij,ntime)
!
! checking for rmse and tcor of pc between successive iterations
!
        call difference(pc, pc2,tol,corr,ntime)
        ij=1
!
! skipping iterations if thresholds are met
!
        if (tol.lt.threshold1 .and. corr.ge.threshold2) go to 200
!
! on to successive iterations
!
100     continue
        ij=ij+1
        pc=pc2
!
! next iteration of eof
!
        call geteof(pc,origo,eof,nij,ntime,n)
!
! recomputing pc
!
        call getpc(origo,eof,pc2,wt,nij,ntime)
!
! checking for tolerance of pc
!
        call difference(pc, pc2, tol,corr,ntime)
        s9=omp_get_wtime()
!
! if tolerance check fails go on to next iteration
!
        if (tol.gt.threshold1 .and. corr.lt.threshold2) go to 100
200     continue
        print '(a,i3,a,i4,a,e13.6,a,e13.6)','eof#=',n, '  iteration=',ij, &
                '  tolerance=',tol,'  correlation= ',corr
        print '(a, f15.3,a)', 'eofpc timing= ', s9-s1,' s'
        return
        end subroutine eofpc
!
! subroutine to compute difference in pc's between successive iterations
!
      subroutine difference(difft, difftT,tol,corr,ntime)
      implicit none
      integer, parameter ::ikind=selected_real_kind(8)
      integer, intent(in):: ntime
      integer in,it
      real(kind=ikind), intent(in) :: difft(ntime),difftT(ntime)
      real(kind=ikind),intent(out) :: tol,corr
      real(kind=ikind) :: r1,r2,d1,d2
      tol=0.0_ikind
      corr=0.0_ikind
      d1=0.0_ikind
      d2=0.0_ikind
      r1=0.0_ikind
      r2=0.0_ikind
      do it=1,ntime
       r1=r1+difft(it)
       r2=r2+difftT(it)
      end do
      r1=r1/float(ntime)
      r2=r2/float(ntime)
      in=0
      do it=1,ntime
       d1=d1+(difft(it)-r1)**2
       d2=d2+(difftT(it)-r2)**2
       tol=tol+(difft(it)-difftT(it))**2
       corr=corr+(difft(it)-r1)*(difftT(it)-r2)
       in=in+1
      end do
      tol=sqrt(tol)/float(ntime)
      corr=corr/sqrt(d1*d2)
      return
      end subroutine difference
!
! subroutine to compute pc given eof and original data
!
      subroutine getpc(difftT,eof,pc,wt,nx,ntime)
      use omp_lib
      implicit none
      integer, parameter ::ikind=selected_real_kind(8)
      integer, parameter :: nthreads = 32
      integer,intent (in) :: ntime,nx
      integer it,i
      real(kind=ikind),intent(out) :: pc(ntime)
      real(kind=ikind) :: num,denom
      real(kind=ikind),intent(in) :: wt(nx),difftT(ntime,nx),eof(nx)
      call omp_set_num_threads(nthreads)

      !$omp parallel shared(difftT,eof,pc,wt) private(i,it,num,denom)
      !$omp do schedule(static)
      do it=1,ntime
       num=0.0_ikind
       denom=0.0_ikind
       do i=1,nx
        num=num+wt(i)*eof(i)*difftT(it,i)
        denom=denom+wt(i)*eof(i)*eof(i)
       end do
       pc(it)=num/denom
      end do
      !$omp end do
      !$omp end parallel
      return
      end subroutine getpc
!
! subroutine to compute eof given original data and pc
!
      subroutine geteof(pc,difft,eof,nx,ntime,n)
      use omp_lib
      implicit none
      integer, parameter ::ikind=selected_real_kind(8)
      integer, parameter :: nthreads = 32
      integer,intent(in) :: nx,ntime,n
      integer it,i
      real(kind=ikind),intent(in) :: difft(ntime,nx),pc(ntime)
      real(kind=ikind), intent(out) :: eof(nx)
      real(kind=ikind) :: pc2,num
      call omp_set_num_threads(nthreads)

       !$omp parallel shared(eof,pc,difft) private(i,num,pc2)
       !$omp do schedule(static)

      do i=1,nx
       num=0.0_ikind
       pc2=0.0_ikind
       do it=1,ntime
        num=num+pc(it)*difft(it,i)
        pc2=pc2+pc(it)*pc(it)
       end do
       eof(i)=num/pc2
      end do
      !$omp end do
      !$omp end parallel
      return
      end subroutine geteof
!
! subroutine to compute domain mean average for first guess of pc
!
      subroutine domainavg(difft,nx,ntime,mean,n)
      implicit none
      integer, parameter ::ikind=selected_real_kind(8)
      integer, intent(in):: nx,ntime,n
      integer ix,it
      real(kind=ikind),intent(in) :: difft(ntime,nx)
      real(kind=ikind),intent(out) :: mean(ntime)
      real(kind=ikind) :: sum
!!
      do it=1,ntime
       sum=0.0_ikind
       do ix=1,nx
        sum=sum+difft(it,ix)
       end do
       sum=sum/float(nx)
       mean(it)=sum
      end do
      return
      end subroutine domainavg
!
      function ran2(idum)
!---------------------------------------------------------------------
!  See NUMERICAL RECIPES for detail of the function under the same
!  name
!  PARAMETERS:
!       idum       : random seed
!---------------------------------------------------------------------
      implicit none

      real :: ran2
      integer, intent(inout) :: idum
      integer, parameter :: im1=2147483563, im2=2147483399, imm1=im1-1
      integer, parameter :: ia1=40014, ia2=40692, iq1=53668, iq2=52774
      integer, parameter :: ir1=12211, ir2=3791, ntab=32, ndiv=1+imm1/ntab
      real, parameter :: am=1./im1, eps=1.2e-7, rnmx=1.-eps
      integer :: idum2=123456789, iy=0, j, k
      integer, dimension(ntab) :: iv
      
      save iv,iy,idum2
      
      do j=1,ntab
        iv(j)=0
      enddo

      if(idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if(idum.lt.0) idum=idum+im1
          if(j.le.ntab) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
      end function ran2
!
!
      FUNCTION GASDEV(IDUM)
!----------------------------------------------------------------------
!  Gaussian white noise generator using uniformly distributed white
!  noise generator
!---------------------------------------------------------------------- 

      implicit none
      real:: gasdev
!      real:: ran2
      integer, intent(inout):: idum
      real :: v1, v2, gv1, gv2
      real, parameter :: pi=3.1415927
      V1=RAN2(IDUM)
      V2=RAN2(IDUM)

      gv1=sqrt(-2.0*log(v1))*cos(2.0*pi*v2)
      gv2=sqrt(-2.0*log(v2))*cos(2.0*pi*v1)

      gasdev=gv1

      END FUNCTION GASDEV
!
      subroutine eemd(LXY,indata,An,Nesb,Nm,idum,rslt)
!-----------------------------------------------------------------------------
! This is a subroutine to decompose data(LXY) in terms of its EEMD component.
! When the added noise amplitude is zero (An=0) and the ensemble number
! is one (Nesb=0), the code degenerates to standard EMD
!
! In this code, the number of the oscillatory components is specified as an
! input, Nm. For most cases, automatic calculation, Nm can be set to log2(LXY)-2
!
! INPUT DATA:
!		indata(LXY)		input array with a length of LXY
!  		An 			amplitude of added noise
!  		Nesb 			number of ensemble members
!  		Nm 			number of mode in each decomposition
!  		idum 			seed for the random number
! OUTPUT DATA:
!  		rslt(LXY, Nm+2) 	output rseults, which contains Nm+2 columns
!					The first column is the origninal indata;
!					The second column is the first component;
!					The third column is the second component;
!					and so on
!					The last column is the residue
!-------------------------------------------------------------------------------
      implicit none
      INTEGER, PARAMETER :: MAXSIZE=1000000
!      INCLUDE 'eemd.common'
      real :: An				! noise amplitude
      integer, intent(in) :: LXY, Nesb, Nm
      integer, intent(inout):: idum
      real, dimension(LXY), intent(in) :: indata
      real, dimension(LXY) :: ximf		! data for sifting
      real, dimension(LXY) :: spmax		! upper envelope, cubic spline
      real, dimension(LXY) :: spmin		! lower envelope, cubic spline
      real, dimension(LXY) :: ave		! the average of spmax and spmin
      real, dimension(LXY) :: remp		! the input data for sifting
      real, dimension(LXY) :: rem		! the remainder (remp-ximf)

      real, dimension(LXY,Nm) :: rslt		! the final output data
      real, dimension(LXY,Nm) :: allmode	! the results of a single EMD decomposition

      integer :: nmax, nmin			! numbers of maximum and minimum
      real, dimension(LXY) :: trend		! linear trend of indata
      real :: std		! standard deviation of the linearly detrended indata
      integer:: i,j,IE,im,ii,Nm1
      real :: fNesb,  gs
!      real :: gasdev


!  initialize the output
      rslt=0.0

!
!  ensemble EMD
!
!  *******************************************************
      do IE=1,Nesb
!  *******************************************************
!
        call standev(LXY,indata,trend, std)


!  inputted data + noise
        if(Nesb.eq.1) then
          ximf=indata
        else
          do i=1,LXY
            gs =gasdev(idum)
            ximf(i)=indata(i)+An*std*gs
          enddo
        endif

!  save modified data
        do i=1,LXY
          allmode(i,1)=ximf(i)
        enddo

!  calculate modes
        Nm1=Nm-2
!       =======================================================
        do im=1,Nm1
!       =======================================================
!
!  leave a copy of the input data before IMF is calculated
          remp=ximf
!
!  Sifting
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          do ii=1,10
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            call min_max(LXY,ximf,spmax,spmin,nmax,nmin)
            call natural_spline(spmax,LXY,nmax)
            call natural_spline(spmin,LXY,nmin)

            ave=(spmax+spmin)/2.0
            ximf=ximf-ave

!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

          rem=remp-ximf

          do i=1,LXY
            allmode(i,im+1)=ximf(i)
          enddo

          ximf=rem

!       =======================================================
        enddo
!       =======================================================

        do i=1,LXY
          allmode(i,Nm)=ximf(i)
        enddo

        do j=1,Nm
          do i=1,LXY
            rslt(i,j)=rslt(i,j)+allmode(i,j)
          enddo
        enddo

!     ---------------------------------------------------------
      enddo
!     ---------------------------------------------------------
      
      fNesb=real(Nesb)
      rslt=rslt/fNesb

      end subroutine eemd
!
!**********************************************************************
!
      subroutine min_max(LEX,ximf,spmax,spmin,nmax,nmin)
!------------------------------------------------------------------
!  This is a routine to define maxima and minima from series ximf.
!  All the extrema are defined as the corresponding values of
!  ximf in spmax and spmin. All non-extrema values in spmax and
!  spmin are defined as 1.0e31.
!------------------------------------------------------------------
      implicit none

      integer, intent(in):: LEX
      real, dimension(LEX), intent(in):: ximf
      real, dimension(LEX), intent(out):: spmax, spmin
      integer, intent(out):: nmax, nmin
      integer:: i
      
      nmax=0
      nmin=0

      spmax(1)=ximf(1)
      spmax(LEX)=ximf(LEX)
      spmin(1)=spmax(1)
      spmin(LEX)=spmax(LEX)

      nmax=2
      nmin=2

      do i=2,LEX-1
        if(ximf(i) > ximf(i-1) .and. ximf(i) >= ximf(i+1)) then
          spmax(i) = ximf(i)
          nmax = nmax+1
        else
          spmax(i)=1.0e31
        endif
        if(ximf(i) < ximf(i-1) .and. ximf(i) <= ximf(i+1)) then
          spmin(i)=ximf(i)
          nmin=nmin+1
        else
          spmin(i)=1.0e31
        endif
      enddo

      call endmax(LEX, spmax, nmax)
      call endmin(LEX, spmin, nmin)

      end subroutine min_max

!****************************************************************

      subroutine endmax(LEX, temp, nmax)
!--------------------------------------------------------------------
! This is a subroutine to determine end values of the upper envolope
! using the method described in Appendix B of Wu and Huang (2009, 
! AADA, Vol. 1, pp1).
!--------------------------------------------------------------------
      implicit none

      integer, intent(in) :: nmax, LEX
      real, dimension(LEX), intent(inout):: temp 
      real, dimension(nmax) :: exmax, X
      integer :: I, J, lend
      real :: slope1, slope2, tmp1, tmp2
     
      lend=nmax

      J=1
      DO I=1, LEX
        IF( temp(I).LT.1.0E30 ) THEN
          X(J)=float(I)
          exmax(J)=temp(I)
          J=J+1
        ENDIF
      ENDDO

      if (nmax >= 4) then
        slope1=(exmax(2)-exmax(3))/(X(2)-X(3))
        tmp1=slope1*(X(1)-X(2))+exmax(2)
        if(tmp1 > exmax(1)) then
          temp(1)=tmp1
        endif
        
        slope2=(exmax(lend-1)-exmax(lend-2))/(X(lend-1)-X(lend-2))
        tmp2=slope2*(X(lend)-X(lend-1))+exmax(lend-1)
        if(tmp2 > exmax(lend)) then
          temp(LEX)=tmp2
        endif
      endif
      
      end subroutine endmax


!****************************************************************

      subroutine endmin(LEX, temp, nmax)
!--------------------------------------------------------------------
! This is a subroutine to determine end values of the lower envolope
! using the method described in Appendix B of Wu and Huang (2009, 
! AADA, Vol. 1, pp1).
!--------------------------------------------------------------------
      implicit none

      integer, intent(in) :: nmax, LEX
      real, dimension(LEX), intent(inout):: temp 
      real, dimension(nmax) :: exmax, X
      integer :: I, J, lend
      real :: slope1, slope2, tmp1, tmp2
     
      lend=nmax

      J=1
      DO I=1, LEX
        IF( temp(I).LT.1.0E30 ) THEN
          X(J)=float(I)
          exmax(J)=temp(I)
          J=J+1
        ENDIF
      ENDDO

      if (nmax >= 4) then
        slope1=(exmax(2)-exmax(3))/(X(2)-X(3))
        tmp1=slope1*(X(1)-X(2))+exmax(2)
        if(tmp1 < exmax(1)) then
          temp(1)=tmp1
        endif
        
        slope2=(exmax(lend-1)-exmax(lend-2))/(X(lend-1)-X(lend-2))
        tmp2=slope2*(X(lend)-X(lend-1))+exmax(lend-1)
        if(tmp2 < exmax(lend)) then
          temp(LEX)=tmp2
        endif
      endif
      
      end subroutine endmin



      SUBROUTINE NATURAL_SPLINE(YA,LEX,N)
!----------------------------------------------------------------------
!  This is a program of cubic spline interpolation. The imported
!  series, YA have a length of LEX, with N numbers of value
!  not equal to 1.0E31. The program is to use the cubic line to
!  interpolate the values for the points other thatn these N
!  numbers.        
!-----------------------------------------------------------------------

      implicit none

      integer, intent(in):: N, LEX      
      real, dimension(LEX), intent(inout):: YA
      real, dimension(N):: Y, Y2
      integer, dimension(N):: LX
      integer:: J, I, KLO, KHI
      real:: YP1, YPN, H, A, B 


!  The following code is to realocate the series of X(N), Y(N)

      J=1
      DO I=1, LEX
        IF( YA(I).LT.1.0E30 ) THEN
          LX(J)=I
          Y(J)=YA(I)
          J=J+1
        ENDIF
      ENDDO

!  The following code is used to calculate the second order derivative, 
!  set the derivatives at both ends for a natural cubic spline.
      YP1=1.0E31
      YPN=1.0E31

      CALL SPLINE(LX,Y,N,YP1,YPN,Y2)

!  calculate the cubic spline
      DO I=2,N
        KLO=LX(I-1)
        KHI=LX(I)
        H=real(KHI-KLO)
        DO J=KLO+1,KHI-1
          A=REAL(KHI-J)/H
          B=REAL(J-KLO)/H
          YA(J)=A*Y(I-1)+B*Y(I)   &
                +((A*A*A-A)*Y2(I-1)+(B*B*B-B)*Y2(I))*(H*H)/6.0
        ENDDO
      ENDDO

      END SUBROUTINE NATURAL_SPLINE


!**************************************************************************


      SUBROUTINE SPLINE(LX,Y,N,YP1,YPN,Y2)
!-----------------------------------------------------------------
!  see NUMEIRCAL RECIPES to find out meaning of each variables
!-----------------------------------------------------------------
      implicit none
    
      integer, intent(in) :: N
      real, intent(in) :: YP1, YPN 
      integer, dimension(N), intent(in):: LX
      real, dimension(N) :: X, U
      real, dimension(N), intent(inout):: Y, Y2 
      integer:: I, K
      real:: SIG, P, QN, UN

      DO I=1,N
        X(I)=real(LX(I))
      ENDDO

      Y2(1)=0.
      U(1)=0.
      
      DO I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
             /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO

      QN=0.
      UN=0.

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)

      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      ENDDO

      END SUBROUTINE SPLINE


      subroutine standev(nsize, indata, trend, std)
!-------------------------------------------------------------------
!  This is a program to obtain standard deviation of the linearly
!  detrended data
!  
!  PARAMETERS:
!  nsize    :    indata size
!  indata   :    input data
!  trend    :    the linear trend of "indata"
!  std      :    standard deviation
!--------------------------------------------------------------------

      implicit none

      integer, intent(in) :: nsize
      real, dimension(nsize) :: indata
      real, dimension(nsize), intent(out):: trend
      real, intent(out) :: std
      real sigmaX, sigmaY, sigmaX2, sigmaXY, real_nsize, Xbar, Ybar, real_i
      integer :: i
      real :: trend_const, trend_slope, temp

      sigmaX = 0.0
      sigmaY = 0.0
      do i = 1, nsize
        sigmaX = sigmaX + real(i)
        sigmaY = sigmaY + indata(i)
      enddo

      real_nsize=real(nsize)

      Xbar=sigmaX/real_nsize
      Ybar=sigmaY/real_nsize

      sigmaX2 = 0.0
      sigmaXY = 0.0
      do i = 1, nsize
        real_i=real(i)
        sigmaX2 = sigmaX2+(real_i-xbar)*(real_i-xbar)
        sigmaXY = sigmaXY+(real_i-xbar)*(indata(i)-ybar)
      enddo

      trend_slope=sigmaXY/sigmaX2
      trend_const=Ybar-trend_slope*Xbar

      do i=1, nsize
        trend(i)= trend_const + trend_slope*real(i)
      enddo

      std=0.0

      do i=1,nsize
        temp= indata(i)-trend(i)
        std=std+temp*temp
      enddo
      std=std/real(nsize)
      std=sqrt(std)

      end subroutine standev
!
end module spatiotemporal_filter
