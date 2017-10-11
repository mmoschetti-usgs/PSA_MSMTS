
! ----------------------------------------------------------------- begin rscalc_interp_acc
      subroutine rscalc_interp_acc(acc, na, omega, damp_in, dt_in,
     :                             rd, rv, aa)
        
!-----------------------------------------------------------------------
! This version does not return response time series.

! Dates: 03/04/10 - Program rsp obtained from B. Chiou.  cmpmax written by
!                   I. Idriss; ucmpmx by R. Youngs.  D. Boore changed the input 
!                   parameters to be equivalent to rdrvaa.for
!        03/05/10 - Substituted rdrvaa for cmpmax
!        03/06/10 - Renamed from rsp_rdrvaa_ucmpmx to rscalc_interp_acc.
!                 - Renamed subroutine ucmpmx to icmpmx ("i" for "interpolate
!                   acceleration) and modified icmpmx.
!        08/13/12 - Norm Abrahamson suggests doing a more exact interpolation 
!                   when period < n*dt (n=10 here).   He suggests interpolating
!                   by adding zeros in the frequency domain, figuring
!                   out at the beginning what will be the shortest period desired,
!                   interpolating the input time series accordingly, and then feeding 
!                   this into rdrvaa.    This requires a major restructuring of this 
!                   subroutine.  I will not do this yet; I just wanted to write
!                   down the suggested revision someplace.
!        10/09/12 - Following up on the last comment, the interpolation is in the driver
!                   smc_interpolate_time_series_using_fork.for, and it uses 
!                   interpolate_time_series_using_fork.for.  I have not incorporated
!                   it yet into programs such as smc2rs, smc2rs2, or blpadflt.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents
!                 - Get rvrdaa.for via an include statement rather than including it in this file
!                   in order to get the most recent version (I had done this before, but for
!                   some reason had copied rvrdaa.for to the bottom of this file, probably
!                   to simplify the distribution of the program).

      implicit none

      real(4) :: acc(*), omega, damp_in, dt_in, 
     :       rd, rv, aa, d0, v0
      
      integer :: na, kg, kug, npr


      real(8) :: ug(:), pr, damp, dt, z(3)

      real(8) :: w, twopi
      
      allocatable :: ug
      
      integer :: i, nn

      dt = dble(dt_in)
      damp = dble(damp_in)
      
      kg = na
      
      allocate(ug(na))

          do i=1,kg
            ug(i) = dble(acc(i))
          enddo
!...
!... Compute response spectra
!...
           kug=kg-1
           
      w = dble(omega)
      
      twopi = 4.0d0*dasin(1.0d0)
      
      pr = twopi/w

      if(dt == 0.0d0 .or. pr < 10.0d0*dt) then
        call icmpmx(kug, ug, dt, pr, w, damp, z)
        rd = sngl(z(1))
        rv = sngl(z(2))
        aa = sngl(z(3))
      else
        d0 = 0.0
        v0 = 0.0
        call rdrvaa(acc,kg,omega,damp_in,dt_in,rd,rv,aa, d0, v0)
      endif
            
      deallocate(ug)

      return
      end

!-----------------------------------------------------------------------
      subroutine icmpmx(kug, ug, dt_in, pr, w, d, z)
      
! z(1) = SD
! z(2) = RV
! z(3) = AA

! Dates: 03/06/10 - Original program ucmpmx (written by Bob Youngs/I. Idriss)
!                   renamed icmpmx and modified to assume equal time spacing
!                   if the original and interpolated acceleration time
!                   series (the "u" in the original name referred to 
!                   "u"nequal spacing).
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents

      implicit none

! Input
      integer :: kug
      real(8) :: ug(*), pr, w, d, dt_in

! Output
      real(8) :: z(*)

! Working variables
      real(8) :: t(3), c(3), x(2,3)
      real(8) :: f1, f2, f3, f4, f5, f6, wd, w2, w3
      real(8) :: dt, e, g1, g2, h1, h2, dug, g, z1, z2, z3, z4
      real(8) :: a, b
      integer :: nn, i, k, ns, is, j
!
      nn=1
      wd=sqrt(1.-d*d)*w
      w2=w*w
      w3=w2*w
      DO i=1,3
        x(1,i)=0.
        z(i)=0.
      END DO
      
      f2=1./w2
      f3=d*w
      f4=1./wd
      f5=f3*f4
      f6=2.*f3
      
      ns= int(10.*dt_in/pr-0.01)+1   !! 05/05/2008
      dt=dt_in/real(ns)
      
      DO k=1,kug
      
        f1=2.*d/w3/dt
        e=dexp(-f3*dt)
        g1=e*dsin(wd*dt)
        g2=e*dcos(wd*dt)
        h1=wd*g2-f3*g1
        h2=wd*g1+f3*g2
        dug=(ug(k+1)-ug(k))/real(ns)
        g=ug(k)
        z1=f2*dug
        z3=f1*dug
        z4=z1/dt
        
        DO is=1,ns
          z2=f2*g
          b=x(1,1)+z2-z3
          a=f4*x(1,2)+f5*b+f4*z4
          x(2,1)=a*g1+b*g2+z3-z2-z1
          x(2,2)=a*h1-b*h2-z4
          x(2,3)=-f6*x(2,2)-w2*x(2,1)
          nn = nn + 1
          DO j=1,3
            c(j)=abs(x(2,j))
            IF (c(j) > z(j)) THEN
              z(j)=c(j)
            END IF
            x(1,j)=x(2,j)
          END DO
          
          g=g+dug
          
        END DO
      END DO
  
      RETURN
      END

!------------------------------
 
!      include '../tspp_subroutines/rdrvaa.for'

! ----------------------------------------------------------------- end rscalc_interp_acc

! NOTE: This program was incorporated into rscalc_interp_acc.for.
! Note that rdrvaa includes initial conditions in the argument
! list, while rscalc_interp_acc does not.

!----------------- BEGIN RDRVAA -----------------------------
      subroutine rdrvaa(acc,na,omega,damp,dt,rd,rv,aa, d0, v0)
! This is a modified version of "Quake.For", originally
! written by J.M. Roesset in 1971 and modified by
! Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
! Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  

!   acc = acceleration time series
!    na = length of time series
! omega = 2*pi/per
!  damp = fractional damping (e.g., 0.05)
!    dt = time spacing of input
!    rd = relative displacement of oscillator
!    rv = relative velocity of oscillator
!    aa = absolute acceleration of oscillator
! d0,v0 = initial displacement and velocity (usually set to 0.0)

! Dates: 02/11/00 - Modified by David M. Boore, based on RD_CALC
!        03/11/01 - Double precision version
!        03/14/01 - Added d0, v0 (note on 05 March 2010: I recommend 
!                   that they not be used, by setting them to 0.0.  
!                   I've kept them as arguments in the subroutine call
!                   so that I do not have to modify programs that use
!                   this subroutine).                   
!        03/14/01 - Changed name back to rdrvaa
!        01/31/03 - Moved implicit statement before the type declarations
!        10/10/07 - Initial variable assignments and iteration loop modified 
!                   to double-precision (Chris Stephens)
!        03/05/10 - Delete old (single precision) lines of code
!        12/22/10 - Remove minus sign in front of the initialization of y, ydot. 
!                   The minus sign was a remnant of an earlier version where I did not
!                   understand the meaning of y and ydot.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents

      real(4) :: acc(*), omega, damp, dt, rd, rv, aa, d0, v0
      
      real(8) :: d2, bom, d3, omd, om2, c1, omt, omdt, c2, c3, c4, ss, 
     :           cc, bomt, ee, s1, s2, s3, a11, a12, a21, a22,
     :           s4, s5, b11, b12, b21, b22, y, ydot, y1, z, z1, z2, ra
     
      integer :: i, na
      
      
      d2=1.d0-dble(damp)*dble(damp)
      d2=dsqrt(d2)
      bom=dble(damp)*dble(omega)
      d3 = 2.d0*bom                 ! for aa
      omd=dble(omega)*d2
      om2=dble(omega)*dble(omega)
      c1=1.d0/om2
      
      omt=dble(omega)*dble(dt)
      omdt=omd*dble(dt)
      c2=2.d0*dble(damp)/(om2*omt)
      c3=c1+c2
      c4=1.d0/(dble(omega)*omt)
      ss=dsin(omdt)
      cc=dcos(omdt)
      bomt=dble(damp)*omt
      ee=dexp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.d0-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      
      rd=0.
      rv = 0.                           ! for rv
      aa = 0.                           ! for aa
      
      y=    dble(d0)
      ydot= dble(v0)    
!      y=0.
!      ydot=0.

      do i=1, na-1
      
        y1=a11*y+a12*ydot+b11*dble(acc(i))+b12*dble(acc(i+1))
        ydot=a21*y+a22*ydot+b21*dble(acc(i))+b22*dble(acc(i+1))
        y=y1    ! y is the oscillator output at time corresponding to index i
        z=dabs(y)
        if (z > rd) rd=z
        z1 = dabs(ydot)                   ! for rv
        if (z1 > rv) rv = z1            ! for rv
        ra = -d3*ydot -om2*y1            ! for aa
        z2 = dabs(ra)                     ! for aa
        if (z2 > aa) aa = z2            ! for aa
        
      end do
      
      return
      end
!----------------- END RDRVAA -----------------------------


