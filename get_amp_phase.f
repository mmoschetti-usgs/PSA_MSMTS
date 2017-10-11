* ----------------------------- BEGIN GET_AMP_PHASE --------------------------
      subroutine get_amp_phase(y, npts, amp, phase)

* Dates: 02/24/99 - written by D. Boore
*        05/21/02 - Used cmplx in filling work, use allocatable array for work

      real y(*), amp(*), phase(*)
      complex work(:)

      allocatable :: work
      allocate (work(npts+2))

      do i = 1, npts
        work(i) = cmplx(y(i), 0.0)
c        write(*,*)y(i)
      end do
      call fork(npts,work,1.0)
c      print *, "Random value", work(251)

      do i = 1, npts
        amp(i) = cabs(work(i))
        if (aimag(work(i)) .eq. 0.0 .and. 
     :       real(work(i)) .eq. 0.0)  then
          phase(i) = 0.0
        else
          phase(i) = atan2( aimag(work(i)), real(work(i)) )
        end if
c        write(*,*)amp(i)
      end do
      
      deallocate (work)
      
      return
      end
* ----------------------------- END GET_AMP_PHASE --------------------------

             

! ----------------------------- BEGIN FORK --------------------------
      SUBROUTINE FORK(LX,CX,SIGNI)
! FAST FOURIER                                  2/15/69
!                          LX
!    CX(K) = SQRT(1.0/LX)* SUM (CX(J)*EXP(2*PI*SIGNI*I*(J-1)*(K-1)/LX))
!                          J=1                        FOR K=1,2,...,LX
!
!  THE SCALING BETWEEN FFT AND EQUIVALENT CONTINUUM OUTPUTS
!  IS AS FOLLOWS.
!
!
!     GOING FROM TIME TO FREQUENCY:
!             F(W)=DT*SQRT(LX)*CX(K)
!
!                  WHERE W(K)=2.0*PI*(K-1)*DF

!                  and    DF = 1/(LX*DT)
!
!
!     GOING FROM FREQUENCY TO TIME, WHERE THE FREQUENCY
!     SPECTRUM IS GIVEN BY THE DIGITIZED CONTINUUM SPECTRUM:
!
!             F(T)=DF*SQRT(LX)*CX(K)
!
!                  WHERE T(K)=(K-1)*DT
!
!
!  THE RESULT OF THE SEQUENCE...TIME TO FREQUENCY,POSSIBLE MODIFICATIONS
!  OF THE SPECTRUM (FOR FILTERING,ETC.), BACK TO TIME...
!  REQUIRES NO SCALING.
!
!
!  THIS VERSION HAS A SLIGHT MODIFICATION TO SAVE SOME TIME...
!  IT TAKES THE FACTOR 3.1415926*SIGNI/L OUTSIDE A DO LOOP (D.BOORE 12/8
!  FOLLOWING A SUGGESTION BY HENRY SWANGER).
!

! Some brief notes on usage:

! "signi" is a real variable and should be called either with the value "+1.0"
! of "-1.0".  The particular value used depends on the conventions being used
! in the application (e.g., see Aki and Richards, 1980, Box 5.2, pp. 129--130).

! Time to frequency:
! In calling routine,
! 
!       do i = 1, lx
!         cx(i) = CMPLX(y(i), 0.0)
!       end do
      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./real(lx))

      do i=1,lx
      
        if (i <= j) then
          ctemp=cx(j)*sc
          cx(j)=cx(i)*sc
          cx(i)=ctemp
        end if
        
        m=lx/2        
        
        DO
          if (j <= m) EXIT
          j=j-m
          m=m/2
          if (m < 1) EXIT
        END DO
        j = j + m
        
      end do
     
      l=1
      DO WHILE (l < lx)
        istep=2*l
        temp= pi * signi/real(l)

        do m=1,l
          carg=(0.,1.)*temp*(m-1)
          cw=cexp(carg)

          do i=m,lx,istep
            ctemp=cw*cx(i+l)
            cx(i+l)=cx(i)-ctemp
            cx(i)=cx(i)+ctemp
          end do
        end do

        l=istep
      END DO

      return
      end
! ----------------------------- END FORK --------------------------
