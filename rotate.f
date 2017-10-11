
! --------------------------------------------- begin rotate
      subroutine Rotate(z1,z2,znull,azm1,azm2,azmr,nz)

! Rotates z1, z2 into azmr (stored in z1 on return) and
! azmr+90.0 (stored in z2 on return).  No assumptions are made 
! about the relation between azm1 and azm2.  Znull is a number that
! represents no data or clipped data; if encountered, both z1 and z2
! are replaced by znull.

! NOTE:  The azimuths can be absolute (clockwise from north), or azmr
!        can be relative to the azimuth for component 1, in which
!        case azm1 is 0.0 and azm2 is relative to azm1 (+-90).

! NOTE NOTE NOTE:  z1, z2, azm1, azm2 are both input and output variables
!  I've been fooled by this, calling rotate several times without resetting 
!  azm1 and azm2

! Conventions: all angles are clockwise from north and are measured
! in degrees.

! See NB#10, p. 22 for development of equations.

! Written by:  Dave Boore
! Dates: 07/07/88 - written
!        07/04/05 - replace isiz in real declaration with "*",
!                   and remove isize from calling list.
!        09/22/12 - Small changes in coding (e.g., remove statement number,
!                   replace .ge. with >=, etc).

      real*4 z1(*), z2(*)

      pi = 4.0*atan(1.0)
      dtor = pi/180.0

      cr1 = cos( (azmr-azm1)*dtor )
      sr1 = sin( (azmr-azm1)*dtor )
      cr2 = cos( (azmr-azm2)*dtor )
      sr2 = sin( (azmr-azm2)*dtor )

      loop over points: DO i = 1, nz

         if (z1(i) == znull .or. z2(i) == znull) then

            z1(i) = znull
            z2(i) = znull
            cycle loop over points

         endif

         cp1 = z1(i)
         cp2 = z2(i)

         z1(i) =  cr1*cp1 + cr2*cp2
         z2(i) = -sr1*cp1 - sr2*cp2
         
      END DO loop over points

      azm1 = azmr
      azm2 = azmr + 90.0
      if (azm2 >= 360.0) azm2 = azm2-360.0

      return
      end
! --------------------------------------------- end rotate 
