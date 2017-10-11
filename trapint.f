
!----------------- BEGIN Trapezoidal_Integration -----------------------------
      subroutine trapint(y, npts, dt, yint_0, yint)


* Integrate a time series y with initial condition yint_0

 
* Dates: 02/07/11 - Written by D.M. Boore, patterned after acc2vd.for
      real y(*), yint(*) 
      double precision cumv, y1, y2, v1, v2,
     : ddt, ddt_2 

 
* Integrate

      ddt     = dble(dt)
      ddt_2   = 0.5d0 * dble(dt)
 
      cumv = 0.0d0
 
      v1 = 0.0d0
      v2 = 0.0d0
      
      yint(1) = yint_0
 
      do j=2,npts
      
        y1 = dble(y(j-1))
        y2 = dble(y(j))
        
        cumv = cumv + (y1 + y2)*ddt_2
        v1 = v2
        v2 = cumv + dble(yint_0)  ! the yint_0 term is the constant of integration
        yint(j) = sngl(v2)
        
      end do

      return
      end
!----------------- END Trapezoidal_Integration -----------------------------

