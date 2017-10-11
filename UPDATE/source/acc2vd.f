
!----------------- BEGIN Acc2VD -----------------------------
      subroutine acc2vd(acc, npts, dt, rmv_trnd, vel0, dis0, vel, dis)


* Compute velocity and displacement time series from acceleration,
* assuming that the acceleration
* is represented by straight lines connecting the digitized values.

* Dates: 02/09/99 - Written by D.M. Boore
*        01/07/00 - Added initial velocity and displacements (v0, d0).
*                   This also requires the addition of a linear trend in
*                   displacement.
*                   Also, Bill Joyner, Chris Stephens, and I considered how to
*                   handle the first point.  Assuming v0 = 0, BAP uses a 
*                   trapezoidal rule such that the first vel(1) = 0.5*dt*a(1),
*                   but in acc2vd, vel(1) = 0.  In effect, BAP starts the
*                   integration at -dt rather than 0.0, with the result that
*                   because vel(1) depends on dt, vel(1) will change if the
*                   digitization interval is changed.  We agreed that this is 
*                   not correct.
*        12/18/02 - Some minor changes, such as specifying 0.0d0 rather 
*                   than 0.0 and using the function "dble" in various places
*        10/15/07 - Rename "d0" to "dis0" so as not to confuse with "d0"
*                   in "0.0d0".  For consistency, change "v0" to "vel0".
*                   More importantly, restructure the calculations so as not to
*                   lose the double-precision estimates of velocity when
*                   accumulating the displacement, and abandon the calculation 
*                   of disp based on integrating straightlines in acceleration, 
*                   using a simple trapezoidal rule instead.  These changes were
*                   made in an attempt to track down a drift in displacement
*                   for the HHE component of the KYTH recording of the
*                   Kythera earthquake recorded on a station of the EGELADOS
*                   velocity-sensor network.   The drift did not appear when I
*                   integrated the acceleration in Excel, and the SD computed 
*                   using SMC2RS also showed no indication of the drift (the
*                   long-period SD equals the peak displacement during the
*                   strong part of the shaking).   Unfortunately, the changes
*                   made today did not eliminate the drift.   My next check is
*                   to see if the 5 digit resolution in the smc files might
*                   be responsible (the data given to me by Andreas Skarlatoudis
*                   were in column files with 7 digits of resolution, with a maximum
*                   exponent of 9).  When I 
*                   opened the asc file made by smc2asc in Excel and computed dis, 
*                   I found the same trend as in the smc2vd file. This tells me
*                   that the problem is in the asc2smc conversion (but why the
*                   SD computed using the smc files shows no sign of the drift
*                   is a mystery to me).  I eventually confirmed that the problem
*                   is the limited resolution of the standard format; I have added
*                   an option to the smc utility programs to use a higher-resolution
*                   format (converting all of the utility programs will take some
*                   time---I am only doing those that I currently need right now).
!        02/06/11 - Added comments to make it clear that now use trapezoidal integration.
!                 - Delete "sngl(cumv)" and "sngl(cumd)" for vel(1) and dis(1)
!        02/07/11 - Call a subroutine to do the integration
!        06/01/11 - Delete the include statement for the integration routine; it should be placed in the programs that
!                   call acc2vd instead.  The reason for this is because
!                   the subroutine acc2vd is also included in \smsim\td_subs, and a compile
!                   error can occur if the \forprogs folder with the integration routine
!                   is not available.  This means that the trapezoidal_integration.for must be available in smsim.
!                   I guarantee this by including it in make_td_subs_file.bat, which should be called
!                   in time domain cl*.bat files.

      real acc(*), vel(*), dis(*) 
      logical rmv_trnd
      double precision cumv, cumd, a1, a2, v1, v2,
     : ddt, ddt_2, ddtdt_6

* compute velocity and displacement

      call Trapezoidal_Integration(acc, npts, dt, vel0, vel)
      call Trapezoidal_Integration(vel, npts, dt, dis0, dis)

      return
      end
!----------------- END Acc2VD -----------------------------


