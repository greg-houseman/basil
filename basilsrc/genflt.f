C*--------------------------------------------------------------------
C*    Basil / Sybil:   genflt.f  1.1  20 August 2003
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE BNDFLTLIMITS(MIN,MAX)
      INCLUDE "limits.parameters"

      MIN = FLTMIN
      MAX = FLTMAX
      END

      SUBROUTINE BNDEXTLIMITS(MIN,MAX)
      INCLUDE "limits.parameters"

      MIN = EXTMIN
      MAX = EXTMAX
      END
