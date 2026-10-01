      program check4
      implicit none
#include "looptools.h"
      double complex dv
      call ltini
      call setmudim(1D0)
      call setversionkey(2*KeyD0)
      write(*,"(A)") "POINT 113 (sorted D0 args, masses squared)"
      dv = D0(0.25D0,1D0,-0.5D0,0D0,-0.25D0,0.25D0, 4D0, 0.5625D0, 7.5625D0, 4D0)
      write(*,"(A,2F24.18)") "RETURNED", dv
      call clearcache
      write(*,"(A)") "POINT 10"
      dv = D0(-0.5D0,-0.25D0,0.75D0,0D0,-1D0,0.25D0, 7.5625D0, 9D0, 7.5625D0, 0.5625D0)
      write(*,"(A,2F24.18)") "RETURNED", dv
      call ltexi
      end
