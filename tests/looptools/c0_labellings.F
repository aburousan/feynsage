      program check3
      implicit none
#include "looptools.h"
      double complex cv
      integer kk
      double precision av(6,6)
      data av / -4D0,-3D0,-1D0, 0D0,0D0,1D0,
     &          -1D0,-3D0,-4D0, 0D0,1D0,0D0,
     &          -3D0,-1D0,-4D0, 1D0,0D0,0D0,
     &          -4D0,-1D0,-3D0, 0D0,0D0,1D0,
     &          -1D0,-4D0,-3D0, 1D0,0D0,0D0,
     &          -3D0,-4D0,-1D0, 0D0,1D0,0D0 /
      call ltini
      call setmudim(1D0)
      call setversionkey(2*KeyC0)
      do kk = 1, 6
        call clearcache
        write(*,"(A,6F5.1)") "LABELLING", av(:,kk)
        cv = C0(av(1,kk),av(2,kk),av(3,kk),av(4,kk),av(5,kk),av(6,kk))
        write(*,"(A,2F22.16)") "RETURNED", cv
      enddo
      call ltexi
      end
