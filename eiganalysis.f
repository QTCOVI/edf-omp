      subroutine eiganalysis (n,nm,nf,a1,a2,nl,it,lw,wr,dow,typ,ei,ov)
      implicit none
      integer(kind=4) n,nm,nf,nl,k,a1,a2,it,lw
      real(kind=8) ov,ei(1:n)
      logical wr,dow
      character(len=5) typ

      nf=nm
      do k=nf,1,-1
        if (ei(k).lt.ov) then
          if (wr) then
            if (dow) then
              write (lw,113) typ(1:5),ov,it
              dow = .false.
            endif
            write (lw,117) k,ei(k)
          else
            exit
          endif
        else
          nm = nm - 1
          nl = nl + 1
          if (dow) then
            write (lw,113) typ(1:5),ov,it
            dow = .false.
          endif
          if (a2.eq.0) write (lw,116) k,ei(k),a1
          if (a2.ne.0) write (lw,116) k,ei(k),a1,a2
        endif
      enddo
      return
 113  format (/,1x,70('-'),/,
     & 2x,a,'center localizations  Critical overlap = ',F12.5,
     & 2x,'ITER = ',I1,/,1x,70('-'))
 116  format (3x,'MO ',I4,' ( eig = ',F12.5,
     &         ' )     localized on atom(s) ',100I4)
 117  format (3x,'MO ',I4,' ( eig = ',F12.5,' ) non localized yet')
      end
