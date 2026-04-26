c
c----------------------------------------------------------------------
c
      subroutine gtogto ()
c
c.....Overlap matrix between primitive Cartesian Gaussian Functions
c
      USE        space_for_wfnbasis
      USE        space_for_primgto
      include   'implicit.inc'
      include   'param.inc'
      include   'wfn.inc'
      include   'primgto.inc'
      real(kind=8)     ax(1:3),bx(1:3)
      real(kind=8),    allocatable,dimension (:,:) :: ab2
      parameter (lamx=12)
      integer it(3),jt(3)
c
c.....ceabx() are the coefficients, except for the factor 
c     EXP(-XMU*R_AB^2), where XMU=a*b/(a+b), that result from the 
c     expansion of the product of two primitive cartesian Gaussian
c
      real(kind=8) ceabx(-1:2*lamx,-1:lamx,-1:lamx,3)
      integer t,u,v
c
      allocate (ab2(ncent,ncent))
      do ica=1,ncent
        do icb=1,ica
          abaux=0d0
          do j=1,3
            abaux=abaux+(xyz(ica,j)-xyz(icb,j))**2
          enddo
          ab2(ica,icb)=abaux
          ab2(icb,ica)=abaux
        enddo
      enddo
c
c.....Compute the electronic molecular electrostatic potential.
c
      do ica=1,ncent
        ax(1:3)=xyz(ica,1:3)
        do ma=1,ngroup(ica)   
          nua=nuexp(ica,ma,1)
          itipa=ityp(nua)
          la=nlm(itipa,1)+nlm(itipa,2)+nlm(itipa,3)
          za=oexp(nua)
          do icb=1,ica
            bx(1:3)=xyz(icb,1:3)
            do mb=1,ngroup(icb)
              nub=nuexp(icb,mb,1)
              itipb=ityp(nub)
              lb=nlm(itipb,1)+nlm(itipb,2)+nlm(itipb,3)
              zb=oexp(nub)
              p=za+zb
              xmu=za*zb/p
              prefactor=exp(-xmu*ab2(ica,icb))
              pioverp=pi/p
              pioverp=sqrt(pioverp*pioverp*pioverp)
              do j=1,3
                call etijcalc (j,lamx,la,lb,ceabx,za,zb,ax(j),bx(j))
              enddo
c
c.............Compute the target functions for all the products of
c             of Gaussian primitives.
c
              do ka=1,nzexp(ica,ma)
                nua=nuexp(ica,ma,ka)
                itipa=ityp(nua)
                i=nlm(itipa,1)
                k=nlm(itipa,2)
                m=nlm(itipa,3)
                do kb=1,nzexp(icb,mb)
                  nub=nuexp(icb,mb,kb)
                  if (nua.ge.nub) then
                    itipb=ityp(nub)
                    j=nlm(itipb,1)
                    l=nlm(itipb,2)
                    n=nlm(itipb,3)
                    prod=ceabx(0,i,j,1)*ceabx(0,k,l,2)*ceabx(0,m,n,3)
                    sprim(nua,nub)=prod*pioverp*prefactor
                    sprim(nub,nua)=sprim(nua,nub)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate (ab2)
 56   format (3(1x,F15.8))
c
c-----Set to .true. that sprim() is already computed.
c
      sprimok = .true.
      return
      end
