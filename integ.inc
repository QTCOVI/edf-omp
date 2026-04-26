       integer(kind=4), parameter :: mprop=80,minter=10
c
       integer(kind=4) :: nprop, nvar, npt, npp, npr
       integer(kind=4) :: intepro, nphi, ntheta, nphib
       integer(kind=4) :: nthetab, ikey, mask, inuc
       real(kind=8) :: cost, sint, sintcosp, sintsinp
       real(kind=8) :: xnuc
       real(kind=8) :: rmin, rmax, the, phi, distmax, epsode
       real(kind=8) :: rmaxsurf, prop, propx, propy, propf 
       real(kind=8) :: angx, angy, angz
c
       common /iinteg/  inuc
       common /finteg/  xnuc(3)
       common /ffinteg/ cost,sint,sintcosp,sintsinp
c
       common /gfinteg/ rmin,rmax,the,phi,distmax,epsode
       common /iprop/   nprop,nvar
       common /iiprop/  mask(mprop) 
       common /fprop/   prop(mprop),propx(mprop),propy(mprop),
     &                  propf(mprop)
c
       common /npp/     npt,npp,npr,intepro,nphi,ntheta, 
     &                  nphib,nthetab,ikey
       common /surlim/  rmaxsurf
       logical          dosurf,dofullp,dophi,dothe,convert,leb,lebb,
     &                  xsurface
       common /nl/      dosurf,dofullp,dophi,dothe,convert,leb,lebb,
     &                  xsurface
       common /rotangu/ angx,angy,angz
c
c$omp threadprivate(/ffinteg/)
