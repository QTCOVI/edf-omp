c
c-----------------------------------------------------------------------
c
      subroutine dosym (ncent,nmo,nprims,neq,uout,modelsym,largwr)
c
c.......................................................................
c
      USE      space_for_wfnbasis
      USE      mod_sym
      USE      mod_periodic
      USE      space_for_sym

      implicit none
      integer(kind=4) ncent,nmo,nprims,neq
      integer(kind=4) i,uout,modelsym,ier
      logical  largwr

      include     'datatm.inc'

      real(kind=8),    allocatable,dimension (:)   :: mass
      real(kind=8),    allocatable,dimension (:)   :: ax,ay,az
      integer(kind=4), allocatable,dimension (:)   :: atnum
c
c-----Symmetry analysis
c
      if (modelsym.eq.2) then
c
c-------Symmetry analysis with O. Beruski and L.N. Vidal code
c       lnvidal@utfpr.edu.br
c
        call allocate_space_for_sym (nmo,ncent)
        call newsym (ncent,uout,neq,xyz,charge,ineq,mult,atnam)

c
c-------Symmetry analysis with Jose Luis Casals Sainz code
c
      elseif (modelsym.eq.3) then
        allocate (ax(ncent),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate ax()'
        allocate (ay(ncent),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate ay()'
        allocate (az(ncent),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate zy()'
        allocate (atnum(ncent),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate atnum()'
        allocate (mass(ncent),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate mass()'
        do i=1,ncent
          ax(i)=xyz(i,1)
          ay(i)=xyz(i,2)
          az(i)=xyz(i,3)
          atnum(i)=nint(charge(i))
          mass(i)=wgatm(atnum(i))
        enddo
        call symheader (uout)
        call syminit (uout,.true.)
        call init_periodic_table()           
        call symrwcoord (uout,ax(:),ay(:),az(:),
     &       atnum(:),mass(:),.true.,.false.,ncent)  
        call symdeter (uout,ax(:),ay(:),az(:),ncent)
        if (mol_linear) then
          call sym1d (uout,ax(:),ay(:),az(:),atnum(:),ncent)
        else if (mol_planar) then
          call sym2d (uout,ax(:),ay(:),az(:),atnum(:),ncent)
        else
          call sym3d (uout,ax(:),ay(:),az(:),atnum(:),ncent)
        endif
        call symprint (uout,.false.)
      else
c
c-------Symmetry analysis with the standard code of promolden
c
        call allocate_space_for_sym (nmo,ncent)
        call sym (uout,largwr)
        call deallocate_space_for_sym ()
      endif
      return
      end
