      recursive subroutine searchstring (arr,key,imin,imax,ifound,found)
      implicit none
      character(len=*) key, arr(*)
      integer imin,imax,imid,ifound,minimo,maximo
      common /maxmin/ minimo,maximo
      logical found
!
      if (imax < imin) then
        found = .false.
        minimo = imax
        maximo = imin
        return
      else
        imid = (imin+imax)/2
             if ( arr(imid) > key ) then
          call searchstring (arr,key,imin,imid-1,ifound,found)
        else if ( arr(imid) < key ) then
          call searchstring (arr,key,imid+1,imax,ifound,found)
        else
          found = .true.
          ifound = imid
          return
        endif
      endif
      end
