c
c-----------------------------------------------------------------------
c
      character*(*) function uppcase (word)
c
c.....uppcase - transform alphabetical characters in word to uppercase.
c     word and uppcase return the same value.
c
      character*(*)     word
      integer           despl
c
      ilowa = ichar('a')
      ilowz = ichar('z')
      iuppA = ichar('A')
      iuppZ = ichar('Z')
      despl = iuppA - ilowa
c
      do 10 i = 1, leng(word)
         ich = ichar(word(i:i))
         if (ich.ge.ilowa .and. ich.le.ilowz) then
            word(i:i) = char(ich+despl)
         endif
 10   continue
      uppcase = word
      end
