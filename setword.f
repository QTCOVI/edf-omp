c
c-----------------------------------------------------------------------
c
      logical function setword (word, line, lp)
c
c.....setword - get next word from line at lp and increase lp
c             - a word is defined as any sequence on nonblanks
c             - setword returns .true. if it finds a non-empty word,
c               false otherwise
c
      character*(*)     word, line
      character*1       blank,newline,null
      integer           lp
      integer           l, ll
      logical           done
c
      blank = char(32)
      newline = char(10)
      null = char(0)
c
c.....test end of string
c
      l = len(word)
      ll = len(line)
      setword = .false.
      word = null
      if (lp.gt.ll) return
c
c.....skip blanks
c
      done = .false.
      do while (lp.le.ll .and. .not.done)
         done = line(lp:lp).ne.blank
         lp = lp + 1
      enddo
      if (.not.done) return
      lp = lp - 1
c
c.....get the word
c
      done = .false.
      i = lp
      do while (lp.le.ll .and. .not.done)
         done = line(lp:lp).eq.blank .or.
     &          line(lp:lp).eq.null  .or.
     &          line(lp:lp).eq.newline
         lp = lp + 1
      enddo
      if (done) lp = lp - 1
      word = line(i:lp-1)
      setword = .true.
      return
      end
