      program testcg
      integer l1,l2,l,m1,m2
      call fctrln(250)
      l1=0
      do while(l1.ge.0) 
         write(6,*) 'Type l1, l2, l, m1, m2:-'
         read(5,*) l1,l2,l,m1,m2
         write(6,*) cgcoef(l1,l2,l,m1,m2)
      end do
      stop
      end
