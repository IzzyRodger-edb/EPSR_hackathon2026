      recursive subroutine minimagec(x,y,z)
      include 'lattice.inc'
      real x,y,z,fa,fb,fc
c
c subroutine to apply the minimum image convention of a point in Cartesian
c coordinates (x,y,z) relative to an arbitrary set of lattice vectors 
c a,b,c. Note that if specified as non-cubic the supplied coordinates do not have to
c be normalised, since the routine normalises them in any case.

      if(cubic) then

         call minimage(x,y,z)
         return

      else
c
c fractional coordinates relative to (a,b,c)
c
         fc=z/clatz
         fb=(y-fc*claty)/blaty
         fa=(x-fb*blatx-fc*clatx)/alatx
c
c minimum image convention
c
         call minimage(fa,fb,fc)
c
c calculate corrected coordinates
c
         z=fc*clatz
         y=fb*blaty+fc*claty
         x=fa*alatx+fb*blatx+fc*clatx
   
      endif
      
      return
      end

      subroutine minimage(fa,fb,fc)

      real fa,fb,fc
c
c subroutine to apply the minimum image convention of a point in Cartesian
c coordinates assuming cubic boundary conditions
c
c
c minimum image convention
c
      fc=fc-nint(fc)
      fb=fb-nint(fb)
      fa=fa-nint(fa)
      return
      end


