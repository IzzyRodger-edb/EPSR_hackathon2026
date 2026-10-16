      parameter (mclass=120)

c Variables used to store parameters read from gnuatoms.txt and gnubonds.txt
c and for plotting atoms and bonds.

c Output filenames

      character*256 fnamesplot(0:mcomp)

c Atom and bond class labels

      character*3 classlab(mclass),classlab1(mpairs),classlab2(mpairs)

c Labels for atoms and bonds to be plotted

      character*3 cbstr1(mpairs),cbstr2(mpairs)

c Atom radii and RGB values

      real classrad(mclass),classrgb(3,mclass)

c Bond minimum and maximum distances and radii

     *,classbond1(mpairs),classbond2(mpairs),bondrad(mpairs)

c Atomic coordinates of atoms to be plotted

     *,xat(matom),yat(matom),zat(matom)

c Parameters to control the units of the axes, and GNUplot fractional position of the xy plane
c and scaler for size of plotting points

     *,unit,posfrac,ptscaler

c Overall scale factor on atom sizes and x and y rotation angles of plot

     *,scaleatoms,xangle,yangle

c Maximum and minimum distances and radii for bonds read in

     *,dstr1(mpairs),dstr2(mpairs),radstr(mpairs)

c RGB values of bonds and box lines, and linewidth of box lines

     *,bondrgb(3),boxrgb(3),boxlwidth
 
c Class colours for atoms

      integer nclassread,classtype(mclass),classcolour(mclass),classcoloura(mclass)

c Reference integer for plotting bonds

     *,nclassbond,atomref(matom)

c Number of components and component numbers to be plotted

     *,nclplt,nplt(mcomp)

c Flag to signal whether each atom will be plotted, also contains the associated class number
c of that atom. Will be zero if this atom is not to be plotted.

     *,iplot(matom)

c Length of output filenames

     *,ilenfnamesplot(0:mcomp)

c Atom types for bonds

     *,ibstr1(mpairs),ibstr2(mpairs)

c No. of classes to be plotted, list of associated class numbers as read from gnuatoms.txt,
c and number and list of atom types to be plotted

     *,nclass,iplotclass(mcomp),ntypplot,itypplot(mcomp)

c Atom at centre of plot

     *,nmolcentre

c Logical variable used in writexyz to signify whether to append the current coordinates to an existing file
c or create a new file with a different version number

      logical appendxyz

      common /plotato_c/fnamesplot,classlab,classlab1,classlab2
     *,cbstr1,cbstr2

      common /plotato_r/classrad,classrgb
     *,classbond1,classbond2,bondrad
     *,xat,yat,zat,unit,posfrac,ptscaler,scaleatoms,xangle,yangle
     *,dstr1,dstr2,radstr,bondrgb,boxrgb,boxlwidth

      common /plotato_i/nclassread,classtype,classcolour,classcoloura
     *,nclassbond,atomref,nclplt,nplt,iplot,ilenfnamesplot,ibstr1,ibstr2
     *,nclass,iplotclass,ntypplot,itypplot,nmolcentre

      common /plotato_l/appendxyz
