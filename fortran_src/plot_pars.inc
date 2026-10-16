      character*256 fname1,fname2,fname1save,fname2save
      character*256 fnameplot
      character*20 ext1,ext2,style,histogram,legend,grid
      character*29 plabel(100)
      character*256 title,xlabel,ylabel,zlabel
      character*256 logx,logy

      integer*4 plottype,nplot,ncolumn,nblock,npcolumn
     *,icolumn(100),ilenfname1,ilenfname2,ilenfnameplot,ilentitle
     *,ilenxlabel,ilenylabel,ilenzlabel,ilenlogx,ilenlogy
     *,ilen1save,ilen2save,nplotsave,ncolumnsave,npcolumnsave
     *,xcolumn,ilenhistogram,ilenstyle,ilenlegend,ilengrid
c
c store the symbol numbers for the plot, filename and block symbols
c
     *,icountp,icountf,icountb,ilenext1,ilenext2
c
c offsets and spacing for each block in a file
c 
     *,boffset,boffset2,bspace

      real*4 xmin,xmax,ymin,yup,ymax,ydel,yshift,yfact1,yfact2
     *,xfrac,yfrac

      common/plot_c/fname1,fname2,fname1save,fname2save
     *,fnameplot,ext1,ext2,style,plabel,title,xlabel,ylabel,zlabel
     *,logx,logy,histogram,legend,grid

      common/plot_i/plottype,plottypeold,nplot,ncolumn,nblock,npcolumn
     *,icolumn,ilenfname1,ilenfname2,ilenfnameplot,ilentitle,ilenxlabel
     *,ilenylabel,ilenzlabel,ilenlogx,ilenlogy,ilen1save,ilen2save
     *,nplotsave,ncolumnsave,npcolumnsave,xcolumn,ilenhistogram
     *,ilenstyle,ilenlegend,ilengrid,icountp,icountf,icountb,ilenext1,ilenext2,boffset
     *,boffset2,bspace

      common/plot_r/xmin,xmax,ymin,yup,ymax,ydel,yshift,yfact1,yfact2
     *,xfrac,yfrac
