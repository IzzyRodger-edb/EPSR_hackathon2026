!****************************************************************************
module gif_util
!  Conversion of raster data to GIF format.
!
!  Version 1.01, August 1999
!  Written by Jos Bergervoet
! 2008 Jan 28: Modified by Clive Page to use stream I/O, array as colourmap.
!
implicit none         !  Check all declarations
private               !  bin_io is used private, no transfer to main program
public  :: writegif   !  Writes GIF89 image, given pixel array and color map
private :: giflzw, slicewrite, InitTable, flushbuffer
integer, parameter, private     :: Bufend=260
character(len=Bufend), private  :: buf
integer, private  ::  ibuf                            ! output buffer vars
integer, parameter, private  ::    maxcode = 4095 
integer, parameter, private  ::    nocode = maxcode+1 ! definitions for LZW

! Define LZW code tables for hashing:
character(len=1), private, dimension(0:maxcode+1)  :: endbyte
integer, private, dimension(0:maxcode)             :: follow, next
  !
  ! For any code P, which codes for a sequence af pixel-values, endbyte(P)
  ! is the last pixel-value, follow(P) points to another code (if it exists)
  ! which codes for this same sequence, but with one more pixel-value
  ! appended.
  !   For each code P, next(P) points to another code which codes for a
  ! similar sequence with only the endbyte different. This is a hashing
  ! pointer, for fast look-up.
  !   All pointers are 'nocode' if they point to nothing
  !
  
integer, private :: ncod, curmaxcode, EOI, CC, P, K, child, &
                    maxbase, skip, slen, blen, accum, nout   ! local vars

contains
!-----------------------------------------------------------------------------
! CHAR2 Converts the two least sig bytes of an integer to a 2-character string
character(len=2) function char2(ival)
integer, intent(in) :: ival
char2 = achar(mod(ival,256)) // achar(mod(ival/256,256))
end function char2
!-----------------------------------------------------------------------------
subroutine flushbuffer(F_unit)
! Flushes up to 255 bytes to output file if buffer contains data, keeping
! rest of data in buffer. If skip>0 there is a partially filled last byte
! in buf[ibuf]. This byte will be written only if ibuf<256. That should be
! the last call to flushbuffer.
integer, intent(in) :: F_unit   ! I/O unit to use
  integer  :: bl    !   number of bytes to write (to be determined)

  if (ibuf > 255) then
    bl = 255        !   we will write buf[1..255]
  else if (skip /= 0) then
    bl = ibuf       !   buf[ibuf] is partially used, write buf[1..ibuf]
  else if (ibuf > 1) then
    bl = ibuf-1     !   write buf[1..ibuf-1], there is no partial byte
  else
    return          !   nothing to write
  end if

  write(F_unit) CHAR(bl)
  write(F_unit) buf(1:bl)
  buf(1:ibuf-bl) = buf(bl+1:ibuf) ! shift down remaining data
  ibuf = ibuf - bl
  return
end subroutine flushbuffer
!-----------------------------------------------------------------------------
subroutine giflzw(F_unit, Pixel)          ! routine for LZW coding
  integer, intent(in)                 :: F_unit
  integer, intent(in), dimension(:,:) :: Pixel
  integer                             :: i, j

  nout=0                          ! for counting the codes going out
  if (blen<2) then
    blen=2                        ! pixel code-length, 2 is minimum for GIF
  end if
  write(F_unit) CHAR(blen)
  maxbase = 2**blen - 1
  call InitTable()
  call slicewrite(F_unit, CC)

  do j=1, ubound(Pixel,2)
   do i=1, ubound(Pixel,1)
    K = modulo(Pixel(i,j), maxbase+1)    ! take next byte, prevent overflow
    if (i==1 .and. j==1) then
      P = K                       ! first raster byte has one-byte code P
      cycle                       ! for the first byte no further action
    end if
                                  ! Now see if code exists for sequence [.P.]K
    child = follow(P)             ! [.P.]K is "string coded by P" followed by K
    childloop: do
      if ((child == nocode) .or. (ichar(endbyte(child)) == K)) then
        exit childloop
      end if
      child = next(child)
    end do childloop

    if (child /= nocode) then     ! If code for [.P.]K was found, store it in P
      P = child
    else                          ! If not: output P and create code for [.P.]K
      call slicewrite(F_unit, P)
      if (ncod > maxcode) then    ! check if a new code can be added
        call slicewrite(F_unit, CC)       ! If not: tell listener to clear table
        call InitTable()          ! and clear our own table
      else
        if (ncod > curmaxcode) then
          slen = slen+1                     ! New codes will be one bit longer
          curmaxcode = curmaxcode * 2 + 1   ! and more codes are possible
        end if
        endbyte(ncod) = char(K)   ! ncod is the new code for [.P.]K
        follow(ncod) = nocode
        next(ncod) = follow(P)    ! include ncod in the hashing list
        follow(P) = ncod          !     of codes with same start-sequence
        ncod = ncod+1
      end if
      P = K
    end if
   end do
  end do
  call slicewrite(F_unit, P)              ! send the last code to buffer
  call slicewrite(F_unit, EOI)            ! send 'end of image' to buffer
  call flushbuffer(F_unit)        ! extra flush, including partial last byte
  return
end subroutine giflzw
!-----------------------------------------------------------------------------
subroutine InitTable()
  integer :: i

  do i=0,maxbase                  ! Start with defining the codes 0..maxbase
    endbyte(i) = char(i)          ! for one-pixel sequences (code=pixelvalue)
  end do                          ! Initially no multi-pixel codes exist
  follow(0:maxbase) = nocode
  next(0:maxbase) = nocode
  CC = maxbase+1                  ! `clear code-tabel', a control code
  EOI = maxbase+2                 ! `end of image', another control code
  ncod = CC + 2                   ! ncod = number of currently defined codes
  slen = blen + 1                 ! current number of bits to write one code
  curmaxcode = 2**slen - 1        ! currently the highest, until slen increases
  return
end subroutine InitTable
!-----------------------------------------------------------------------------
subroutine open_for_write(Fname, Funit)
! Creates a new Stream I/O file returning I/O unit used   
! CGP 2009 Jan 28
character(len=*), intent(in)  :: Fname
integer, intent(out)          :: Funit
!
logical :: exists, open
! Get free I/O unit number
do Funit = 90, 7, -1
    inquire(unit=Funit, exist=exists, opened=open)
    if(exists .and. .not. open) EXIT
end do
if(Funit < 7) STOP 'open_for_write failed - no free I/O units'
open (unit=Funit, file=Fname, access="STREAM", status="REPLACE")
end subroutine open_for_write
!-----------------------------------------------------------------------------
subroutine slicewrite(F_unit, code)       ! add some bits (a 'slice') to output buffer
  integer, intent(in)  :: F_unit
  integer, intent(in)  :: code

  if (nout == 0) then             ! initiate output buffer
    ibuf = 1
    skip = 0
    accum = 0
  end if
  nout = nout+1

  accum = accum + code * 2**skip  ! add bits at correct position in accum
  skip = skip + slen              ! slen is current slice length, in bits

  shiftout: do
    buf(ibuf:ibuf) = char(modulo(accum, 256))
    if (skip<8) then
      exit shiftout
    end if
    ibuf = ibuf+1                 ! last written buffer-byte is now permanent
    accum = accum / 256           ! remove that byte from accum
    skip = skip-8                 ! skip points to next bit to write in accum
  end do shiftout

  if (ibuf>255) then
    call flushbuffer(F_unit)            ! won't write unfinished byte in buf[ibuf]
  end if
  return                          ! at most 255 bytes will be left in buffer
end subroutine slicewrite
!-----------------------------------------------------------------------------
subroutine writegif (FileName, Pixel, ColorMap, Transparent)
!
! Codes pixel-map with palette into GIF format. Optional transparent color
!
character(len=*), intent(in)            :: FileName ! file to create or replace
integer, intent(in), dimension(:,:)     :: Pixel    ! Pixel values 0 to ncol
integer, intent(in), dimension(:,0:)    :: ColorMap ! RGB 0:255 for colours 0:ncol
integer, intent(in), optional           :: Transparent ! Optional
  
  character(len=256) :: s
  integer            :: InfoByte, nx, ny, Cblen, HasMap, maxincol,  &
                        maxgifcol, Background, i, F_unit
  
  call open_for_write (FileName, F_unit)
  nx = ubound(Pixel, 1)
  ny = ubound(Pixel, 2)
  maxincol = ubound(ColorMap,2) 
!!  print *,'image size', nx, ny, ' colours', maxincol
  do i=1,8                           ! find the bitsize, blen, for pixels
    blen = i
    maxgifcol = 2**blen - 1          ! Number of colors has to be power of 2 
    if (maxgifcol >= maxincol) then
      exit                           ! now blen and maxgifcol are correct
    end if                           ! only op to 256 colors can be 
  end do
   write(F_unit) "GIF89a" 
!  Create information for screen descriptor
  Background = 0
  if (present(Transparent)) then
    Background = Transparent
  end if
  HasMap = 1
  Cblen = blen
  InfoByte = HasMap * 128 + (Cblen-1) * 16 + blen-1
!  Write the screen descriptor
  write(F_unit) char2(nx), char2(ny), CHAR(InfoByte), CHAR(Background), CHAR(0)
  do i=0,maxgifcol                                 ! write global colormap
    write(F_unit) CHAR(colormap(1,min(i,maxincol))), &
                  CHAR(colormap(2,min(i,maxincol))), &
                  CHAR(colormap(3,min(i,maxincol)))
  end do
  if (present(Transparent)) then
    write(unit=*,fmt=*) "Transparent color: ", Transparent
    s = "!" // char(249) // char(4) // char(1) // char(0) // char(0)  &
            // char(Transparent) // char(0)
    write(F_unit) s(1:8)                            ! GIF transparent extension
  end if
   write(F_unit) ","                                ! Announce image
!    Now create and write image descriptor
  HasMap = 0
  InfoByte = HasMap * 128 + blen-1                 ! add 64, if interlaced
!    x_margin, y_margin (not used), image dimensions
  write(F_unit) char2(0), char2(0), char2(nx), char2(ny), CHAR(InfoByte) 
  call giflzw (F_unit, Pixel)                       ! now the raster data
  write(F_unit) CHAR(0), ';'                    ! Terminating 0-block ; for GIF
  close(unit=F_unit)
  return
end subroutine writegif
!-----------------------------------------------------------------------------
end module gif_util
!*****************************************************************************
!*GIDRIV -- PGPLOT GIF drivers
! Code converted using TO_F90 by Alan Miller  Date: 2009-01-26  Time: 23:50:47
! Substantially updated gidriv.f using above module.  Clive Page 2009 Jan 28
!+
MODULE gidriv_mod
USE gif_util
IMPLICIT NONE
PRIVATE
PUBLIC :: gidriv
CONTAINS

SUBROUTINE gidriv (ifunc, rbuf, nbuf, chr, lchr, mode)
INTEGER, INTENT(IN OUT)                  :: ifunc
REAL, INTENT(OUT)                        :: rbuf(*)
INTEGER, INTENT(OUT)                     :: nbuf
CHARACTER (LEN=*), INTENT(OUT)           :: chr
INTEGER, INTENT(OUT)                     :: lchr
INTEGER, INTENT(IN)                  :: mode
! PGPLOT driver for Graphics Interchange Format (GIF) files.
!***********************************************************************
!                           CAUTION                                    *
!                                                                      *
! The GIF specification incorporates the Lempel-Zev-Welch (LZW)        *
! compression technology which is the subject of a patent awarded to   *
! Unisys. Use of this technology, and in particular creation of GIF    *
! format files using this PGPLOT device driver, may require a license  *
! from Unisys.                                                         *
!***********************************************************************
! Supported device: GIF87a file format
! Device type codes: /GIF or /VGIF
! Default device name: pgplot.gif.
! If you have more than one image to plot (i.e. use PGPAGE) with this
! device, subsequent pages will be named: pgplot2.gif, pgplot3.gif,
! etc, disrespective of the device name you specified.
! You can however bypass this by specifying a device name including a
! number sign (#), which will henceforth be replaced by the pagenumber.
! Example: page#.gif will produce files page1.gif, page2.gif, ...,
! page234.gif, etc.
! Default view surface dimensions are:
! - GIF  : 850 x 680 pixels (translates to 10.0 x  8.0 inch).
! - VGIF : 680 x 850 pixels (translates to  8.0 x 10.0 inch).
! with an assumed scale of 85 pixels/inch.
! Default width and height can be overridden by specifying environment
! variables
! PGPLOT_GIF_WIDTH  (default 850)
! PGPLOT_GIF_HEIGHT (default 680)

! Color capability:
! Indices 0 to 255 are supported. Each of these indices can be assigned
! one color. Default colors for indices 0 to 15 are implemented.

! Obtaining hardcopy: Use a GIF viewer or converter.
!=
!  1-Aug-1994 - Created by Remko Scharroo
!  9-Aug-1994 - New scheme for line plotting
! 16-Aug-1994 - Provide multi-image plotting.
!  8-Sep-1994 - Add opcode 29 [TJP].
!  5-Nov-1994 - Adjust size of bitmap if necessary [TJP].
! 18-Jan-1995 - Attempt to prevent integer overflow on systems where
!               BYTE is signed [TJP].
! 28-Dec-1995 - prevent concurrent access [TJP].
! 29-Apr-1996 - use GRCTOI to decode environment variables [TJP].
!  7-Jul-2000 - Changed BYTE to INTEGER*2 [cgp]
! 28-Jan-2009 - Substantial re-write using Jos Bergervoet's GIF utilities
!-----------------------------------------------------------------------
INTEGER ::   bx, by
CHARACTER (LEN=*), PARAMETER :: ltype='GIF   (Graphics Interchange Format file, landscape orientation)'
CHARACTER (LEN=*), PARAMETER :: ptype='VGIF  (Graphics Interchange Format file, portrait orientation)'
CHARACTER (LEN=*), PARAMETER :: defnam='pgplot.gif'
INTEGER, PARAMETER :: dwd=850
INTEGER, PARAMETER :: dht=680
REAL, PARAMETER :: xres=85.
REAL, PARAMETER :: yres=xres
INTEGER :: UNIT, ic, npict, maxidx, state
INTEGER :: ctable(3,0:255), cdeflt(3,0:15)
INTEGER :: ier, i, l, ll, ix0, iy0, ix1, iy1, userw, userh
INTEGER :: grgmem, grfmem, grofil, grctoi
CHARACTER (LEN=80) :: msg, instr, filenm, filename
! Note: for 64-bit operating systems, change the following
! declaration to INTEGER*8:
!*
!* Changes for NAG f95 compiler: make PIXMAP allocatable
!*
INTEGER, allocatable, SAVE :: pixmap(:,:)

SAVE UNIT, ic, ctable, npict, maxidx, bx, by, filenm, filename
SAVE cdeflt, state
! Interchanged colour definitions below so default is white background, black lines
!  DATA cdeflt /000,000,000, 255,255,255, 255,000,000, 000,255,000,  &
DATA cdeflt /255,255,255, 000,000,000, 255,000,000, 000,255,000,  &
    000,000,255, 000,255,255, 255,000,255, 255,255,000,  &
    255,128,000, 128,255,000, 000,255,128, 000,128,255,  &
    128,000,255, 255,000,128, 085,085,085, 170,170,170/
DATA state /0/
!
SELECT CASE(ifunc)
CASE(1)
!--- IFUNC = 1, Return device name -------------------------------------
    IF (mode == 1) THEN
        chr = ltype
        lchr = LEN(ltype)
    ELSE IF (mode == 2) THEN
        chr = ptype
        lchr = LEN(ptype)
    ELSE
        CALL grwarn('Requested MODE not implemented in GIF driver')
    END IF
!--- IFUNC = 2, Return physical min and max for plot device, and range
!               of color indices ---------------------------------------
!     (Maximum size is set by GIF format to 2**16 pixels)
CASE(2)
    rbuf(1) = 0
    rbuf(2) = 65535
    rbuf(3) = 0
    rbuf(4) = 65535
    rbuf(5) = 0
    rbuf(6) = 255
    nbuf = 6
!--- IFUNC = 3, Return device resolution -------------------------------
CASE(3)
    rbuf(1) = xres
    rbuf(2) = yres
    rbuf(3) = 1
    nbuf = 3
!--- IFUNC = 4, Return misc device info --------------------------------
!    (This device is Hardcopy, supports rectangle fill, pixel
!     primitives, and query color rep.)
CASE(4)
    chr = 'HNNNNNPNYN'  ! CGP: changed to disable rectangle fill
    lchr = 10
!--- IFUNC = 5, Return default file name -------------------------------
CASE(5)
    chr = defnam
    lchr = LEN(defnam)
!--- IFUNC = 6, Return default physical size of plot -------------------
CASE(6)
    rbuf(1) = 0
    rbuf(2) = bx-1
    rbuf(3) = 0
    rbuf(4) = by-1
    nbuf = 4
!--- IFUNC = 7, Return misc defaults -----------------------------------
CASE(7)
    rbuf(1) = 1
    nbuf=1
!--- IFUNC = 8, Select plot --------------------------------------------
CASE(8)
    CONTINUE
!--- IFUNC = 9, Open workstation ---------------------------------------
CASE(9)
!     -- check for concurrent access
    IF (state == 1) THEN
        CALL grwarn('a PGPLOT GIF file is already open')
        rbuf(1) = 0
        rbuf(2) = 0
        RETURN
    END IF
!     -- dimensions of plot buffer
    userw = 0
    userh = 0
    CALL grgenv('GIF_WIDTH', instr, l)
    ll = 1
    IF (l > 0) userw = grctoi(instr(:l),ll)
    CALL grgenv('GIF_HEIGHT', instr, l)
    ll = 1
    IF (l > 0) userh = grctoi(instr(:l),ll)
    IF (mode == 1) THEN
!-- Landscape
        bx = dwd
        IF (userw >= 8) bx = userw
        by = dht
        IF (userh >= 8) by = userh
    ELSE
!-- Portrait
        bx = dht
        IF (userh >= 8) bx = userh
        by = dwd
        IF (userw >= 8) by = userw
    END IF
    npict=1
    maxidx=0
!-- Initialize color table
    ctable(1:3,0:15)   = cdeflt(1:3,0:15)
    ctable(1:3,16:255) = 128
    filenm = chr(:lchr)
    CALL grgi10 (filenm, npict, filename)
!    UNIT = grofil (msg)   ! Now a module procedure using Stream I/O
    rbuf(1) = UNIT
    rbuf(2) = 1
    state = 1
!--- IFUNC=10, Close workstation ---------------------------------------
CASE(10)
    state = 0
!--- IFUNC=11, Begin picture -------------------------------------------
CASE(11)
    bx = nint(rbuf(1))+1
    by = nint(rbuf(2))+1
    allocate(pixmap(bx,by))
!-- initialize to zero (background color)
    pixmap = 0
    IF (npict > 1) THEN
        CALL grgi10 (filenm, npict, filename)
!        UNIT = grofil(msg)
    END IF
!--- IFUNC=12, Draw line -----------------------------------------------
CASE(12)
    ix0=nint(rbuf(1))+1
    ix1=nint(rbuf(3))+1
    iy0=nint(rbuf(2))+1
    iy1=nint(rbuf(4))+1
    CALL grgi01(ix0, iy0, ix1, iy1, ic, bx, by, pixmap)
!--- IFUNC=13, Draw dot ------------------------------------------------
CASE(13)
    ix0=nint(rbuf(1))+1
    iy0=nint(rbuf(2))+1
    CALL grgi01(ix0, iy0, ix0, iy0, ic, bx, by, pixmap)
!--- IFUNC=14, End picture ---------------------------------------------
CASE(14)
    if(state == 1) then
       write(9,*) 'writegif: ', trim(filename)
       pixmap = pixmap(:,by:1:-1)
       call writegif(filename, pixmap, ctable)
! following line ok with g95, failed with gfortran
!       call writegif(filename, pixmap(:,by:1:-1), ctable)
    end if
!    state = 0
    npict = npict+1
    deallocate(pixmap)
!--- IFUNC=15, Select color index --------------------------------------
CASE(15)
    ic = rbuf(1)
    maxidx = MAX(maxidx, ic)
!--- IFUNC=16, Flush buffer. -------------------------------------------
!    (Not used.)
CASE(16)
    CONTINUE
!--- IFUNC=17, Read cursor. --------------------------------------------
!    (Not implemented: should not be called)
!--- IFUNC=18, Erase alpha screen. -------------------------------------
!    (Not implemented: no alpha screen)
CASE(18)
!--- IFUNC=19, Set line style. -----------------------------------------
!    (Not implemented: should not be called)
!--- IFUNC=20, Polygon fill. -------------------------------------------
!    (Not implemented: should not be called)
!--- IFUNC=21, Set color representation. -------------------------------
CASE(21)
    i = rbuf(1)
    ctable(1, i) = nint(rbuf(2)*255)
    ctable(2, i) = nint(rbuf(3)*255)
    ctable(3, i) = nint(rbuf(4)*255)
!--- IFUNC=22, Set line width. -----------------------------------------
!    (Not implemented: should not be called)
!--- IFUNC=23, Escape --------------------------------------------------
!    (Not implemented: ignored)
CASE(23)
    CONTINUE
!--- IFUNC=24, Rectangle fill (not implemented) -------------------------
CASE(24)
!--- IFUNC=25, Not implemented -----------------------------------------
CASE(25)
!--- IFUNC=26, Line of pixels ------------------------------------------
CASE(26)
    CALL grgi04(nbuf, rbuf, bx, by, pixmap, maxidx)
!--- IFUNC=27, Not implemented -----------------------------------------
CASE(27)
    CONTINUE
!--- IFUNC=28, Not implemented -----------------------------------------
CASE(28)
    CONTINUE
!--- IFUNC=29, Query color representation. -----------------------------
CASE(29)
    i = rbuf(1)
    rbuf(2) = ctable(1,i)/255.0
    rbuf(3) = ctable(2,i)/255.0
    rbuf(4) = ctable(3,i)/255.0
    nbuf = 4
CASE DEFAULT
    WRITE (msg,'(I10)') ifunc
    CALL grwarn('Unimplemented function in GIF device driver:' //msg)
    nbuf = -1
END SELECT
RETURN
END SUBROUTINE gidriv
!========================================================================
!*GRGI01 -- PGPLOT GIF driver, draw line
!+
SUBROUTINE grgi01 (ix0, iy0, ix1, iy1, icol, bx, by, pixmap)
INTEGER, INTENT(IN) :: ix0, iy0, ix1, iy1
INTEGER, INTENT(IN) :: icol
INTEGER, INTENT(IN) :: bx, by
INTEGER, INTENT(INOUT) :: pixmap(bx,by)
! Draw a straight-line segment from absolute pixel coordinates
! (IX0, IY0) to (IX1, IY1).
! Arguments:
!  ICOL            (input): Color index
!  PIXMAP   (input/output): The image data buffer.
!
INTEGER :: ix, iy, is
REAL :: d
! CGP: more error-trapping below to check (ix,iy) always in valid range
!! write(9,*) 'draw line', ix0, iy0, ix1, iy1, icol
IF (ix0 == ix1 .AND. iy0 == iy1) THEN
    IF(ix0 > 0 .and. ix0 <= bx .AND. iy0 > 0 .and. iy0 <= by) &
                        pixmap(ix0,iy0)=icol
ELSE IF (ABS(iy1-iy0) > ABS(ix1-ix0)) THEN
    d= (ix1-ix0) / REAL(iy1-iy0)
    is=1
    IF (iy1 < iy0) is=-1
    DO  iy=iy0,iy1,is
        ix=nint(ix0+(iy-iy0)*d)
        IF(ix > 0 .and. ix <= bx .and. iy > 0 .and. iy <= by) pixmap(ix,iy)=icol
    END DO
ELSE
    d=(iy1-iy0)/REAL(ix1-ix0)
    is=1
    IF (ix1 < ix0) is=-1
    DO  ix=ix0,ix1,is
        iy=nint(iy0+(ix-ix0)*d)
        IF(ix > 0 .and. ix <= bx .and. iy > 0 .and. iy <= by) pixmap(ix,iy)=icol
    END DO
END IF
END SUBROUTINE grgi01
!========================================================================
!*GRGI04 -- PGPLOT GIF driver, fill image line
SUBROUTINE grgi04(nbuf,rbuf,bx,by,pixmap,maxidx)
INTEGER, INTENT(IN)      :: nbuf
REAL, INTENT(IN)         :: rbuf(nbuf)
INTEGER, INTENT(IN)  :: bx, by
INTEGER, INTENT(INOUT) :: pixmap(bx,by)
INTEGER, INTENT(INOUT)   :: maxidx
INTEGER :: i,j, n,ic
!! write(9,*) 'fill image line', nbuf, rbuf(1:2), ' maxidx=', maxidx
i = nint(rbuf(1))+1
j = nint(rbuf(2))+1
DO  n=3,nbuf
    ic=rbuf(n)
    maxidx=MAX(maxidx,ic)
    pixmap(i+n-3,j)=ic
END DO
END SUBROUTINE grgi04
!========================================================================
!*GRGI10 -- Replace # in filename by picture number
SUBROUTINE grgi10 (name1, np, name2)
CHARACTER (LEN=*), INTENT(IN)  :: name1
INTEGER, INTENT(IN)            :: np
CHARACTER (LEN=*), INTENT(OUT) :: name2
! CGP: much updated to preserve .gif extension with picture number before it
CHARACTER (LEN=10) :: tmp
INTEGER :: idx
!
write(tmp, '(i0)') np
idx = INDEX(name1,'#')
IF (idx > 0) THEN
! name contains a #-character, replace it with the page number
    name2 = name1(1:idx-1) // trim(tmp) // name1(idx+1:)
ELSE IF (np == 1) THEN
! if this is the first page, use the supplied name
    name2 = name1
ELSE
! append an underscore and the page number to name before the dot
    idx = index(name1, '.', .true.)
    name2 = name1(1:idx-1) // "_" // trim(tmp) // name1(idx:)
END IF
CALL grwarn ('Writing new GIF image as: '//trim(name2))
END SUBROUTINE grgi10

END MODULE gidriv_mod
!****************************************************************************
!*GREXEC -- PGPLOT device handler dispatch routine
! Code converted using TO_F90 by Alan Miller
! Date: 2009-01-28  Time: 14:34:28
!+  Author: Clive Page 2009 Jan 28
SUBROUTINE grexec(idev,ifunc,rbuf,nbuf,chr,lchr)
use gidriv_mod
INTEGER, INTENT(IN OUT)                  :: idev
INTEGER, INTENT(IN OUT)                  :: ifunc
REAL, INTENT(OUT)                        :: rbuf(*)
INTEGER, INTENT(OUT)                     :: nbuf
CHARACTER (LEN=*), INTENT(IN OUT)        :: chr
INTEGER, INTENT(IN OUT)                  :: lchr
!---
INTEGER, PARAMETER :: ndev=13
CHARACTER (LEN=10) :: msg
!
SELECT CASE ( idev )
CASE(0)
  rbuf(1) = ndev
  nbuf = 1
  CASE ( 1)
    CALL nudriv(ifunc,rbuf,nbuf,chr,lchr)
  CASE ( 2)
    CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,1)
  CASE ( 3)
    CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,2)
  CASE ( 4)
    CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,3)
  CASE ( 5)
    CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,4)
  CASE ( 6)
    CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,1)
  CASE ( 7)
    CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,2)
  CASE ( 8)
   CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,3)
  CASE ( 9)
   CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,4)
  CASE (10)
    CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,5)
  CASE (11)
    CALL gwdriv(ifunc,rbuf,nbuf,chr,lchr,6)
  CASE (12)
    CALL gidriv(ifunc, rbuf, nbuf, chr, lchr, 1)
  CASE (13)
    CALL gidriv(ifunc, rbuf, nbuf, chr, lchr, 2)
CASE DEFAULT
  WRITE (msg,'(I10)') idev
  CALL grwarn('Unknown device code in GREXEC: '//msg)
END SELECT
END SUBROUTINE grexec

