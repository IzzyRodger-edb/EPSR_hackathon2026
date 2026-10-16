      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call changeato()

      end      

      subroutine changeato()
c
c Make changes to a .ato file
c
      use jmk_readline_module
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'symbols.inc'
      include 'setup_changeato.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      call get_history_list()
cc
c setup the changeato symbols
c
      call setup_symbols_changeato()

c
c setup the changeato input file (i.e. the .ato file) for reading
c
      call setup_input_file(10,6,ierrout)
      if(ierrout.ne.0) return
c
c not wanted to leave this file open
c
      close(10)

c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .ato file - need to do this before setting the symbol tables
c
      sizefactor=1.0
      write(6,*) fnameinp(1:ilenfnameinp)
      call minitc(fnameinp,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
      if(ierr.ne.0) then
         write(6,599) fnameinp(1:ilenfnameinp)
599   format(/'Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
      write(6,32) ntype
32    format(/'changeato> There are ',i2,' types of atom in this file'/)
      do it=1,ntype
            write(6,33) it,type(it)
33    format(1x,'Atom type ',i2,' has label ',a3)
      end do

c
c set the number of "files" in this input file
c

      nfiles=ntype

c
c set up the initial values of symbol_value based on the values read from the .ato file
c

      call write_symbols_changeato(1,nsymbol)

c
c read commands from the keyboard and process them
c

      call process_symbols(5,6,ierrout,0)

      call save_history_list()
      return
      end

