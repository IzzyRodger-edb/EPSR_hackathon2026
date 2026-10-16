      program epsrshell
c
c runs the EPSR commands and command script. Also contains some general purpose routines.
c
c 2006-11-03 Several mods - including ability to end shell from a script, causing the 
c script to pause before it exits, so that if the shell is restarted again, it 
c automatically starts executing the script again.
c
c 2006-11-05 Modified to allow character type and size for data plots within GNUplot to be
c specified in the systems_commands.txt
c
c 2006-11-13 Included version in welcome line. Removed convertcluster - points to splotato
c instead.
c
c 2006-11-15 Deprecated splotato and convertcluster. These are now replaced by the single command
c 'plotato'
c
c 2007-08-24 New command 'voids' included to calculate void distribution in .ato file
c
c 2007-10-25 cd ..\xxxx now takes you to the correct directory

      use jmk_readline_module
      include 'system_commands.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'ran1com.inc'
      include 'version.inc'
      character*256 output_string
      integer nchartext

      call get_history_list()
      version=date
      ilenversion=ilendate
      write(6,101) version(1:ilenversion)
101   format(/'EPSRshell> Welcome to EPSR version ',a,
     *//'Type "help" or "?" for a list of commands')
c
c get the system commands for the particular operating system in use
c
      call get_system_commands()

c Add on the path separator to the system_binaries.

      system_binaries=system_binaries(1:ilensystem_binaries)
     *//system_termination(1:ilensystem_termination)
      ilensystem_binaries=ilensystem_binaries+ilensystem_termination

      write(6,600) system_binaries(1:ilensystem_binaries)
600   format(/'Binaries folder: ',a)
c
c set initial directory name - this is called the home directory
c
      call get_directory_name()
      ncharhomedir=nchardir
      home_directory=directory_name
c
c generate a file runflag.txt which will decide whether the run
c will read commands from input stream or from a script file
c
c initialize the text stream
c
      sometext=' '
c
c get the global script status, and start the script if it is not running
c
      ncharscript=0
      call start_script(1)

      do while (runflag.eq.'interactive'.or.runflag.eq.'script')
c
c initialize sometext to spaces
c
         sometext=' '
c
c read a line of text from the appropriate input stream
c
         if(runflag.eq.'interactive') then
            write(6,*)
!            write(6,102,advance='no')
!102   format('EPSRshell> ')
            call readline(sometext,'EPSRshell>') ! read editable input line
!            read(5,100,iostat=ierr) sometext
!            if(ierr.ne.0) then
!               write(6,*) 'Error condition ',ierr,' detected on read'
!               stop
!            endif
         else if(runflag.eq.'script') then
            read(ncommand_stream,100,iostat=ierr) sometext
c
c parse the text to see if it contains any words
c
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c if no. of words is zero rewind the input file and try again
c if this fails, abandon the run
c
            if(nwords.le.0) then
               rewind(ncommand_stream)
               read(ncommand_stream,100,iostat=ierr1) sometext
               if(ierr1.ne.0) then 
                  write(6,*) 'Error condition ',ierr1,' after rewind'
                  stop
               endif
            endif
         endif
100   format(a)
c
c now parse it into separate words
c
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c see if the first character of the first word contains a # - if so 
c we will ignore this line
c
         if(nwords.gt.0) then
            ilenhash=index(sometext(ncf(1):ncl(nwords)),'#')
         else
            ilenhash=1
         endif
c
c if running a script, reconstruct the original text and write it out
c
         if(runflag.eq.'script') then
            nlast=nwords
            call concatenate(1,nlast,output_string,nchartext
     *,nwords,ncf,ncl,sometext)
            write(6,*)
            write(6,*) output_string(1:nchartext)
         endif
         if(nwords.gt.0.and.ilenhash.eq.0) then
c
c find which command was called, if any, and act on it
c
            call process_command()
         endif
c     write(6,*) runflag,runflagglobal
         call check_script_status()
c     write(6,*) runflag,runflagglobal
c
c if the global runflag has been set to interactive, this means the script is to be
c abandoned.
c
         if(runflagglobal.eq.'interactive') then
            if(runflag.eq.'script') then
               close(ncommand_stream)
               ncharscript=0
               runflag='interactive'
            endif
         else
c
c if the current program is running the script and script_status has been set to 
c pausing, stop running the script and go to interactive mode, but leaving open
c the possibility that the script can be restarted again
c
            if(runflag.eq.'script') then
               if(script_status.eq.'pausing'.or.script_status.eq.'ps') then
                  script_status='notrunning'
                  ilenscript_status=10
                  open(10,file='runflag.txt',status='unknown',iostat=ierr)
                  write(10,*,iostat=ierr1) 'script',' ',script_status(1:ilenscript_status)
     *,' ',script_fname(1:ncharscript)
                  close(10)
c
c close the script file
c
                  close(ncommand_stream)
                  runflag='interactive'
                  if(runflagglobal.eq.'exit') runflag='exit'
               else if(runflagglobal.eq.'exit') then
                  open(10,file='runflag.txt',status='unknown',iostat=ierr)
                  write(10,*,iostat=ierr1) 'interactive'
                  close(10)
                  runflag='exit'
               endif
            endif
         endif
      end do
      end

      subroutine process_command()

      use String_Utility
      use jmk_readline_module
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
c     character*2 texttype
      sometext(ncf(1):ncl(1))=StrLowCase(sometext(ncf(1):ncl(1)))
      command_text=sometext(ncf(1):ncl(1))
      ilencommand_text=ncl(1)-ncf(1)+1
c
c save the name of the command being executed
c
      symbol_setup=command_text
      ilensetup=len_trim(command_text) 
      if(command_text.eq.'help'.or.command_text.eq.'?') then
         call help()
      else if(command_text.eq.'system') then
         call system_command()
      else if(command_text.eq.'cd'.or.command_text.eq.'pwd') then
         call set_directory_name()
      else if(command_text.eq.'ls'.or.command_text.eq.'dir') then
         call list_directory()
      else if(command_text.eq.'md') then
         call make_directory()
      else if(command_text.eq.'ed'.or.command_text.eq.'edit') then
         call edit_file()
      else if(command_text.eq.'del') then
         call delete_file()
      else if(command_text.eq.'ss') then
         call start_script(2)
      else if(command_text.eq.'ps') then
         call pause_script()
      else if(command_text.eq.'es') then
         call end_script()
      else if(command_text.eq.'setup') then
         call run_program()
      else if(command_text.eq.'makeato') then
         call run_program()
      else if(command_text.eq.'readgpr') then
         call run_program()
      else if(command_text.eq.'readjmol') then
         call run_program()
      else if(command_text.eq.'readmole') then
         call run_program()
      else if(command_text.eq.'readxyz') then
         call run_program()
      else if(command_text.eq.'makemole') then
         call run_program()
      else if(command_text.eq.'fmole') then
         call run_program()
      else if(command_text.eq.'mixato') then
         call run_program()
      else if(command_text.eq.'addato') then
         call run_program()
      else if(command_text.eq.'makelattice') then
         call run_program()
      else if(command_text.eq.'makelatticeato') then
         call run_program()
      else if(command_text.eq.'dockato') then
         call run_program()
      else if(command_text.eq.'changeato') then
         call run_program()
      else if(command_text.eq.'randomise') then
         call run_program()
      else if(command_text.eq.'epsrwts') then
         call run_program()
      else if(command_text.eq.'nwts') then
         call run_program()
      else if(command_text.eq.'epsrwtsx') then
         call run_program()
      else if(command_text.eq.'xwts') then
         call run_program()
      else if(command_text.eq.'bonds') then
         call run_program()
      else if(command_text.eq.'epsr') then
         call run_program()
      else if(command_text.eq.'readpdb') then
         call run_program()
      else if(command_text.eq.'sharm') then
         call run_program()
      else if(command_text.eq.'sdf') then
         call run_program()
      else if(command_text.eq.'sdfcube') then
         call run_program()
      else if(command_text.eq.'clusters') then
         call run_program()
      else if(command_text.eq.'partials') then
         call run_program()
      else if(command_text.eq.'coord') then
         call run_program()
      else if(command_text.eq.'chains') then
         call run_program()
      else if(command_text.eq.'rings') then
         call run_program()
      else if(command_text.eq.'lsi') then
         call run_program()
      else if(command_text.eq.'voids') then
         call run_program()
      else if(command_text.eq.'bye'.or.command_text.eq.'e'
     *.or.command_text.eq.'q') then
c
c if the current program is running the script, set the script status
c to 'notrunning'. This means when EPSR starts up next time it will automatically
c start the script running again.
c
         if(runflag.eq.'script') then
            script_status='notrunning'
            ilenscript_status=10
            open(10,file='runflag.txt',status='unknown',iostat=ierr)
            write(10,*,iostat=ierr1) 'script',' ',script_status(1:ilenscript_status)
     *,' ',script_fname(1:ncharscript)
            close(10)
c
c close the script file
c
            close(ncommand_stream)
            runflag='interactive'
         endif
         call save_history_list()
         call quit()
      else if(command_text.eq.'plot') then
         call run_program()
      else if(command_text.eq.'plot2d') then
         call run_program()
      else if(command_text.eq.'splot2d') then
         call run_program()
      else if(command_text.eq.'plot3d') then
         call run_program()
      else if(command_text.eq.'plot3djmol') then
         call run_program()
      else if(command_text.eq.'plotato') then
         call run_program()
      else if(command_text.eq.'triangles') then
         call run_program()
      else if(command_text.eq.'torangles') then
         call run_program()
      else if(command_text.eq.'fluctuations') then
         call run_program()
      else if(command_text.eq.'writexyz') then
         call run_program()
      else if(command_text.eq.'runjmol') then
         call run_program()
      else if(command_text.eq.'mapgrentropy') then
         call run_program()
      else if(command_text.eq.'mapgr') then
         call run_program()
      else
         write(6,102) sometext(ncf(1):ncl(1))
 102     format(/'process_command> Command "',a,'" not recognised.'
     *//'  -- Will attempt to run command anyway')
         call run_program()
      endif
      return
      end

      subroutine help()
      include 'commands.inc'
      if(nwords.lt.2) then
         nfirst=1
         if(nwords.eq.0) nfirst=1
         write(6,*)
         write(6,*) 'The list of valid commands is:- '
         write(6,100)
100   format(
     */' setup'
     */' makeato',
     */' runjmol',
     */' readgpr',
     */' readjmol',
     */' readmole',
     */' readxyz',
     */' readpdb',
     */' makemole',
     */' fmole',
     */' bonds',
     */' mixato',
     */' addato',
     */' makelattice',
     */' makelatticeato',
     */' mapgr',
     */' mapgrentropy',
     */' dockato',
     */' changeato',
     */' randomise',
     */' nwts',
     */' xwts',
     *//' setup',
     *//' epsr',
     */' partials',
     */' triangles',
     */' torangles',
     */' chains',
     */' rings',
     */' clusters',
     */' voids',
     */' coord'
     */' sharm',
     */' sdf',
     */' sdfcube',
     */' writexyz',
     *//' plot  -  Starts the plot menu',
     */' plotato  -  Plots the specified .ato file using GNUplot, PGplot, or ATOMS',
     *//' plot2d - plots a 2D density surface using PGPLOT/PGXTAL',
     */' splot2d - plots a 2D density surface using GNUplot',
     */' plot3d - plots a 3D density surface using PGPLOT/PGXTAL',
     */' plot3djmol - Plots a density surface from SHARM/SDF using Jmol',
     *//' system  -  Invoke the specified system command.',
     */' cd  -  Sets the working folder.',
     */' pwd  -  Shows the current working folder.'
     */' ls or dir  -  Lists the contents of the working folder.',
     */' md  -  Makes a new folder.',
     */' ed or edit  -  Invokes the specified editor.',
     */' del  -  Delete a file or files.',
     */' ss  -  Starts an EPSR script file.',
     */' ps  -  Pauses a script if one is running.',
     */' es  -  Ends a script if one is running.',
     *//' bye or e',
     */)
         write(6,99)
99    format(' For more help on a particular command ',
     *'type: help command_name.',
     *//' <> means optional parameter')
      else
c
c evaluate which command is being sought
c
         command_text=sometext(ncf(2):ncl(2))
         ilen=index(command_text,' ')-1   
         if(command_text.eq.'help'.or.command_text.eq.'?') then
            write(6,101)
101   format(/'Lists the valid commands. Type: help or ?'/)
         else if(command_text.eq.'cd') then
            write(6,102)
102   format(/'Sets the working folder. Type: cd <FOLDER_NAME>'/)
         else if(command_text.eq.'pwd') then
            write(6,1021)
1021  format(/'Shows the current working folder. Type: pwd'/)
         else if(command_text.eq.'ls'.or.command_text.eq.'dir') then
            write(6,103)
103   format(/'Lists the contents of the working folder.'
     *,' Type: ls (or dir) <FILENAME>'/)
         else if(command_text.eq.'ss') then
            write(6,104)
104   format(/'Starts an EPSR script file. Type: ss FILENAME'/)
         else if(command_text.eq.'ps') then
            write(6,105)
105   format(/'Pauses a script if one is running. Type: ps'/)
         else if(command_text.eq.'es') then
            write(6,106)
106   format(/'Ends a script if one is running. Type: es'/)
         else if(command_text.eq.'epsr') then
            write(6,107)
107   format(/'Runs EPSR for the specified input files. Type: epsr FILENAME'/)
         else if(command_text.eq.'plot') then
            write(6,108)
108   format(/'Sets up the plot file for GNUPlot. Type PLOT <PLOTTYPE>'/)
         else if(command_text.eq.'bye'.or.command_text.eq.'e') then
            write(6,109)
109   format(/'Exits from EPSR. Type: bye or e'/)
         else if(command_text.eq.'randomise') then
            write(6,110)
110   format(/'Randomizes rotations and positions of molecules.'
     *//'Type randomise <.ato FILENAME>'/)
         else if(command_text.eq.'fmole') then
            write(6,112)
112   format(/'Molecular "shake" to generate molecules. '
     *//'Type fmole <.ATO FILENAME> <NTIMES> '/)
         else if(command_text.eq.'mixato') then
            write(6,113)
113   format(/'Mix .ato files. '
     *//'Type mixato <No. of .ATO files to be mixed> '/)
         else if(command_text.eq.'changeato') then
            write(6,114)
114   format(/'Change parameters in a .ato file. '
     *//'Type changeato <.ato FILENAME> '/)
         else if(command_text.eq.'nwts'.or.command_text.eq.'epsrwts') then
            write(6,115)
115   format(/'Sets up the neutron .wts files from an input filename'
     *//'Type nwts <FILENAME> '/)
         else if(command_text.eq.'makeato') then
            write(6,116)
116   format(/'Makes a .ato file (a single molecule or headgroup). '
     *//'Type makeato <.ato FILENAME> '/)
         else if(command_text.eq.'readgpr') then
            write(6,1161)
1161  format(/'Converts a ghemical .gpr file to .mol file. '
     *//'Type readgpr <.gpr FILENAME> '/)
         else if(command_text.eq.'readjmol') then
            write(6,1162)
1162  format(/'Converts a Jmol .jmol file to .mol file. '
     *//'Type readjmol <.jmol FILENAME> '/)
         else if(command_text.eq.'readmole') then
            write(6,1163)
1163  format(/'Converts a ghemical .gpr or Jmol .jmol file to .mol file. '
     *//'Type readmole <FILENAME>, readmole <.EXT> <FILENAME>, where .EXT=.gpr or .jmol '/)
         else if(command_text.eq.'readxyz') then
            write(6,1164)
1164  format(/'Converts a .xyz file to a .ato file. '
     *//'Type readxyz <FILENAME>, readxyz <FILENAME> <vanderWaalsFactor> '/)
         else if(command_text.eq.'makemole') then
            write(6,117)
117   format(/'Makes an .ato file for a single molecule from a template .mol file'
     *//'Type makemole <.mol FILENAME> ')
         else if(command_text.eq.'setup') then
            write(6,118)
118   format(/'Setup an input file for EPSR or one of the auxiliary programs'
     *//'Type setup <PROGRAM NAME> ')
         else if(command_text.eq.'bonds') then
            write(6,119)
119   format(/'Calculate the bond lengths for each molecule in a .ato file.'
     *//'Type bonds <.ATO FILENAME> '/)
         else if(command_text.eq.'sharm') then
            write(6,121)
121   format(/'Runs SHARM for the specified input file.'
     *//'Type: sharm <FILENAME>'/)
         else if(command_text.eq.'sdf') then
            write(6,1211)
1211  format(/'Runs SDF for the specified input file.'
     *//'Type: sdf <FILENAME>'/)
         else if(command_text.eq.'sdfcube') then
            write(6,488)
488   format(/'Runs SDFCUBE for the specified input file.'
     *//'Type: sdfcube <FILENAME>'/)
         else if(command_text.eq.'clusters') then
            write(6,122)
122   format(/'Runs clusters for the specified input file.'
     *//' Type: clusters <FILENAME>'/)
         else if(command_text.eq.'voids') then
            write(6,1221)
1221  format(/'Runs voids for the specified input file.'
     *//' Type: voids <FILENAME>'/)
         else if(command_text.eq.'plot3d') then
            write(6,123)
123   format(/'Runs plot3d for the specified input file.'
     *//' Type: plot3d <FILENAME>'/)
         else if(command_text.eq.'plot3djmol') then
            write(6,1231)
1231   format(/'Runs plot3djmol for the specified input file.'
     *//' Type: plot3djmol <FILENAME>'/)
         else if(command_text.eq.'triangles') then
            write(6,124)
124   format(/'Runs tri-angles for the specified input file.'
     *//' Type: triangles <FILENAME>'/)
         else if(command_text.eq.'coord') then
            write(6,125)
125   format(/'Runs coord for the specified input file.'
     *//' Type: coord <FILENAME>'/)
         else if(command_text.eq.'chains') then
            write(6,126)
126   format(/'Runs chains for the specified input file.'
     *//' Type: chains <FILENAME>'/)
         else if(command_text.eq.'torangles') then
            write(6,127)
127   format(/'Runs torangles for the specified input file.'
     *//' Type: torangles <FILENAME>'/)
         else if(command_text.eq.'partials') then
            write(6,128)
128   format(/'Runs partials for the specified input file.'
     *//' Type: partials <FILENAME>'/)
         else if(command_text.eq.'md') then
            write(6,129)
129   format(/'Creates a new directory or folder with specified name.'
     *//' Type: md <FOLDERNAME>'/)
         else if(command_text.eq.'ed'.or.command_text.eq.'edit') then
            write(6,130)
130   format(/'Invokes the editor specified in "system_commands.txt".'
     *//' Type: edit <FILENAME>'/)
         else if(command_text.eq.'system') then
            write(6,131)
131   format(/'Invokes the system command in the subsequent words.'
     *//' Type: system <COMMAND> '/)
         else if(command_text.eq.'del') then
            write(6,132)
132   format(/'Delete the specified files.'
     *//' Type: del <FILENAME>'/)
         else if(command_text.eq.'plot2d') then
            write(6,133)
133   format(/'Runs plot2d for the specified input file.'
     *//' Type: plot2d <FILENAME>'/)
         else if(command_text.eq.'splot2d') then
            write(6,551)
551   format(/'Runs splot2d for the specified input file.'
     *//' Type: splot2d <FILENAME>'/)
         else if(command_text.eq.'plotato') then
            write(6,134)
134   format(/'Plots the specified .ato file using GNUplot, PGplot or ATOMS.'
     *//' Type: plotato <FILENAME>'/)
         else if(command_text.eq.'xwts'.or.command_text.eq.'epsrwtsx') then
            write(6,135)
135   format(/'Sets up X-ray form factor coefficients for a .ato file. '
     *//'Type xwts <FILENAME> '/)
         else if(command_text.eq.'rings') then
            write(6,136)
136   format(/'Runs rings for specified input file. '
     *//'Type rings <FILENAME> '/)
         else if(command_text.eq.'dockato') then
            write(6,137)
137   format(/'Docks the first .ato file onto the second. '
     *//'Type dockato '/)
         else if(command_text.eq.'runjmol') then
            write(6,138)
138   format(/'Runs Jmol.'
     *//'Type runjmol '/)
         else if(command_text.eq.'readpdb') then
            write(6,139)
139   format(/'Limited read of .pdb file. Type: readpdb FILENAME'/)
         else if(command_text.eq.'writexyz') then
            write(6,140)
140   format(/'Writes xyz for specified ato file. Type: writexyz FILENAME'/)
         else if(command_text.eq.'mapgrentropy') then
            write(6,141)
141   format(/'Calculates 6D entropy from .SHARM or .SDF files. Type: mapgrentropy FILENAME'/)
         else if(command_text.eq.'mapgr') then
            write(6,142)
142   format(/'Shows g(r) in specified directions (views) from .SHARM or .SDF files. Type: mapgr FILENAME'/)
         else
            write(6,*) 'This command is not implemented.'
            write(6,*)
            write(6,*) 'For a list of valid commands type: help'
            write(6,*)
         endif
      endif
      return
      end

      subroutine system_command()
c
c shows the files in the current input directory
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,text
      integer nchartext
      if(nwords.gt.2) then
c
c concatenate the remaining words into the filename list
c
         nlast=nwords
         call concatenate(3,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
         command=sometext(ncf(2):ncl(2))
     *//' '
     *//text(1:nchartext)
      else if(nwords.gt.1) then

         command=sometext(ncf(2):ncl(2))

      else

         write(6,100)
100   format(/'system_command> No system command specified. ')
         return

      endif

      call system(command,ierr)

      return
      end

      subroutine start_script(nfirst)
c
c loads an EPSR command script from the folder where this program is running 
c (not necessarily from the directory defined by directory_name)
c
      include 'commands.inc'
      integer nfirst
c
c first find out whether a script is already running
c
      call check_script_status()

      if(script_status.ne.'notrunning') then
         write(6,102)
102   format(/'start script> A script may already be running - '
     *,'will run in interactive mode')
         runflag='interactive'
         return
      endif
c
c generate the script file name if one is not defined so far
c
      if(nfirst.gt.1.and.ncharscript.eq.0) call get_script_fname(nfirst)
c
c if no file was generated exit now
c
      if(ncharscript.eq.0) then
         if(nfirst.gt.1) write(6,101)
101   format(/'start script> Script file name not specified')
         runflag='interactive'
         return
      endif
c
c open the current script file
c
      ncommand_stream=19
      open(ncommand_stream,file=script_fname,status='old',iostat=ierr)
      if(ierr.ne.0) then
         write(6,103) script_fname(1:ncharscript)
103   format(/'start_script> Problems with script file ',a)
         ncharscript=0
         runflag='interactive'
         return
      else
c
c set the script running
c
         write(6,104) script_fname(1:ncharscript)
104   format(/'start script> Running script file: ',a)
         script_status='running'
         ilenscript_status=7
c
c set the runflag, and write it to the runflag.txt file
c
         runflag='script'
         runflagglobal='script'
         open(10,file='runflag.txt',status='unknown',iostat=ierr)
         write(10,*,iostat=ierr1) 'script',' ',script_status(1:ilenscript_status)
     *,' ',script_fname(1:ncharscript)
         close(10)
         return
      endif
      end

      subroutine check_script_status()
c
c opens the runflag.txt and determines the global status of EPSR within
c the operating directory 
c
      include 'commands.inc'
c
c need distinct variables for this routine to avoid overwriting sometext
c
      character*256 testtext
      integer nwordst,ncft(100),nclt(100)

c
c first make sure runflag.txt exists - if not create it
c
      ierr=1
      do while(ierr.ne.0)
         open(10,file='runflag.txt',status='old',iostat=ierr)
c
c if there is a problem with this file, create the runflag.txt file
c as interactive
c
         if(ierr.ne.0) then
c     write(6,*) ierr
            open(10,file='runflag.txt',status='unknown',iostat=ierr)
            write(10,*) 'interactive'
            rewind(10)
            ncharscript=0
         endif
      end do
c
c now read some text. We expect to find either 1, 2 or 3 words, depending
c on the value of the first. Any problems with this, then the runflagglobal
c is set to interactive
c
      testtext=' '
      read(10,100,iostat=ierr) testtext
100   format(a)
      close(10)
c
c parse this text to see if contains any words. If 3 or more words,
c and the first is 'script', then set up the script status and script file 
c name.
c
      call parse(testtext,nwordst,ncft,nclt,5,256,100)
      if(nwordst.gt.0) then
         runflagglobal=testtext(ncft(1):nclt(1))
! If runflagglobal is 'exit' then script_status='ps' or 'pausing' will pause the script if runflag='script'
         if(runflagglobal.eq.'exit') then
            if(nwordst.gt.1) script_status=testtext(ncft(2):nclt(2))
         else if(runflagglobal.eq.'script'.and.nwordst.gt.2) then
            script_status=testtext(ncft(2):nclt(2))
            call concatenate(3,nwordst,script_fname,ncharscript
     *,nwordst,ncft,nclt,testtext)
         else
            runflagglobal='interactive'
            ncharscript=0
         endif
      else
         runflagglobal='interactive'
         ncharscript=0
      endif
c
c if no script file was set up - for whatever reason, the local run status
c is set to 'interactive'
c
      if(runflagglobal.eq.'interactive') then
         runflag=runflagglobal
         script_status='notrunning'
         open(10,file='runflag.txt',status='unknown',iostat=ierr)
         write(10,*) 'interactive'
         close(10)
      endif
      return
      end

      subroutine get_script_fname(nfirst)
c
c generates the name of the script file from the supplied words
c ncharscript must have been set to zero before this call will work
c
      include 'commands.inc'
      integer nfirst,nlast
c
c generate the script file name from the supplied words if possible
c
      if(nwords.lt.nfirst) then
         ncharscript=0
         return
      else
         nlast=nwords
         call concatenate(nfirst,nlast,script_fname,ncharscript
     *,nwords,ncf,ncl,sometext)
      endif
      return
      end

      subroutine pause_script()
c
c sets the script_status so that the input commands are read from the terminal
c
      include 'commands.inc'
c
c get the current status of any runs in progress
c
      call check_script_status()
c
c set the runflag, and write it to the runflag.txt file
c
      if(script_status.ne.'running') then
         write(6,*)
         write(6,*) 'No need to pause script'
         return
      endif
      script_status='pausing'
      ilenscript_status=7
      open(10,file='runflag.txt',status='unknown')
      write(10,*) 'script',' ',script_status(1:ilenscript_status)
     *,' ',script_fname(1:ncharscript)
      close(10)
      return
      end

      subroutine end_script()
c
c ends the current script so that it will not be able to start again
c without re-specifying the script file
c
      include 'commands.inc'

      open(10,file='runflag.txt',status='unknown')
      write(10,*) 'interactive'
      close(10)

      return
      end

      subroutine quit()
c
c quits the current program.
c
      include 'commands.inc'
c
c if a script is running in the current program, set it to not running before
c quiting
c
      if(runflag.eq.'script'.and.script_status.eq.'running') then
         call pause_script()
      endif

      stop
      end

      subroutine run_program()
c
c runs the specified program in the directory system_binaries
c
      use String_Utility
      use jmk_readline_module
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,commandsuffix,text,programcall
      integer ncharcommand,ncharcommandsuffix,nchartext,ncharprogramcall
c      write(6,*) 'run_program> number of words read in (including command): ',nwords

      call save_history_list() !Save the current history list for access by other programs
      if(nwords.gt.0) then

c Get the program to execute. This should always be lower case

         command=StrLowCase(command_text)
         ncharcommand=ilencommand_text
         commandsuffix=' '
         ncharcommandsuffix=0
         nchartext=0

c For the command 'setup' Vista gets upset

         if(command(1:ncharcommand).eq.'setup') command='upset'
c Translate other commands
         if(command(1:ncharcommand).eq.'readgpr') then
            command='readmole'
            ncharcommand=len_trim(command)
            commandsuffix='.gpr'
            ncharcommandsuffix=len_trim(commandsuffix)
         else if(command(1:ncharcommand).eq.'readjmol') then
            command='readmole'
            ncharcommand=len_trim(command)
            commandsuffix='.jmol'
            ncharcommandsuffix=len_trim(commandsuffix)
         else if(command(1:ncharcommand).eq.'readms') then
            command='readmole'
            ncharcommand=len_trim(command)
            commandsuffix='.ms'
            ncharcommandsuffix=len_trim(commandsuffix)
         else if(command(1:ncharcommand).eq.'readmsi') then
            command='readmole'
            ncharcommand=len_trim(command)
            commandsuffix='.msi'
            ncharcommandsuffix=len_trim(commandsuffix)
         else if(command(1:ncharcommand).eq.'mapgr') then
            command='mapgrquad'
            ncharcommand=len_trim(command)
         endif
c
c concatenate the all the words into the command line
c
         nlast=nwords
         if(nlast.gt.1) then
            call concatenate(2,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
c            write(6,'(1x,i4,1x,a)') nchartext,text(1:nchartext)
         else
            nchartext=0
         endif

c Add on the directory where the binaries are stored

         programcall=system_binaries(1:ilensystem_binaries)
     *//command(1:ncharcommand)
         ncharprogramcall=ilensystem_binaries+ncharcommand
         if(ilensystem_bin_ext.gt.0) then
            programcall=programcall(1:ncharprogramcall)//system_bin_ext(1:ilensystem_bin_ext)
            ncharprogramcall=ncharprogramcall+ilensystem_bin_ext
         end if

c Add on the arguments to the command line.

         command=programcall(1:ncharprogramcall)//' '
     *//directory_name(1:nchardir)//' '//command(1:ncharcommand)
         ncharcommand=ncharprogramcall+1+nchardir+1+ncharcommand
         if(ncharcommandsuffix.gt.0) then
            command=command(1:ncharcommand)//' '//commandsuffix(1:ncharcommandsuffix)
            ncharcommand=ncharcommand+1+ncharcommandsuffix
         endif
         if(nchartext.gt.0) then
            command=command(1:ncharcommand)//' '//text(1:nchartext)
            ncharcommand=ncharcommand+1+nchartext
         endif

         write(6,100) command(1:ncharcommand)
100   format(/a)

         call system(command,ierr)

      endif
      call get_history_list() !Get the current list - possibly modified by program.
      return
      end




