      parameter (ncommand_max=100,ncommand_line_max=200,nworddim=201,nstsize=3020)
c
c list of commands
c
      character*20 command_text,runflag,runflagglobal
c
c filename to read script commands from. This file is assumed to be in the same directory where
c EPSRrun is running - this may not be the same as that defined by the current_directory. On start
c up commands are read from the input stream channel 5
c
      character*256 script_fname
c
c script status: 1 means script is executing, 0 means not executing
c
      character*20 script_status
c
c text to store commands in and to interrogate
c
      character*3020 sometext

c no. of characters in script file name. Zero signfies no valid script file has been opened

      integer ncharscript
 
c integer used for reading command stream from a file or input

      integer ncommand_stream

c used to specify words in sometext

      integer nwords,ncf(nworddim),ncl(nworddim),ilencommand_text

c common blocks

      common /commands_c/command_text,runflag,runflagglobal,sometext

      common /commands_i/ncommand,nwords,ncf,ncl,ilencommand_text

      common /script_c/script_fname,script_status

      common /script_i/ncharscript

      common /command_line_i/ncommand_stream
