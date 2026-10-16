c
c stores the parameters of different runs in character representation
c
      parameter (msymbol=mpairs+36+10*(mfile-1))

c
c list of symbol names
c
      character*11 symbol(msymbol)

c
c list of symbol values, symbol default values, and save variable used
c to reset a value if an error occurs on trying to set this value
c
      character*256 symbol_value(msymbol),symbol_text(msymbol),symbol_default(msymbol)
     *,symbol_save
c
c short and full filename for inp file
c
     *,fname_inp,fnameinp
c
c short and long run name read from input file - this should be the same as fname_inp
c but without the .extension for all files
c except 2d and 3d plot files where it contains the name of file containing the coefficients, and
c therefore will need to include the .extension
c
     *,fname_run,fnamerun
c
c short and long name for .ato files
c
     *,fname_ato,fnameato
c
c short and long name for .pcof files
c
     *,fname_pcof,fnamepcof

c short and long name for template files

     *,fname_template,fnametemplate

c
c symbol_text gives a brief description of each symbol
c prompt is the shown when reading symbols from the terminal
c dataprompt indicates the subsequent symbols correspond to a particular dataset or distribution
c datatext shows which file or dataset is being displayed
c dataprompttext gives a brief description of this command
c
      character*80 prompt,dataprompt
     *,datatext,dataprompttext

c
c extinp is used to define the input file. It may include a class (e.g. cls for clusters, shm for 
c sharm, etc. EPSR files have no specific class name - the extdot is .inp for these files which is
c a unique qualifier, so there is no need for a specific class.) followed by a .ext 
c (such as .inp for EPSR, .dat for many others, .ato for changeato, 
c and .txt for plot files
c
c extdot is used to show the extension beyond and including the dot.
c
c symbol_setup is used to carry forward which type of file(s) is/are being setup - it is determined
c from the setup command
c

      character*20 extinp,extdot,symbol_setup

c
c type of symbol: 'c' or 'r' or 'i' (character, real, or integer) - this is used to check that
c any input variable to be changed has the specified type.
c
      character*1 symbol_type(msymbol) 

      integer nsymboluptodata,nsymboldata
     *,nsymbol,nsymbolplus
     *,iflag(msymbol),isymboltype
c isymboltype tells setup what kind of setup is being performed.
c EPSR setup has isymboltype=1
     *,ilenfname_inp,ilenfnameinp
     *,ilenfname_run,ilenfnamerun
     *,ilenfname_ato,ilenfnameato
     *,ilenfname_pcof,ilenfnamepcof
     *,ilenfname_template,ilenfnametemplate
     *,nfiles,ilenprompt,ilendataprompt
     *,ilendatatext,ilenextinp
c ilenextinpdot gives the position of the . in extinp
     *,ilenextinpdot,ilenextdot,ilensetup
c variables used in process commands to indicate whether current values have been saved
c or a command is not recognised
     *,isave,inotrecog,irefrho

      common/symbol_c/symbol,symbol_value,symbol_default
     *,symbol_save,symbol_type,symbol_text
     *,fname_inp,fnameinp
     *,fname_run,fnamerun
     *,fname_ato,fnameato
     *,fname_pcof,fnamepcof
     *,fname_template,fnametemplate
     *,prompt,dataprompt,datatext,dataprompttext
c symbol_setup defines the type of symbols being set up
     *,extinp,extdot,symbol_setup

      common/symbol_i/nsymboluptodata,nsymboldata
     *,nsymbol,nsymbolplus
     *,iflag,isymboltype
     *,ilenfname_inp,ilenfnameinp
     *,ilenfname_run,ilenfnamerun
     *,ilenfname_ato,ilenfnameato
     *,ilenfname_pcof,ilenfnamepcof
     *,ilenfname_template,ilenfnametemplate
     *,nfiles,ilenprompt,ilendataprompt,ilendatatext,ilenextinp
     *,ilenextinpdot,ilenextdot,ilensetup
     *,isave,inotrecog,irefrho
