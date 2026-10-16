      subroutine setup_symbols_rmc()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'iresetnumbers.inc'
c
c sets the symbols and their default values for input to EPSR.
c
c
c set up the symbol table for the RMC.inp file
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if ireset >= iflag
c iflag=1 means all accumulators are reset
c iflag=2 means potential coefficients are zerod as well
c
      
      iref=0
      
      iref=iref+1
      symbol(iref)='feedback'
      symbol_type(iref)='r'
      symbol_text(iref)='Confidence factor - should be < 1. [0.8]'
      iflag(iref)=-1
      symbol_default(iref)='0.9'
      
      iref=iref+1
      symbol(iref)='potfac'
      symbol_type(iref)='r'
      symbol_text(iref)='1.0 to enable repulsive potential refinement, 0.0 to inhibit'
      iflag(iref)=0
      symbol_default(iref)='1.0'
      
      iref=iref+1
      symbol(iref)='ereq'
      symbol_type(iref)='r'
      symbol_text(iref)='Overall requested energy amplitude - overrules efilereq'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='sizefactor'
      symbol_type(iref)='r'
      symbol_text(iref)='Multiplying factor for box dimension. [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'
      
      iref=iref+1
      symbol(iref)='nq'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of Q values. [400]'
      iflag(iref)=-2
      symbol_default(iref)='400'
      
      iref=iref+1
      symbol(iref)='qstep'
      symbol_type(iref)='r'
      symbol_text(iref)='Size of Q step [1/A]. [0.05]'
      iflag(iref)=-2
      symbol_default(iref)='0.05'
      
      iref=iref+1
      symbol(iref)='ireset'
      symbol_type(iref)='i'
      symbol_text(iref)='Sets the Empirical Potential to zero'
      iflag(iref)=0
      symbol_default(iref)='1'

      iresetref=iref
      
      iref=iref+1
      symbol(iref)='iinit'
      symbol_type(iref)='i'
      symbol_text(iref)='Sets accumulators to zero. Recalculates r and Q. [1]'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iinitref=iref

      iref=iref+1
      symbol(iref)='ntimes'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of MC cycles between potential refinements. [5]'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='niter'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of potential refinements before exitting. [1]'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iref=iref+1
      symbol(iref)='nsumt'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations already accumulated. [-1 with reset]'
      iflag(iref)=0
      symbol_default(iref)='-1'
      
      iref=iref+1
      symbol(iref)='intra'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of molecule moves between molecule shakes. [100]'
      iflag(iref)=0
      symbol_default(iref)='100'
      
      iref=iref+1
      symbol(iref)='inter'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations in running averages. [5]'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='rho'
      symbol_type(iref)='r'
      symbol_text(iref)='Atomic number density - will be derived from .ato file'
      iflag(iref)=1
      symbol_default(iref)='0.1'
      
      iref=iref+1
      symbol(iref)='cellst'
      symbol_type(iref)='r'
      symbol_text(iref)='Size of r step [A]. [0.03]'
      iflag(iref)=-2
      symbol_default(iref)='0.03'
      
      iref=iref+1
      symbol(iref)='fwhm'
      symbol_type(iref)='r'
      symbol_text(iref)='Resolution width - Q independent term. [0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='fwhmq'
      symbol_type(iref)='r'
      symbol_text(iref)='Resolution width - Q dependent term. [0.02 for SLS]'
      iflag(iref)=0
      symbol_default(iref)='0.02'
      
      iref=iref+1
      symbol(iref)='nsmoop'
      symbol_type(iref)='i'
      symbol_text(iref)='1 means background subtraction is ON, 0 means OFF'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iref=iref+1
      symbol(iref)='fnameato'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of .ato file'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='fnamepcof'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of potential coefficients file.'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='qmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum value of Q used for potential fits. [0.05]'
      iflag(iref)=0
      symbol_default(iref)='0.05'
      
      iqminref=iref

      iref=iref+1
      symbol(iref)='roverlap'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum allowed intermolecular separation between two atoms'
      iflag(iref)=-2
      symbol_default(iref)='0.9'
      iroverlap=iref
      
      iref=iref+1
      symbol(iref)='rminfac'
      symbol_type(iref)='r'
      symbol_text(iref)='Factor to set the minimum separation between pairs. [0.6]'
      iflag(iref)=-2
      symbol_default(iref)='0.6'
      irminfac=iref     

      iref=iref+1
      symbol(iref)='rminpt'
      symbol_type(iref)='r'
      symbol_text(iref)='Radius at which potential truncation gets to 0.5'
      iflag(iref)=-2
      symbol_default(iref)='9'
      
      iref=iref+1
      symbol(iref)='rmaxpt'
      symbol_type(iref)='r'
      symbol_text(iref)='Radius at which potential truncation goes to 0.0'
      iflag(iref)=-2
      symbol_default(iref)='12'
      
      iref=iref+1
      symbol(iref)='rbroad'
      symbol_type(iref)='r'
      symbol_text(iref)='Weighting on chisq in energy function'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='expecf'
      symbol_type(iref)='r'
      symbol_text(iref)='Potential decay for short distance repulsive term. [0.3]'
      iflag(iref)=-2
      symbol_default(iref)='0.3'
      
      iref=iref+1
      symbol(iref)='rcharge'
      symbol_type(iref)='r'
      symbol_text(iref)='Calculates energy due to molecular polarisation. [0.0]'
      iflag(iref)=-2
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='power'
      symbol_type(iref)='r'
      symbol_text(iref)='Repulsive power in Lennard-Jones potential. [12]'
      iflag(iref)=-2
      symbol_default(iref)='12'
      
      iref=iref+1
      symbol(iref)='ncoeffp'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of coefficients used to define the EP. [1]'
      iflag(iref)=-2
      symbol_default(iref)='0'
      incoeffp=iref
      
      iref=iref+1
      symbol(iref)='pdmax'
      symbol_type(iref)='r'
      symbol_text(iref)='Maximum distance of Empirical Potential'
      iflag(iref)=1
      symbol_default(iref)='12.0'
      
      iref=iref+1
      symbol(iref)='pdstep'
      symbol_type(iref)='r'
      symbol_text(iref)='Spacing between coefficients in r. Set by program'
      iflag(iref)=1
      symbol_default(iref)='0.12'
      
      iref=iref+1
      symbol(iref)='npitss'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of steps for refining the potential. [1000]'
      iflag(iref)=0
      symbol_default(iref)='1000'
      
      iref=iref+1
      symbol(iref)='paccept'
      symbol_type(iref)='r'
      symbol_text(iref)='Acceptance factor for potential refinement. [0.0005]'
      iflag(iref)=0
      symbol_default(iref)='0.0005'
      
      iref=iref+1
      symbol(iref)='psigma2'
      symbol_type(iref)='r'
      symbol_text(iref)='Width of Poisson curves. [0.003]'
      iflag(iref)=1
      symbol_default(iref)='0.003'
      
      iref=iref+1
      symbol(iref)='ndata'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of data files to be fit by EPSR'
      iflag(iref)=-1
      symbol_default(iref)='1'
      
      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='datafile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of data file to be fit'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='wtsfile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of weights file for this data set'
      iflag(iref)=-1
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='nrtype'
      symbol_type(iref)='i'
      symbol_text(iref)='Data type - see User Manual for more details'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='rshmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum radius [A] - used for background subtraction'
      iflag(iref)=0
      symbol_default(iref)='0.7'
      
      iref=iref+1
      symbol(iref)='szeros'
      symbol_type(iref)='r'
      symbol_text(iref)='Zero limit - 0 means use first data point for Q=0'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='tweak'
      symbol_type(iref)='r'
      symbol_text(iref)='Scaling factor for this data set. [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'
      
      iref=iref+1
      symbol(iref)='efilereq'
      symbol_type(iref)='r'
      symbol_text(iref)='Requested energy amplitude for this data set [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'

      nsymboldata=iref
      nsymbol=iref
      nsymbolplus=nsymboldata-nsymboluptodata

c
c write the defaults to symbol_value
c
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
      end do

      prompt='Setup rmc>'
      ilenprompt=index(prompt,'>')

c
c set data prompt text
c
      dataprompt='data'
      ilendataprompt=4
      dataprompttext='jumps to the specified dataset'

c
c setup the extension for the input file
c
      extinp='.RMCSR.inp'
      ilenextinp=index(extinp,' ')-1

c
c default number of files
c
      nfiles=1

      return
      end

