import os
import numpy as np
import sys


class linearBRDFBase():
  '''
  Base class of utilities for linearBRDF()
  '''
  def __init__(self,*args,**kwargs):
    '''
    Set base defaults:

    self.verbose  = False
    self.idir     = 'idir'
    self.odir     = 'odir'
    '''
    self.verbose = ('verbose' in kwargs and kwargs['verbose']) or False
    self.idir   =  ('idir' in kwargs and kwargs['idir']) or 'idir'
    self.odir   =  ('odir' in kwargs and kwargs['odir']) or 'odir'
    # check exist
    self.idir_check = False
    self.odir_check = False

  def err(self,*args,**kwargs):
     '''
     Error message: assumed fatal if fatal=True

     All args are printed as error msg

     e.g.

     self.err('that was not correct','so do it again',fatal=True)

     gives:

     that was not correct
     so do it again

     fatal error
     An exception has occurred, use %tb to see the full traceback.

     SystemExit

     To exit: use 'exit', 'quit', or Ctrl-D.

     '''
     import sys
     try:
       print sys.exc_info()[0]
     except:
       pass
     for i in args:
       print i
     if 'fatal' in kwargs and kwargs['fatal']:
         print '\nfatal error'
         raise SystemExit

  def mkdir(self,odir):
    '''
    make directories if dont exist
    '''
    try:
      os.makedirs(odir)
    except:
      pass

