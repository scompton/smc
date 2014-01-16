#############################################################################
# Utilities for GBC data

from imports import *

#############################################################################

baseDir = "/data/scompton/gbc/dat/"
n1,d1,f1 = 2000,0.002,0.0
n2,d2,f2 =  150,0.033531,0.0
n3,d3,f3 =  145,0.033531,0.0
global s1,s2,s3
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)

def getSamplings():
  return s1,s2,s3

def getBaseDir():
  return baseDir

###############################################################################
# I/O

def getGbcImage(datDir,filename,n1=n1,n2=n2,n3=n3):
  f = readImage(datDir,filename,n1,n2,n3)
  zm = ZeroMask(f)
  gain(100,f)
  WarpUtils.normalize(f,-1.5,1.5)
  zm.apply(0.0,f)
  return f

def readImage(datDir,fileName,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datDir+fileName+".dat")
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datDir,fileName,x):
  aos = ArrayOutputStream(datDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()

def makeSubset2(datDir,filename,n1s,n1e,n2s,n2e):
  n1S = n1e-n1s+1
  n2S = n2e-n2s+1
  print "makeSubset: n1=%d, n2=%d"%(n1S,n2S)
  f = readImage(datDir,filename,n1,n2)
  g = zerofloat(n1S,n2S)
  copy(n1S,n2S,n1s,n2s,f,0,0,g)
  writeImage(datDir,filename+"_%d_%d"%(n1S,n2S),g)

def makeSubset(datDir,filename,n1s,n1e,n2s,n2e,n3s,n3e):
  n1S = n1e-n1s+1
  n2S = n2e-n2s+1
  n3S = n3e-n3s+1
  print "makeSubset: n1=%d, n2=%d, n3=%d"%(n1S,n2S,n3S)
  f = readImage(datDir,filename,n1,n2,n3)
  g = zerofloat(n1S,n2S,n3S)
  copy(n1S,n2S,n3S,n1s,n2s,n3s,f,0,0,0,g)
  writeImage(datDir,filename+"_%d_%d_%d"%(n1S,n2S,n3S),g)

###############################################################################
# Utility functions

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  return div(f,sqrt(g),f)

def toFloats2(f):
  n2 = len(f)
  n1 = len(f[0])
  ff = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      ff[i2][i1] = float(f[i2][i1])
  return ff

def toFloats3(f):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  ff = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        ff[i3][i2][i1] = float(f[i3][i2][i1])
  return ff

def stretch3(c,f):
  """ stretch (supersample) by specified factor c time sampling of image f """
  si = SincInterp()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  ci = 1.0/c
  g = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      stretch(si,n1,f[i3][i2],ci,g)

def stretch2(c,f):
  """ stretch (supersample) by specified factor c time sampling of image f """
  si = SincInterp()
  n1,n2 = len(f[0]),len(f)
  ci = 1.0/c
  g = zerofloat(n1)
  for i2 in range(n2):
    stretch(si,n1,f[i2],ci,g)

def stretch1(c,f):
  """ stretch (supersample) by specified factor c time sampling of image f """
  si = SincInterp()
  n1 = len(f)
  g = zerofloat(n1)
  stretch(si,n1,f,1.0/c,g)

def stretch(si,n1,f,ci,g):
  si.interpolate(n1,1.0,0.0,f,n1,ci,0.0,g)
  copy(g,f)

def checkShifts(u):
  if not isMonotonic(u):
    print "Interpolated shifts are not monotonic!"

def checkShifts2(u):
  n2 = len(u)
  results = []
  monotonic = True
  for i2 in range(n2):
    if not isMonotonic(u[i2]):
      monotonic = False
      results.append("i2=%d: Interpolated shifts are not monotonic!"%i2)
  if monotonic:
    results.append("Interpolated shifts are monotonic")
  return results

def checkShifts3(u):
  n3 = len(u)
  n2 = len(u[0])
  results = []
  monotonic = True
  for i3 in range(n3):
    for i2 in range(n2):
      if not isMonotonic(u[i3][i2]):
        monotonic = False
        results.append("i2,i3=%d,%d: Interpolated shifts are not monotonic!"%(
          i2,i3))
  if monotonic:
    results.append("Interpolated shifts are monotonic")
  return results

###############################################################################
# Plotting

def displayGBC(ppName="pp",ps1Name="ps1",ps2Name="ps2"):
  pp = getGbcImage( baseDir, ppName)
  ps1 = getGbcImage(baseDir,ps1Name)
  ps2 = getGbcImage(baseDir,ps2Name)
  he0 = 320
  lc = Color.YELLOW
  plotPP3(pp,title=ppName,s1=s1,s2=s2,s3=s3,label1="PP time (s)",he0=he0,
          lineColor=lc)
  plotPP3(ps1,title=ps1Name,s1=s1,s2=s2,s3=s3,label1="PS1 time (s)",he0=he0,
          lineColor=lc)
  plotPP3(ps2,title=ps2Name,s1=s1,s2=s2,s3=s3,label1="PS2 time (s)",he0=he0,
          lineColor=lc)

###############################################################################
# Run the function main on the Swing thread
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main))
