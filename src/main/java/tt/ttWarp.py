import sys
import math
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.dsp.DynamicWarping import ErrorExtrapolation
from edu.mines.jtk.interp.CubicInterpolator import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.mosaic.PixelsView import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from fault import Warp1

from utils import *
from viewer import *
from warp import *

##############################################################################

umax = 25
umin = -umax

##############################################################################

def main(args):
  #makeG()
  warpTT()
  #warpTTImages()
  #warpFG()

def warpFG():
  f = getData1("Book1.txt")
  g,s = warp(f)
  plot1(f,"f")
  plot1(g,"g")
  plot1(s,"s")
  dw = DynamicWarping(umin,umax)
  dw.setStrainMax(0.5)
  dw.setErrorExtrapolation(ErrorExtrapolation.REFLECT)
  e = dw.computeErrors(f,g)
  d = dw.accumulateForward(e)
  u = dw.backtrackReverse(d,e)
  et = Transpose.transpose12(e)
  su = Sampling(umax-umin+1,1.0,umin)
  plot2(et,s2=su,title="alignment errors")
  plot2(et,p1=u,p2=s,s2=su,title="alignment errors")

def warpTT():
  f = getData1("Book1.txt")
  g,s = warp(f)
  plot1(f,"f")
  plot1(g,"g")
  plot1(s,"s")
  ftt = getData2("tt.txt")
  gtt = getData2("tt2.txt")
  plot2(ftt,title="ftt")
  plot2(gtt,title="gtt")
  plotPointSlices(ftt,g=gtt,title="tt",i2Label="tau")
  dw = DynamicWarping(-25,25)
  dw.setErrorExtrapolation(ErrorExtrapolation.REFLECT)
  dw.setStrainMax(0.5)
  dw.setShiftSmoothing(2.0)
  gttw = zerofloat(len(gtt[0]),len(gtt))
  ftta = zerofloat(len(gtt[0]),len(gtt))
  for it in range(len(ftt)):
    fa = abs(ftt[it])
    ga = abs(gtt[it])
    u = dw.findShifts(fa,ga)
    #gttw[it] = dw.applyShifts(u,gtt[it])
    gttw[it] = dw.applyShifts(u,ga)
    ftta[it] = fa
  plotPointSlices(ftta,g=gttw,title="tt warped",i2Label="tau")
  plot2(gttw,title="gtt warped")

def warpTTImages():
  f = getData1("Book1.txt")
  g,s = warp(f)
  plot1(f,"f")
  plot1(g,"g")
  plot1(s,"s")
  ftt = getData2("tt.txt")
  gtt = getData2("tt2.txt")
  plot2(ftt,title="ftt")
  plot2(gtt,title="gtt")
  dw = DynamicWarping(-25,25)
  dw.setStrainMax(0.2,0.2)
  dw.setErrorSmoothing(2)
  dw.setShiftSmoothing(2.0,2.0)
  e = dw.computeErrors(ftt,gtt)
  et = Transpose.transpose312(e)
  plot3(et,p=u)
  u = dw.findShifts(ftt,gtt)
  #d = dw.accumulateForward(u)
  #uu = dw.backtrackReverse(d,u)
  h = dw.applyShifts(u,gtt)
  plot2(u,title="u")
  plot2(h,title="gtt warped")
  #plot1(uu,"uu")

def makeG():
  f = getData1("Book1.txt")
  g,s = warp(f)
  plot1(f,"f")
  plot1(g,"g")
  plot1(s,"s")
  writeData(g,"Book2.txt")
 
def warp(f):
  n = len(f)
  w = Warp1.sinusoid(20,n)
  g = w.warp(f)
  s = zerofloat(n)
  for i in range(n):
    s[i] = w.ux(i)
  return g,s

def getData1(fileName):
  tar = TextArrayReader(fileName)
  return tar.getArray()

def getData2(fileName):
  tar = TextArrayReader2(fileName)
  return tar.getArray()

def writeData(f,fileName):
  TextArrayWriter(f,fileName)

##############################################################################
# Plotting

x1rx2u = PlotPanel.Orientation.X1RIGHT_X2UP

def plot1(f,title):
  sp = SimplePlot()
  sp.addPoints(f)
  sp.setTitle(title)

def plot2(f,p1=None,p2=None,s1=None,s2=None,title=""):
  if s1==None:
    s1 = Sampling(len(f[0]))
  if s2==None:
    s2 = Sampling(len(f))
  sp = SimplePlot()
  sp.addPixels(s1,s2,f)
  sp.setTitle(title)
  sp.addColorBar()
  if p1:
    pt = sp.addPoints(s1,p1)
    pt.setStyle("w-")
  if p2:
    pt = sp.addPoints(s1,p2)
    pt.setStyle("w--")

def plot3(f,p=None):
  v = Viewer2D(x1rx2u)
  v.addPixels(f,"f")
  if p:
    pt = v.addPoints(p,"u")
    pt.setStyle("w-")
    pt.setLineWidth(2.0)
  v.show()

def plotPointSlices(f,g=None,title="",i2Label=None):
  v = Viewer2D(x1rx2u)
  v.setTitle(title)
  v.addPoints(f,"f",i2Label)
  if g:
    ptg = v.addPoints(g,"g",i2Label)
    ptg.setStyle("b-")
  v.setSize(800,800)
  v.show()

##############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
