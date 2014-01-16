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

from utils import *
from viewer import *
from fault import *

##############################################################################

def main(args):
  goTest()

def goTest():
  r1 = Ray(Point(-8.0,-10.0,0.0),Point(0.0,-5.0,0.0))
  r2 = Ray(Point( 3.0,-1.0,0.0),Point( -4.0,5.0,0.0))
  ip = r1.getIntersection2(r2)
  ir = r1.getIntersection(r2)
  #print ip
  #print ir
  d = ir.getDistanceSqrd()
  print "distance squared: %g"%d
  plot([r1,r2],ip=ip,ir=ir)

def plot(rList,ip=None,ir=None):
  sp = SimplePlot()
  sp.setHLimits(-25.0,25.0)
  sp.setVLimits(-25.0,25.0)
  for r in rList:
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    x1[0] = r.p1.x ; x1[1] = r.p2.x
    x2[0] = r.p1.y ; x2[1] = r.p2.y
    sp.addPoints(x1,x2)
  if ip:
    ptv = sp.addPoints([ip.x],[ip.y])
    ptv.setStyle("bO")
  if ir:
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    x1[0] = ir.p1.x ; x1[1] = ir.p2.x
    x2[0] = ir.p1.y ; x2[1] = ir.p2.y
    ptv = sp.addPoints(x1,x2)
    ptv.setStyle("r-")
    ptv.setLineWidth(4.0)

##############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
