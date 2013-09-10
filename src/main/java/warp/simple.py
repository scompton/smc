import sys

from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from viewer import *

from utils import *
from warp import *

##############################################################################

def main(args):
  goConstant()

def goConstant():
  u1c = 0.014
  uSc = 0.008
  fpeak = 0.06 # cycles/sample

  # Create samplings and data
  sf = Sampling(101,0.002,0.0)
  sg = Sampling(101,0.002,0.0)
  sh = Sampling(111,0.002,0.0)
  fe = Synthetic.makeRandomEvents(sf.getCount(),1)
  si = SincInterp()
  ge = zerofloat(sg.getCount())
  he = zerofloat(sh.getCount())
  si.interpolate(
    sf.getCount(),sf.getDelta(),sf.getFirst(),fe,
    sg.getCount(),sg.getDelta(),-u1c,ge)
  si.interpolate(
    sg.getCount(),sg.getDelta(),sg.getFirst(),ge,
    sh.getCount(),sh.getDelta(),-uSc,he)
  f = Synthetic.addRickerWavelet(fpeak,fe)
  g = Synthetic.addRickerWavelet(fpeak,ge)
  h = Synthetic.addRickerWavelet(fpeak,he)
  # plot1(fe,title="fe")
  # plot1(ge,title="ge")
  # plot1(he,title="he")
  plot1(f,title="f")
  plot1(g,title="g")
  plot1(h,title="h")
  # Sampling for warping
  n1w = sg.getCount()-int(ceil(u1c/sf.getDelta()))+2
  s1 = Sampling(n1w,sf.getDelta(),sf.getFirst())

  # Make min/max shifts
  s1Min = 0.0
  s1Max = u1c+2.0*s1.getDelta()
  sSMin = 0.0
  # sSMax = 0.0
  sSMax = uSc+2.0*s1.getDelta()
  print "s1Max=%f, sSMax=%f"%(s1Max,sSMax)
  dw = DynamicWarpingC3(s1Min,s1Max,sSMin,sSMax,s1)
  su1 = dw.getShiftSampling1()
  suS = dw.getShiftSamplingS()
  dg = 1
  g1i = Subsample.subsample(n1w,dg)
  ng = len(g1i)
  g1 = zerofloat(ng)
  d1 = s1.getDelta()
  for ig in range(ng):
    g1[ig] = g1i[ig]*d1
  u1S = dw.findShifts(sf,f,sg,g,sh,h,g1)
  e = dw.computeErrors(sf,f,sg,g,sh,h)
  DynamicWarpingC.normalizeErrors(e)
  e = DynamicWarpingC.transposeLag(e)
  e = DynamicWarpingC.transposeLag12(e)
  # plot3(e,title="Alignment Errors",s1=s1,s2=su1,hLabel="PP time (s)",
  #       vLabel="Time shift 1 (s)",width=900,height=600,o=x1rx2u)
  plot3D(s1,su1,suS,e)
  print "u1:"; dump(u1S[0])
  print "uS:"; dump(u1S[1])

##############################################################################
# Plotting

x1rx2u = PlotPanel.Orientation.X1RIGHT_X2UP

def plot1(f,g=None,s=None,title=None,hlabel=None,vlabel=None,vLimits=None,
          hLimits=None,w=None,h=None):
  pp = PlotPanel()
  if s==None:
    s = Sampling(len(f))
  pp.addPoints(s,f)
  if g:
    pv2 = pp.addPoints(s,g)
    pv2.setLineColor(Color.BLUE)
  if title:
    pp.setTitle(title)
  if hlabel:
    pp.setHLabel(hlabel)
  if vlabel:
    pp.setVLabel(vlabel)
  if vLimits:
    pp.setVLimits(vLimits[0],vLimits[1])
  if hLimits:
    pp.setHLimits(hLimits[0],hLimits[1])
  pf = PlotFrame(pp)
  if w and h:
    pf.setSize(w,h)
  pf.setVisible(True)

def plot3(f,g=None,p=None,p2=None,p3=None,x11=None,x21=None,x13=None,x23=None,
          s1=None,s2=None,s3=None,title="",hLabel="Crossline (km)",vLabel="",
          cbar="Amplitude",cmap1=None,cmap2=None,clips1=None,clips2=None,
          width=600,height=900,hLimits=None,vLimits=None,o=None,interp=None):
  v = Viewer2D(o)
  if s1==None:
    s1 = Sampling(len(f[0][0]))
  if s2==None:
    s2 = Sampling(len(f[0]))
  if s3==None:
    s3 = Sampling(len(f))
  pv1 = v.addPixels(s1,s2,s3,f,"f")
  if clips1:
    pv1.setClips(clips1[0],clips1[1])
  if cmap1:
    pv1.setColorModel(cmap1)
  if interp:
    pv1.setInterpolation(interp)
  if g:
    pv2 = v.addPixels(g,"g")
    if clips2:
      pv2.setClips(clips2[0],clips2[1])
    if cmap2:
      pv2.setColorModel(cmap2)
  if p:
    pt = v.addPoints(p[0],p[1])
    pt.setStyle(p[2])
    pt.setLineWidth(p[3])
  if p2:
    pt2 = v.addPoints(p2[0],p2[1])
    pt2.setStyle(p2[2])
    pt2.setLineWidth(p2[3])
  if p3:
    pt3 = v.addPoints(p3[0],p3[1])
    pt3.setStyle(p3[2])
    pt3.setLineWidth(p3[3])
  if x11 and x21:
    ptx11x21 = v.addPoints(x11,x21,"x1")
    ptx11x21.setStyle("rO")
    ptx11x21.setMarkSize(10.0)
  if x13 and x23:
    ptx13x23 = v.addPoints3(x13,x23)
    ptx13x23.setStyle("rO")
    ptx13x23.setMarkSize(10.0)
  v.setTitle(title)
  v.setHLabel(hLabel)
  v.setVLabel(vLabel)
  if cbar:
    v.addColorBar(cbar)
  if width and height:
    v.setSize(width,height)
  if hLimits:
    v.setHLimits(hLimits[0],hLimits[1])
  if vLimits:
    v.setVLimits(vLimits[0],vLimits[1])
  v.show()

def plot3D(s1,s2,s3,f):
  world = World()
  frame = SimpleFrame(world)
  ipg = ImagePanelGroup(s1,s2,s3,f)
  ipg.setClips(0.0,0.5)
  ipg.setColorModel(ColorMap.JET)
  world.addChild(ipg)
  frame.setVisible(True)
  
##############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
