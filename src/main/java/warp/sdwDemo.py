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
from edu.mines.jtk.interp.CubicInterpolator import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.mosaic.PixelsView import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from utils import *
from viewer import *
from warp import *
from warp.ShiftInterp import Method

################################################################################
npp,dpp,fpp = 2000,0.002,0.0
spp = Sampling(npp,dpp,fpp)

# Vp/Vs curve is a cosine function scaled +/- m.
m = 2.0  
a = 0.8  # scale for VpVs curve
b = 0.09 # frequency in Hz for VpVs curve

fpeak = 0.05 # frequency for ricker wavelet

# Constraints for first dimension (converted to r1min/r1max later).
vpvsMin = 1.3 
vpvsMax = 3.0

# Constraints for second dimension (trivial for this example).
r2min,r2max,dg2 = -0.15,0.15,15

# Arrays of coarse sampling intervals for 1D and 2D demo.
dga1D = [1,10,25,50,100,200,400]
dga2D = [1,200]

# Some plotting variables
vLabel = "PP time (s)"
w,h = 800,900

##############################################################################

def main(args):
  go1D()
  #go2D()

def go1D():
  pp,ps,sps,kvpvs,ku,vpvsAvg = makeSequences(m,a,b,spp)
  #pp = Synthetic.addNoise(0.5,fpeak,pp,10+1)
  #ps = Synthetic.addNoise(0.5,fpeak,ps,10+1)

  dw = DynamicWarpingC.fromVpVs(spp,sps,vpvsAvg+0.1,0.0)
  # LINEAR interpolation is more representative of the algorithm.
  # MONOTONIC satisfies the time shifts assumptions for PP-PS registration
  # and are smoother shifts.
  dw.setInterpolationMethod(Method.LINEAR)
  goSparse1(spp,pp,sps,ps,dw,dga1D,ku,kvpvs)

def go2D():
  pp2,ps2,sps,s2,kvpvs2,ku2,vpvsAvg = make2D()
  dw = DynamicWarpingC.fromVpVs(spp,sps,vpvsAvg+0.1,0.0,s2)
  setPlotVars(dw)
  def plot(f,title,s,cbar,cmap):
    plot2(f,title=title,s1=s,vLabel=vLabel,width=w,height=h,cmap1=cmap,
          cbar=cbar,cbw=100)
  plot(pp2,"PP",spp,"Amplitude",gray)
  plot(ps2,"PS",sps,"Amplitude",gray)
  plot(kvpvs2,"Known Vp/Vs",spp,"Vp/Vs ratio",jet)
  plot(ku2,"Known shifts",spp,"Time shift (s)",jet)
  # LINEAR interpolation is more representative of the algorithm.
  # MONOTONIC satisfies the time shifts assumptions for PP-PS registration
  # and are smoother shifts.
  goSparse2(spp,pp2,sps,ps2,dw,dga2D,ku2,Method.LINEAR)

def goSparse1(sf,f,sg,g,dw,dga,ku,kvpvs):
  se = dw.getSampling1()
  ne = se.getCount()
  e = dw.computeErrors(sf,f,sg,g)
  dw.setStrainLimits(0.0,1.0) # Can uncomment strain limits below
  du = {}
  dg1 = {}
  for dg in dga:
    # This is how I would set strain limits. For dg==1 we can't
    # have strain less than 100% in this implementation, otherwise
    # the slope is zero.
    #if dg==1:
    #  rmin,rmax = 0.0,1.0
    #else:
    #  rmin,rmax = (vpvsMin-1.0)/2,(vpvsMax-1.0)/2
    #dw.setStrainLimits(rmin,rmax)
    g1 = Subsample.subsample(ne,dg)
    g1 = Subsample.indicesToSampling(se,g1)
    dg1[dg] = g1
    du[dg] = dw.findShifts(sf,e,g1)
  plot(f,g,dw,ku,kvpvs,du,e,dg1=dg1)

def goSparse2(sf,f,sg,g,dw,dga,ku,interp):
  du = {}
  g1D = {}
  n2 = s2.getCount()
  g2 = Subsample.subsample(n2,dg2)
  g2 = Subsample.indicesToSampling(s2,g2)
  for dg1 in dga:
    if dg1==1:
      dw.setInterpolationMethod(Method.LINEAR)
      r1min,r1max = 0.0,1.0
    else:
      r1min,r1max = (vpvsMin-1.0)/2,(vpvsMax-1.0)/2
      dw.setInterpolationMethod(interp)
    g11 = Subsample.subsample(ne,dg1)
    g11 = Subsample.indicesToSampling(se,g11)
    g1 = zerofloat(ne,n2)
    for i2 in range(n2):
      g1[i2] = copy(g11)
    # print "g1:"; dump(g11)
    # print "g2:"; dump(g2)
    du[dg1] = dw.findShifts(sf,f,sg,g,g1,g2)
    g1D[dg1] = g1
  e = dw.computeErrors2(sf,f,sg,g)
  e = Transpose.transpose312(e)
  WarpUtils.normalizeErrors(e)
  for dg1 in sorted(du):
    desc = " dg1=%d"%dg1
    g1,u = g1D[dg1],du[dg1]
    if dg1==1:
      sigma = 32 # for smoothing of fine grid shifts
      ref = RecursiveExponentialFilter(sigma)
      ref.apply(u,u)
    if dg1!=1:
      plotGrid(f,g1,g2,desc)
    plotShifts(u,desc)
    plotVpvs(f,u,desc)
    plotWarped(sg,g,u,dw,desc)
    # Plot alignment errors with time shift curves. Dashed white line
    # shows the the known shifts, solid white line shows computed shifts.
    p1 = [copy(ne,n2,u),"u","w-",2.0]
    p2 = [copy(ne,n2,ku),"ku","w--",2.0]
    plot3(e,p1=p1,p2=p2,title="Alignemnt Errors"+desc,s1=se,s2=su,
          hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Error",
          clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
          o=x1rx2u)

def makeSequences(m,a,b,spp):
  npp = spp.getCount()
  seed = 1000
  vpvs = VpVsWarp.getVpVs(m,a,b,spp)
  u = VpVsWarp.getU(vpvs,dpp)
  vpvsAvg = VpVsWarp.getAverage(vpvs)
  print "vpvsMin=%g, vpvsMax=%g, vpvsAvg=%g"%(min(vpvs),max(vpvs),vpvsAvg)
  nps = int(npp*(vpvsAvg+1.0)/2.0)
  sps = Sampling(nps,dpp,fpp)
  ps = Synthetic.makeRandomEvents(nps,seed)
  pp = VpVsWarp.warp(ps,u,sps,spp)
  pp = Synthetic.addRickerWavelet(fpeak,pp)
  ps = Synthetic.addRickerWavelet(fpeak,ps)
  return pp,ps,sps,vpvs,u,vpvsAvg

def make2D():
  pp,ps,sps,kvpvs,ku,vpvsAvg = makeSequences(m,a,b,spp)
  n2 = 100
  s2 = Sampling(n2)
  pp2 = zerofloat(npp,n2)
  kvpvs2 = zerofloat(npp,n2)
  ku2 = zerofloat(npp,n2)
  ps2 = zerofloat(sps.getCount(),n2)
  for i2 in range(n2):
    pp2[i2] = copy(pp)
    ps2[i2] = copy(ps)
    kvpvs2[i2] = copy(kvpvs)
    ku2[i2] = copy(ku)
  return pp2,ps2,sps,s2,kvpvs2,ku2,vpvsAvg

def getSparseGridCoords(g1,g2):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = zerofloat(ng1,ng2)
  x2 = zerofloat(ng1,ng2)
  for i2 in range(ng2):
    ig2 = s2.indexOfNearest(g2[i2])
    for i1 in range(ng1):
      x1[i2][i1] = g1[ig2][i1]
      x2[i2][i1] = g2[i2]
  return x1,x2

##############################################################################
# Utilities

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  return div(f,sqrt(g),f)

def normalize(f,nmin,nmax):
  n1 = len(f)
  vmin = min(f)
  vmax = max(f)
  r = vmax-vmin
  nr = nmax-nmin
  for i1 in range(n1):
    vi = f[i1]
    f[i1] = nr*(vi-vmin)/r + nmin

##############################################################################
# Plotting

nearest = Interpolation.NEAREST
linear = Interpolation.LINEAR
x1rx2u = PlotPanel.Orientation.X1RIGHT_X2UP
gray = ColorMap.GRAY
jet = ColorMap.JET

def setPlotVars(dw):
  global ne,nu,n2,se,su,s2,el,ul,xl
  se = dw.getSampling1() # error sampling
  su = dw.getSamplingU() # shift sampling
  s2 = dw.getSampling2()
  ne = se.getCount()
  nu = su.getCount()
  n2 = s2.getCount()
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting

def plot(pp,ps,dw,ku,kvpvs,du,e,dg1):
  se = dw.getSampling1()
  su = dw.getSamplingU()
  ulimits = [su.getFirst(),su.getLast()]
  tlimits = [se.getFirst(),se.getLast()]
  ne = se.getCount()
  pp = copy(ne,pp)
  ku = copy(ne,ku)
  kvpvs = copy(ne,kvpvs)
  ulabel = "Time shift (s)"
  tlabel = "PP time (s)"
  clabel = "Alignment error"
  eclips = [0.0,0.3]
  kua = [ku,"w--",2.0]
  kvpvsa = [kvpvs,"r--",2.0]
  et = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(et)
  normalize(pp,-1.0,1.0)
  plot2P(et,[pp,"k-",1.0],se,su,tlimits,ulimits,vi1=1.0,
         title="Errors")
  plot2P(et,[pp,"k-",1.0],se,su,tlimits,ulimits,vi1=1.0,p01=kua,
         title="Known shifts")
  plot2P(et,kvpvsa,se,su,tlimits,ulimits,vi1=1.0,vl1=[1.0,3.0],
         vLabel1="Vp/Vs",p01=kua,title="Known shifts and Vp/Vs")
  for dg in sorted(du):
    g1,u = dg1[dg],du[dg]
    if dg==1:
      sigma = 32 # for smoothing of fine grid shifts
      ref = RecursiveExponentialFilter(sigma)
      uSmooth = copy(u)
      ref.apply(uSmooth,uSmooth)
      vpvsSmooth = WarpUtils.vpvsBd(spp,uSmooth)
      uSmootha = [copy(ne,uSmooth),"w-",2.0]
      vpvsSmootha = [copy(ne,vpvsSmooth),"k-",2.0]
      plot2P(et,vpvsSmootha,se,su,tlimits,ulimits,vi1=1.0,vl1=[1.0,3.0],
             vLabel1="Vp/Vs",p01=uSmootha,p02=kua,p12=kvpvsa,
             title="dg=%d smoothed"%dg)
    vpvs = WarpUtils.vpvsBd(spp,u)
    ua = [copy(ne,u),"w-",2.0]
    vpvsa = [copy(ne,vpvs),"k-",2.0]
    plot2P(et,vpvsa,se,su,tlimits,ulimits,vi1=1.0,vl1=[1.0,3.0],
           vLabel1="Vp/Vs",p01=ua,p02=kua,p12=kvpvsa,title="dg=%d"%dg)

def plotGrid(pp,g1,g2,desc):
  x1,x2 = getSparseGridCoords(g1,g2)
  plot2(pp,x1=x1,x2=x2,title="PP"+desc,s1=spp,s2=s2,vLabel=vLabel,
        cbar="Amplitude",hLimits=xl,vLimits=el,cbw=100,width=w,height=h)
  
def plotShifts(u,desc):
  us = copy(ne,n2,u)
  plot2(us,title="Time shifts"+desc,s1=se,s2=s2,vLabel=vLabel,
        cbar="Time shift (s)",cmap1=jet,vLimits=el,cbw=100,width=w,
        height=h)
    
def plotVpvs(pp,u,desc):
  vpvs = WarpUtils.vpvs(spp,u)
  vpvs = copy(ne,n2,vpvs)
  print "vpvs %s: min=%g, max=%g"%(desc,min(vpvs),max(vpvs))
  plot2(copy(ne,n2,pp),g=vpvs,title="Vp/Vs"+desc,s1=se,s2=s2,vLabel=vLabel,
        cbar="Vp/Vs ratio",cmap2=jet,hLimits=xl,vLimits=el,cbw=100,width=w,
        height=h)

def plotWarped(sps,ps,u,dw,desc):
  psw = WarpUtils.applyShifts(spp,u,sps,ps)
  # d = sub(psw,pp)
  plot2(psw,title="PS1 warped"+desc,s1=spp,s2=s2,vLabel=vLabel,
        vLimits=el,cbar="Amplitude",cbw=100,width=w,height=h)
  # plot2(d,title="PP-PS1 warped"+desc,s1=spp,s2=s2,vLabel=vLabel,
  #       vLimits=el,cbar="Amplitude",cbw=100,width=w,height=h)
    
def plot2(f,g=None,p=None,p2=None,x1=None,x2=None,s1=None,s2=None,title="",
          hLabel="",vLabel="",cbar=None,cmap1=None,cmap2=None,clips1=None,
          clips2=None,width=None,height=None,hLimits=None,vLimits=None,o=None,
          interp=None,cbw=None):
  if s1==None:
    s1 = Sampling(len(f[0]))
  if s2==None:
    s2 = Sampling(len(f))
  v = Viewer2D(o)
  pv1 = v.addPixels(s1,s2,f,"f")
  if cmap1:
    pv1.setColorModel(cmap1)
  if clips1:
    pv1.setClips(clips1[0],clips1[1])
  if interp:
    pv1.setInterpolation(interp)
  if g:
    pv2 = v.addPixels(s1,s2,g,"g")
    if clips2:
      pv2.setClips(clips2[0],clips2[1])
    if cmap2:
      pv2.setColorModel(cmap2)
  if p:
    pt1 = v.addPoints(p,"p1")
    pt1.setStyle("w-")
    pt1.setLineWidth(2.0)
  if p2:
    pt2 = v.addPoints(p2,"p2")
    pt2.setStyle("w--")
    pt2.setLineWidth(2.0)
  if x1 and x2:
    ptx = v.addPoints(x1,x2,"x1x2")
    ptx.setStyle("rO")
    ptx.setMarkSize(10.0)
  if cbar:
    v.addColorBar(cbar)
  if cbw:
    v.setColorBarWidthMinimum(cbw)
  if width and height:
    v.setSize(width,height)
  if hLimits:
    v.setHLimits(hLimits[0],hLimits[1])
  if vLimits:
    v.setVLimits(vLimits[0],vLimits[1])
  v.setTitle(title)
  v.setHLabel(hLabel)
  v.setVLabel(vLabel)
  v.show()

def plot3(f,g=None,p1=None,p2=None,p3=None,x11=None,x21=None,x13=None,x23=None,
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
  if p1:
    pt = v.addPoints3(p1[0],p1[1])
    pt.setStyle(p1[2])
    pt.setLineWidth(p1[3])
  if p2:
    pt2 = v.addPoints3(p2[0],p2[1])
    pt2.setStyle(p2[2])
    pt2.setLineWidth(p2[3])
  if p3:
    pt3 = v.addPoints3(p3[0],p3[1])
    pt3.setStyle(p3[2])
    pt3.setLineWidth(p3[3])
  if x11 and x21:
    ptx11x21 = v.addPoints3(x11,x21,"x1")
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

def plot2P(e,f,s1,s2,hl,vl0,vl1=None,vLabel1=None,vi0=0.5,vi1=None,vf=None,
           p01=None,p02=None,p12=None,x1x2a=None,x1x2b=None,cmap=jet,
           interp=linear,width=1000,height=840,he0=100,he1=25,we0=100,
           cbw=None,wm=0,title=None):
  panel = PlotPanel(2,1,x1rx2u)
  panel.mosaic.setHeightElastic(0,he0)
  panel.mosaic.setHeightElastic(1,he1)
  panel.mosaic.setWidthElastic(0,we0)
  panel.mosaic.setWidthMinimum(0,wm)
  panel.setHLabel("PP time (s)")
  panel.setVLabel(0,"Time shift (s)")
  panel.setVLabel(1,vLabel1)
  if title:
    panel.setTitle(title)
  if vf:
    panel.setVFormat(vf)
  cbar = panel.addColorBar("Alignment error")
  cbar.setInterval(0.1)
  if cbw:
    cbar.setWidthMinimum(cbw)
  if vi0:
    panel.setVInterval(0,vi0)
  if vi1:
    panel.setVInterval(1,vi1)
  panel.setHLimits(hl[0],hl[1])
  panel.setVLimits(0,vl0[0],vl0[1])
  if vl1:
    panel.setVLimits(1,vl1[0],vl1[1])
  pv = panel.addPixels(0,0,s1,s2,e)
  pv.setClips(0.0,0.3)
  pv.setColorModel(cmap)
  pv.setInterpolation(interp)
  pt11 = panel.addPoints(1,0,s1,f[0])
  pt11.setStyle(f[1])
  pt11.setLineWidth(f[2])
  if p01:
    pt01 = panel.addPoints(0,0,s1,p01[0])
    pt01.setStyle(p01[1])
    pt01.setLineWidth(p01[2])
  if x1x2a:
    ptx1x2a = panel.addPoints(0,0,x1x2a[0],x1x2a[1])
    ptx1x2a.setStyle(x1x2a[2])
    ptx1x2a.setMarkSize(x1x2a[3])
  if p02:
    pt02 = panel.addPoints(0,0,s1,p02[0])
    pt02.setStyle(p02[1])
    pt02.setLineWidth(p02[2])
  if x1x2b:
    ptx1x2b = panel.addPoints(0,0,x1x2b[0],x1x2b[1])
    ptx1x2b.setStyle(x1x2b[2])
    ptx1x2b.setMarkSize(x1x2b[3])
  if p12:
    pt12 = panel.addPoints(1,0,s1,p12[0])
    pt12.setStyle(p12[1])
    pt12.setLineWidth(p12[2])
  frame = PlotFrame(panel)
  frame.setSize(width,height)
  frame.setVisible(True)
  
################################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
