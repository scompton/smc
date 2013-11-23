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

gray = ColorMap.GRAY
jet = ColorMap.JET
bwr = ColorMap.BLUE_WHITE_RED
x1right_x2up = PlotPanel.Orientation.X1RIGHT_X2UP
x1down_x2right = PlotPanel.Orientation.X1DOWN_X2RIGHT
nearest = Interpolation.NEAREST
linear = Interpolation.LINEAR

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

def plot1(fsA,s,x1x2=None,title=None,hLabel=None,vLabel=None,
          o=x1right_x2up,hLimits=None,vLimits=None,vFormat=None,hFormat=None,
          hInterval=None,vInterval=None,width=None,height=None,pngDir=None,
          png1=None,png2=None,pngS=None,ptw=None,cwp=True,fw=1.0,fh=1.0):
  if o==x1right_x2up:
    sp = SimplePlot()
  else:
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if not s:
    s = Sampling(len(f))
  for i in range(len(fsA)):
    fs = fsA[i]
    pv = sp.addPoints(s,fs[0])
    pv.setStyle(fs[1])
  if x1x2:
    pv12 = sp.addPoints(x1x2[0],x1x2[1])
    pv12.setStyle(x1x2[2])
    pv12.setMarkSize(x1x2[3])
  if title:
    sp.setTitle(title)
  if hLabel:
    sp.setHLabel(hLabel)
  if vLabel:
    sp.setVLabel(vLabel)
  if hLimits:
    sp.setHLimits(hLimits[0],hLimits[1])
  if vLimits:
    sp.setVLimits(vLimits[0],vLimits[1])
  if hInterval:
    sp.setHInterval(hInterval)
  if vInterval:
    sp.setVInterval(vInterval)
  if vFormat:
    sp.setVFormat(vFormat)
  if width and height:
    sp.setSize(width,height)
  if pngDir and png1:
    if not ptw:
      if cwp:
        ptw = 222.0
      else:
        ptw = 240.0
    print "ptw=%g"%ptw
    sp.setFontSizeForPrint(8.0,ptw)
    sp.paintToPng(720.0,3.08,pngDir+png1+".png")
  if pngDir and png2:
    if not ptw:
      if cwp:
        ptw = 469.0
      else:
        ptw = 504.0
    print "ptw=%g"%ptw
    sp.setFontSizeForPrint(8.0,ptw)
    sp.paintToPng(720.0,6.51,pngDir+png2+".png")
  if pngDir and pngS:
    sp.setFontSizeForSlide(fw,fh)
    sp.paintToPng(720.0,3.33,pngDir+pngS+".png")

def plot2(f,g=None,pA=None,x12SingleA=None,x12=None,x22=None,s1=None,s2=None,
          title=None,hLabel="Crossline (km)",vLabel="",cbar="Amplitude",
          cbw=None,cmap1=None,cmap2=None,clips1=None,clips2=None,width=600,
          height=800,hLimits=None,vLimits=None,hInterval=None,vInterval=None,
          vFormat=None,o=None,pngDir=None,png1=None,png2=None,pngS=None,
          ptw=None,cwp=True,fw=1.0,fh=1.0,interp=None):
  v = Viewer2D(o)
  if s1==None:
    s1 = Sampling(len(f[0]))
  if s2==None:
    s2 = Sampling(len(f))
  pv1 = v.addPixels(s1,s2,f,"f")
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
  if pA:
    for i in range(len(pA)):
      p = pA[i]
      pt = v.addPoints(p[0],p[1])
      pt.setStyle(p[2])
      pt.setLineWidth(p[3])
  if x12SingleA:
    for i in range(len(x12SingleA)):
      x12Single = x12SingleA[i]
      pt = v.addPoints(x12Single[0],x12Single[1],x12Single[2])
      pt.setStyle(x12Single[3])
      pt.setMarkSize(x12Single[4])
  if x12 and x22:
    ptx12x22 = v.addPoints(x12,x22,"x2")
    ptx12x22.setStyle("rO")
    ptx12x22.setMarkSize(10.0)
  v.setTitle(title)
  v.setHLabel(hLabel)
  v.setVLabel(vLabel)
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
  if hInterval:
    v.setHInterval(hInterval)
  if vInterval:
    v.setVInterval(vInterval)
  if vFormat:
    print "got vFormat"
    v.setVFormat(vFormat)
  if pngDir and png1:
    if not ptw:
      if cwp:
        ptw = 222.0
      else:
        ptw = 240.0
    print "ptw=%g"%ptw
    v.setFontSizeForPrint(8.0,ptw)
    v.paintToPng(720.0,3.08,pngDir+png1+".png")
  if pngDir and png2:
    if not ptw:
      if cwp:
        ptw = 469.0
      else:
        ptw = 504.0
    print "ptw=%g"%ptw
    v.setFontSizeForPrint(8.0,ptw)
    v.paintToPng(720.0,6.51,pngDir+png2+".png")
  if pngDir and pngS:
    v.setFontSizeForSlide(fw,fh)
    v.paintToPng(720.0,3.33,pngDir+pngS+".png")
  v.show()

def plot2P(e,f,s1,s2,hl,vl0,vl1=None,vLabel1=None,vi0=0.5,vi1=None,vf=None,
           p11lw=1.0,p01=None,p02=None,p12=None,x1x2a=None,x1x2b=None,
           clips=[0.0,0.3],cmap=gray,interp=linear,width=1300,height=906,
           he0=100,he1=25,wm=0,cbw=None,pngDir=None,png=None):
  panel = PlotPanel(2,1,x1right_x2up)
  panel.mosaic.setHeightElastic(0,he0)
  panel.mosaic.setHeightElastic(1,he1)
  panel.mosaic.setWidthMinimum(0,wm)
  panel.setHLabel("PP time (s)")
  panel.setVLabel(0,"Time shift (s)")
  panel.setVLabel(1,vLabel1)
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
  pv.setClips(clips[0],clips[1])
  pv.setColorModel(cmap)
  pv.setInterpolation(interp)
  pt11 = panel.addPoints(1,0,s1,f)
  pt11.setLineWidth(p11lw)
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
  frame.setFontSizeForSlide(1.0,0.9)
  frame.setSize(width,height)
  frame.setVisible(True)
  if pngDir and png:
    frame.paintToPng(720.0,3.33,pngDir+png+".png")

def plot3(f,g=None,pA=None,x12SliceA=None,x13=None,x23=None,
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
  if pA:
    for i in range(len(pA)):
      p = pA[i]
      pt = v.addPoints3(p[0],p[1])
      pt.setStyle(p[2])
      pt.setLineWidth(p[3])
  if x12SliceA:
    for i in range(len(x12SliceA)):
      x12Slice = x12SliceA[i]
      ptx11x21 = v.addPoints3(x12Slice[0],x12Slice[1],x12Slice[2])
      ptx11x21.setStyle(x12Slice[3])
      ptx11x21.setMarkSize(x12Slice[4])
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

def plotPP3(f,g=None,x1x2=None,s1=None,s2=None,s3=None,title="",label1="",
            label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
            cmap1=None,cmap2=None,clips1=None,clips2=None,width=900,height=1000,
            limits1=None,limits2=None,limits3=None,vInterval0=None,
            vInterval1=None,hInterval0=None,hInterval1=None,o=None,slices=None,
            pngDir=None,png1=None,png2=None,pngS=None,cbw=None,coordMap=None,
            ptw=None,cwp=True,fw=1.0,fh=1.0,lineColor=None,we0=None,we1=None,
            he0=None,he1=None):
  v = Viewer3P(s1,s2,s3,f,o)
  if we0:
    # v.setWidthElastic(0,we0)
    v.setWidthMinimum(0,we0)
  if we1:
    # v.setWidthElastic(0,we1)
    v.setWidthMinimum(0,we1)
  if he0:
    # v.setHeightElastic(0,he0)
    v.setHeightMinimum(0,he0)
  if he1:
    # v.setHeightElastic(0,he1)
    v.setHeightMinimum(0,he1)
  if slices:
    v.setSlices(slices[0],slices[1],slices[2])
  if g:
    v.addPixels(g)
  if x1x2:
    # v.addPoints(x1x2[0],x1x2[1],x1x2[2],x1x2[3])
    v.addPoints(x1x2[0],x1x2[1])
  if coordMap:
    v.addPoints(coordMap[0],coordMap[1],coordMap[2])
  v.setTitle(title)
  v.setLabel1(label1)
  v.setLabel2(label2)
  v.setLabel3(label3)
#  v.setVFormat(0,"%.0f")
  if vInterval0:
    v.setVInterval(0,vInterval0)
  if vInterval1:
    v.setVInterval(1,vInterval1)
  if hInterval0:
    v.setHInterval(0,hInterval0)
  if hInterval1:
    v.setHInterval(1,hInterval1)
  if clips1:
    v.setClips1(clips1[0],clips1[1])
  if clips2:
    v.setClips2(clips2[0],clips2[1])
  if cbar:
    v.addColorBar(cbar)
  if cbw:
    v.setColorBarWidthMinimum(cbw)
  if cmap1:
    v.setColorModel1(cmap1)
  if cmap2:
    v.setColorModel2(cmap2)
  if lineColor:
    v.setLineColor(lineColor)
  if width and height:
    v.setSize(width,height)
  if limits1:
    v.setLimits1(limits1[0],limits1[1])
  if limits2:
    v.setLimits2(limits2[0],limits2[1])
  if limits3:
    v.setLimits3(limits3[0],limits3[1])
  if pngDir and png1:
    if not ptw:
      if cwp:
        ptw = 222.0 # CWP
      else:
        ptw = 240.0 # Geophysics
    print "ptw=%g"%ptw
    v.setFontSizeForPrint(8.0,ptw)
    v.paintToPng(720.0,3.08,pngDir+png1+".png")
  if pngDir and png2:
    if not ptw:
      if cwp:
        ptw = 469.0 # CWP
      else:
        ptw = 504.0 # Geophysics
    print "ptw=%g"%ptw
    v.setFontSizeForPrint(8.0,ptw)
    v.paintToPng(720.0,6.51,pngDir+png2+".png")
  if pngDir and pngS:
    v.setFontSizeForSlide(fw,fh)
    v.paintToPng(720.0,3.33,pngDir+pngS+".png")
  v.show()

def addImageToWorld(world,image):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  # ipg.setColorModel1(ColorMap.getGray())
  # ipg.setColorModel2(ColorMap.getJet(0.3))
  world.addChild(ipg)
  return ipg

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  # zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  # view.setAxesScale(1.0,1.0,zscale)
  view.setAxesScale(0.75,0.75,1.5)
  view.setScale(1.0)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1250,900)
  frame.setVisible(True)
  return frame

###############################################################################
# Run the function main on the Swing thread
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main))
