#############################################################################
# Figures for CWP slides

from gbcUtils import *

#############################################################################
pngDir = "./pngSlides/"
# pngDir = None
baseDir = getBaseDir()
threeDDir = baseDir+"3D/"
s1,s2,s3 = getSamplings()
n1 = s1.getCount()
n2 = s2.getCount()
n3 = s3.getCount()

# subset for PS data
n1ps = 1501
s1ps = Sampling(n1ps,s1.getDelta(),s1.getFirst())

vpvsAvg = 1.9
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.3,2.5]
iMap = bwr
# iMap = gray
# k1,k2,k3 = 450,75,72
k1,k2,k3 = 348,30,97
slices = [k1,k2,k3]
sigma = 16 # for smoothing of fine grid shifts
plot = True

vpvsMin,vpvsMax,dg1 = 1.4,2.5,80
r2min,r2max,dg2 = -0.15,0.15,15
r3min,r3max,dg3 = -0.15,0.15,15
# Compute slope parameters from vpvsMin/Max.
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

#############################################################################

def main(args):
  pp,ps1 = getImages3DSmooth()
  dw = DynamicWarpingC(pp,ps1,vpvsAvg)
  setPlotVars(dw)
  printConstraints()
  # makePpPs1(pp,ps1)
  # makePs2()
  # makeSAG(dw,pp)
  # make1D(dw,pp,ps1)
  # make1DErrors(pp,ps1)
  # makeFlatPlots(dw,pp)
  # makeREGGridAndVpvs(dw,pp)
  # makeSAGGridAndVpvs(dw,pp)
  makeEnvelopeImage()
  # makeFNEvsSAGWarp(dw,pp,ps1)
  # makeFNEVpvs(dw,pp)
  # makeInterpolate(pp,ps1,dw)

def getImages3DSmooth():
  pp = getGbcImage(baseDir,"pp_smooth")
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  return pp,ps1

def printConstraints():
  k1min = int(math.ceil( r1min*dg1))
  k1max = int(math.floor(r1max*dg1))
  k2min = int(math.ceil( r2min*dg2))
  k2max = int(math.floor(r2max*dg2))
  k3min = int(math.ceil( r3min*dg3))
  k3max = int(math.floor(r3max*dg3))
  info = """Constraints:
  r1min=%g, r1max=%g, dg1=%g, k1min=%g, k1max=%g
  r2min=%g, r2max=%g, dg2=%g, k2min=%g, k2max=%g
  r3min=%g, r3max=%g, dg3=%g, k3min=%g, k3max=%g"""%(r1min,r1max,dg1,k1min,
                                                     k1max,r2min,r2max,dg2,
                                                     k2min,k2max,r3min,r3max,
                                                     dg3,k3min,k3max)
  print info

def makePpPs1(pp,ps1):
  ppSlices = [348,75,72]
  psSlices = [474,75,72]
  cbw=182
  w = 1009
  h = 1024
  he0=320
  plotPP3(pp,s1=s1,s2=s2,s3=s3,clips1=iClips,slices=ppSlices,width=w,height=h,
          label1="PP time (s)",limits1=el,limits2=xl,limits3=yl,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cbar="Amplitude",
          pngDir=pngDir,pngS="pp",he0=he0,cbw=cbw,lineColor=Color.YELLOW)
  plotPP3(ps1,s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=psSlices,width=w,
          height=h,label1="PS time (s)",vInterval1=1.0,vInterval0=2.0,
          hInterval0=2.0,hInterval1=2.0,cbar="Amplitude",pngDir=pngDir,
          pngS="ps1",he0=he0,limits1=[0.0,2.8],limits2=xl,limits3=yl,cbw=cbw,
          lineColor=Color.YELLOW)

def makePs2():
  ps2 = getGbcImage(baseDir,"ps2")
  psSlices = [481,75,72]
  cbw=182
  w = 1009
  h = 1024
  he0=320
  plotPP3(ps2,s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=psSlices,width=w,
          height=h,label1="PS time (s)",vInterval1=1.0,vInterval0=2.0,
          hInterval0=2.0,hInterval1=2.0,cbar="Amplitude",pngDir=pngDir,
          pngS="ps2",he0=he0,limits1=[0.0,2.8],cbw=cbw,lineColor=Color.YELLOW)

def makeSAG(dw,pp):
  label = "sag_ng"
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  dg1 = 50
  ng = 13
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  u = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  print "u: min=%g, max=%g:"%(min(u),max(u))
  u = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                  Interp.MONOTONIC,Interp.LINEAR)
  print "checking shifts:",checkShifts3(u)
  print "u: min=%9.3g, max=%5.3g:"%(min(u),max(u))
  pswSAG = dw.applyShifts(u)
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  zm.apply(0.00,vpvs)
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  u = mul(u,d1)
  u = DynamicWarpingC.extrapolate(u,len(pp[0][0]))
  zm.apply(0.00,u)
  w,h,he0,cbw = 1009,1024,320,182
  slices,vClips,uClips = [348,75,72],[1.3,2.3],None
  plotPP3(u,s1=s1,s2=s2,s3=s3,clips1=uClips,label1="PP time (s)",
          cbar="Time shift (s)",cbw=cbw,cmap1=jet,width=w,height=h,
          limits1=el,slices=slices,vInterval1=1.0,vInterval0=2.0,
          hInterval0=2.0,hInterval1=2.0,pngDir=pngDir,pngS="uSAG",he0=he0)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cbw=cbw,cmap1=jet,width=w,height=h,
          limits1=el,slices=slices,vInterval1=1.0,vInterval0=2.0,
          hInterval0=2.0,hInterval1=2.0,pngDir=pngDir,pngS="vpvsSAG",he0=he0)
  plotPP3(pswSAG,s1=s1,s2=s2,s3=s3,clips1=iClips,cbar="Amplitude",cbw=cbw,
          label1="PP time (s)",limits1=el,slices=slices,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,width=w,height=h,
          pngDir=pngDir,pngS="pswSAG",he0=he0,lineColor=Color.YELLOW)

def make1D(dw,pp,ps1):
  label = "fne"
  uFNE = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  uFNEs = copy(uFNE)
  ref.apply(uFNEs,uFNEs)
  label = "sag_ng"
  # k2,k3 = 75,72
  k2,k3 = 30,78
  dg1 = 50
  ng = 13
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  uSAG = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  uSAG = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,uSAG,x1m,
                                     Interp.MONOTONIC,Interp.LINEAR)
  # uSAG = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,uSAG,x1m,
  #                                    Interp.LINEAR,Interp.LINEAR)
  pswFNE = dw.applyShifts(uFNEs)
  pswSAG = dw.applyShifts(uSAG)
  pswSAG = pswSAG[k3][k2]
  pswFNE = pswFNE[k3][k2]
  pp = pp[k3][k2]
  ps1 = ps1[k3][k2]
  u = mul(uSAG[k3][k2],d1)
  ur = mul(uFNE[k3][k2],d1)
  urs = mul(uFNEs[k3][k2],d1)
  g = [pswSAG,"r-",1.0]
  h = [pswFNE,"b-",1.0]
  p = [u,"SAG","w-",4.0]
  print len(pp),len(ps1),len(pswSAG)
  print "shift at index 342",u[342]
  dw1 = DynamicWarpingC(pp,ps1,vpvsAvg)
  e1 = dw1.computeErrors()
  e1 = DynamicWarpingC.transposeLag(e1)
  DynamicWarpingC.normalizeErrors(e1)
  vs,vl,vf,vi,hi,w,h = "Amplitude",[-1.3,1.3],None,1.0,0.5,1096,404
  # panel = PlotPanel(2,1,x1right_x2up)
  # panel.mosaic.setHeightElastic(0,100)
  # panel.mosaic.setHeightElastic(1,35)
  # panel.setHLabel("PP time (s)")
  # panel.setVLabel(0,"Time shift (s)")
  # panel.setVLabel(1,"f")
  # panel.setVInterval(0,0.2)
  # panel.setHLimits(el[0],el[1])
  # panel.setVLimits(0,ul[0],ul[1])
  # pv = panel.addPixels(0,0,se,su,e1)
  # pv.setClips(0.0,0.15)
  # panel.addPoints(1,0,se,copy(ne1,pp))
  # frame = PlotFrame(panel)
  # frame.setFontSizeForSlide(1.0,0.9)
  # frame.setSize(1372,866)
  # frame.setVisible(True)
  # frame.paintToPng(720.0,7.0,pngDir+"e1f.png")
  # panel2 = PlotPanel(2,1,x1right_x2up)
  # panel2.mosaic.setHeightElastic(0,100)
  # panel2.mosaic.setHeightElastic(1,35)
  # panel2.setHLabel("PP time (s)")
  # panel2.setVLabel(0,"Time shift (s)")
  # panel2.setVLabel(1,"f")
  # panel2.setVInterval(0,0.2)
  # panel2.setHLimits(el[0],el[1])
  # panel2.setVLimits(0,ul[0],ul[1])
  # pv2 = panel2.addPixels(0,0,se,su,e1)
  # pv2.setClips(0.0,0.15)
  # pt2 = panel2.addPoints(0,0,se,u)
  # pt2.setLineWidth(2.0)
  # pt2.setLineColor(Color.WHITE)
  # panel2.addPoints(1,0,se,copy(ne1,pp))
  # frame2 = PlotFrame(panel2)
  # frame2.setFontSizeForSlide(1.0,0.9)
  # frame2.setSize(1372,866)
  # frame2.setVisible(True)
  # frame2.paintToPng(720.0,7.0,pngDir+"e1fu.png")
  # panel3 = PlotPanel(2,1,x1right_x2up)
  # panel3.mosaic.setHeightElastic(0,100)
  # panel3.mosaic.setHeightElastic(1,35)
  # panel3.setHLabel("PP time (s)")
  # panel3.setVLabel(0,"Time shift (s)")
  # panel3.setVLabel(1,"f")
  # panel3.setVInterval(0,0.2)
  # panel3.setHLimits(el[0],el[1])
  # panel3.setVLimits(0,ul[0],ul[1])
  # pv3 = panel3.addPixels(0,0,se,su,e1)
  # pv3.setClips(0.0,0.15)
  # pt3 = panel3.addPoints(0,0,se,ur)
  # pt3.setLineWidth(2.0)
  # pt3.setLineColor(Color.RED)
  # panel3.addPoints(1,0,se,copy(ne1,pp))
  # frame3 = PlotFrame(panel3)
  # frame3.setFontSizeForSlide(1.0,0.9)
  # frame3.setSize(1372,866)
  # frame3.setVisible(True)
  # frame3.paintToPng(720.0,7.0,pngDir+"e1fur.png")
  # panel4 = PlotPanel(2,1,x1right_x2up)
  # panel4.mosaic.setHeightElastic(0,100)
  # panel4.mosaic.setHeightElastic(1,35)
  # panel4.setHLabel("PP time (s)")
  # panel4.setVLabel(0,"Time shift (s)")
  # panel4.setVLabel(1,"f")
  # panel4.setVInterval(0,0.2)
  # panel4.setHLimits(el[0],el[1])
  # panel4.setVLimits(0,ul[0],ul[1])
  # pv4 = panel4.addPixels(0,0,se,su,e1)
  # pv4.setClips(0.0,0.15)
  # pt4 = panel4.addPoints(0,0,se,urs)
  # pt4.setLineWidth(2.0)
  # pt4.setLineColor(Color.CYAN)
  # panel4.addPoints(1,0,se,copy(ne1,pp))
  # frame4 = PlotFrame(panel4)
  # frame4.setFontSizeForSlide(1.0,0.9)
  # frame4.setSize(1372,866)
  # frame4.setVisible(True)
  # frame4.paintToPng(720.0,7.0,pngDir+"e1furs.png")
  plot1(pp,s=s1,vLabel="PP time (s)",vLimits=el,vInterval=hi,o=None,
        hInterval=vi,vFormat=vf,hFormat=" ",width=432,height=w,pngDir=pngDir,
        pngS="pp1D")
  plot1(ps1,s=s1ps,vLabel="PS time (s)",vLimits=[0.0,2.8],vInterval=hi,o=None,
        hInterval=vi,vFormat=vf,hFormat=" ",width=432,height=w,pngDir=pngDir,
        pngS="ps11D")
  plot1(pp,g=g,s=s1,vLabel="PP time (s)",vLimits=el,vInterval=hi,o=None,
        hInterval=vi,vFormat=vf,hFormat=" ",width=432,height=w,pngDir=pngDir,
        pngS="pp_psw1D")
  # plot1(u,s=se,hLabel="PP time (s)",vLabel="Amplitude",hLimits=el,vInterval=vi,
  #       hInterval=hi,vFormat=vf,width=w,height=h,pngDir=pngDir,pngS="u1D")

def make1DErrors(pp,ps1):
  dg1 = 80
  k2,k3 = 28,72
  pp1,ps11 = pp[k3][k2],ps1[k3][k2]
  dw = DynamicWarpingC(pp1,ps11,vpvsAvg)
  e = dw.computeErrors()
  e = DynamicWarpingC.transposeLag(e)
  DynamicWarpingC.normalizeErrors(e)
  g1REG = Subsample.subsample(ne1,dg1)
  ngREG = len(g1REG)
  pp1 = copy(ne1,pp1)
  env = getEnvelope(pp1,ne1);
  g1SAG = Subsample.subsample(env,50,13)
  # g1SAG = Subsample.subsample(env,dg1)
  ngSAG = len(g1SAG)
  uREG = dw.findSparseShifts(r1min,r1max,g1REG)
  uREG = DynamicWarpingC.interpolate(ne1,g1REG,uREG,Interp.LINEAR)
  x1REG = zerofloat(ngREG)
  x2REG = zerofloat(ngREG)
  for i1 in range(ngREG):
    x1REG[i1] = g1REG[i1]
    x2REG[i1] = uREG[g1REG[i1]]
  x2REG = mul(x2REG,d1)
  x1REG = mul(x1REG,d1)
  uREG = mul(uREG,d1)
  uSAG = dw.findSparseShifts(r1min,r1max,g1SAG)
  uSAG = DynamicWarpingC.interpolate(ne1,g1SAG,uSAG,Interp.LINEAR)
  x1SAG = zerofloat(ngSAG)
  x2SAG = zerofloat(ngSAG)
  for i1 in range(ngSAG):
    x1SAG[i1] = g1SAG[i1]
    x2SAG[i1] = uSAG[g1SAG[i1]]
  x2SAG = mul(x2SAG,d1)
  x1SAG = mul(x1SAG,d1)
  uSAG = mul(uSAG,d1)
  psag = [uSAG,"y--",3.0]
  preg = [uREG,"c--",3.0]
  x12SAG = [x1SAG,x2SAG,"yS",14.0]
  x12REG = [x1REG,x2REG,"cO",14.0]
  p12 = [env,"r-",2.0]
  w,h = 1622,866
  vl1 = [-1.5,1.5]
  # vLabel1 = "f & a"
  vLabel1 = None
  plot2P(e,pp1,se,su,el,ul,vl1=vl1,vi0=0.2,vi1=2.0,p11lw=2.0,width=w,height=h,
         pngDir=pngDir,png="ef")
  plot2P(e,pp1,se,su,el,ul,vl1=vl1,vi0=0.2,vi1=2.0,p11lw=2.0,width=w,height=h,
         p01=preg,x1x2a=x12REG,pngDir=pngDir,png="efr")
  plot2P(e,pp1,se,su,el,ul,vLabel1=vLabel1,vl1=vl1,vi0=0.2,vi1=2.0,p11lw=2.0,
         width=w,height=h,p01=preg,x1x2a=x12REG,p12=p12,pngDir=pngDir,
         png="efra")
  plot2P(e,pp1,se,su,el,ul,vLabel1=vLabel1,vl1=vl1,vi0=0.2,vi1=2.0,p11lw=2.0,
         width=w,height=h,p01=preg,x1x2a=x12REG,p02=psag,x1x2b=x12SAG,p12=p12,
         pngDir=pngDir,png="efrsa")

def makeFlatPlots(dw,pp):
  ff = readImage(threeDDir,"ff",ne1,n2,n3) # flattened image
  x1m = readImage(threeDDir,"x1",ne1,n2,n3) # flattened mappings x1(y1,y2,y3)
  fs = readImage(threeDDir,"fs",ne1,n2,n3) # flattening shifts
  f = copy(ne1,n2,n3,pp) # pp with only ne1 samples
  zm = ZeroMask(f)
  zm.apply(0.0,fs)
  fs = mul(fs,d1*1000)
  dg1 = 50
  ng = 13
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  g1f = zeroint(ne1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      g1f[i3][i2] = copy(g1Flat)
  coordMapFlat = Viewer3P.getSparseCoordsMap(g1f,g2,g3,d1,d2,d3)
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  cmf = [coordMapFlat,"rO",6.0]
  cm = [coordMap,"rO",6.0]
  slicesFlat = [348,83,96]
  slices = [345,83,96]
  cbw=182
  w = 1009
  h = 1024
  he0=320
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,
          hInterval1=2.0,lineColor=Color.YELLOW,cbw=cbw,width=w,height=h,
          he0=he0,pngDir=pngDir,pngS="ppForFlat")
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,
          hInterval1=2.0,coordMap=cm,lineColor=Color.YELLOW,cbw=cbw,width=w,
          height=h,he0=he0,pngDir=pngDir,pngS="pp_sag")
  plotPP3(ff,s1=se,s2=s2,s3=s3,label1="PP time (s)",slices=slicesFlat,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,
          limits1=el,limits2=xl,limits3=yl,lineColor=Color.YELLOW,cbw=cbw,
          width=w,height=h,he0=he0,pngDir=pngDir,pngS="ppflat")
  plotPP3(ff,s1=se,s2=s2,s3=s3,label1="PP time (s)",slices=slicesFlat,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,
          limits1=el,limits2=xl,limits3=yl,coordMap=cmf,lineColor=Color.YELLOW,
          cbw=cbw,width=w,height=h,he0=he0,pngDir=pngDir,pngS="ppflat_grid")
  plotPP3(fs,s1=se,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,clips1=[-50,50],cbar="Flattening shift (ms)",slices=slices,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cmap1=jet,
          cbw=cbw,width=w,height=h,he0=he0,pngDir=pngDir,pngS="ppflatshifts")

def makeREGGridAndVpvs(dw,pp):
  label = "reg"
  g1 = zeroint(ne1,n2,n3)
  g11 = Subsample.subsample(ne1,dg1);
  g2 = Subsample.subsample(n2,dg2);
  g3 = Subsample.subsample(n3,dg3);
  for i3 in range(n3):
    for i2 in range(n2):
      g1[i3][i2] = copy(g11)
  u = readImage(threeDDir,"u_"+label,len(g11),len(g2),len(g3))
  print "u: min=%g, max=%g:"%(min(u),max(u))
  # u = DynamicWarpingC.interpolate(ne1,n2,n3,g11,g2,g3,u,False)
  u = DynamicWarpingC.interpolate(ne1,n2,n3,g11,g2,g3,u,None,
                                  Interp.MONOTONIC,Interp.LINEAR)
  print "checking shifts:",checkShifts3(u)
  print "u: min=%9.3g, max=%5.3g:"%(min(u),max(u))
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  cm = [coordMap,"rO",6.0]
  cbw=182
  w = 1009
  h = 1024
  he0=320
  slices = [345,83,96]
  # slices = [345,83,48]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,limits1=el,slices=slices,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cbw=cbw,
          width=w,height=h,he0=he0,pngDir=pngDir,pngS="vpvsREG")
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,
          hInterval1=2.0,coordMap=cm,lineColor=Color.YELLOW,cbw=cbw,width=w,
          height=h,he0=he0,pngDir=pngDir,pngS="pp_reg")

def makeSAGGridAndVpvs(dw,pp):
  label = "sag_ng"
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  dg1 = 50
  ng = 13
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  u = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  print "u: min=%g, max=%g:"%(min(u),max(u))
  uc = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                  Interp.MONOTONIC,Interp.LINEAR)
  print "checking uc:",checkShifts3(uc)
  print "uc: min=%9.3g, max=%5.3g:"%(min(uc),max(uc))
  ul = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                  Interp.LINEAR,Interp.LINEAR)
  print "checking uc:",checkShifts3(ul)
  print "uc: min=%9.3g, max=%5.3g:"%(min(ul),max(ul))
  cbw=182
  w = 1009
  h = 1024
  he0=320
  slices = [345,83,96]
  # slices = [345,83,48]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  vpvsC = DynamicWarpingC.vpvs(uc)
  vpvsC = DynamicWarpingC.extrapolate(vpvsC,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvsC),max(vpvsC))
  zm.apply(0.00,vpvsC)
  plotPP3(vpvsC,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,limits1=el,slices=slices,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cbw=cbw,
          width=w,height=h,he0=he0,pngDir=pngDir,pngS="vpvsSAG")
  # vpvsL = DynamicWarpingC.vpvs(ul)
  # vpvsL = DynamicWarpingC.extrapolate(vpvsL,len(pp[0][0]))
  # print "vpvsSAG: min=%g, max=%g"%(min(vpvsL),max(vpvsL))
  # zm.apply(0.00,vpvsL)
  # plotPP3(vpvsL,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
  #         cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
  #         limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsSAGvREG_L")

def makeEnvelopeImage():
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  e = getEnvelopeImage(ff)
  es = getEnvelopeSum(ff)
  slicesFlat = [348,83,96]
  cbw=182
  w = 1009
  h = 1024
  he0=320
  plotPP3(e,s1=se,s2=s2,s3=s3,label1="PP time (s)",slices=slicesFlat,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,
          limits1=el,limits2=xl,limits3=yl,cbw=cbw,cmap1=jet,clips1=[0.0,1.25],
          lineColor=Color.WHITE,width=w,height=h,he0=he0,pngDir=pngDir,
          pngS="envelope3D")
  s = Sampling(len(es),s1.getDelta(),s1.getFirst())
  plot1(es,s=s,vLabel="PP time (s)",vLimits=el,o=None,
        width=432,height=w,pngDir=pngDir,pngS="envelopeStack")

  
def makeFNEvsSAGWarp(dw,pp,ps1):
  ptw = 222.0
  k2,k3 = 30,78
  slices = [190,30,78]
  label = "fne"
  uFNE = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  uFNEs = copy(uFNE)
  ref.apply(uFNEs,uFNEs)
  label = "sag_ng"
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  dg1 = 50
  ng = 13
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  uSAG = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  uSAG = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,uSAG,x1m,
                                     Interp.MONOTONIC,Interp.LINEAR)
  print "checking shifts:",checkShifts3(uSAG)
  print "uSAG: min=%9.3g, max=%5.3g:"%(min(uSAG),max(uSAG))
  pswFNE = dw.applyShifts(uFNEs)
  pswSAG = dw.applyShifts(uSAG)
  l1,l2,l3 = [0.0,1.0],[0.5,4.0],[0.5,4.0]
  cbw=182
  w = 896
  h = 1024
  he0=380
  plotPP3(pp,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
          label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,slices=slices,cbw=cbw,
          he0=he0,width=w,height=h,pngDir=pngDir,pngS="pp_zoom")
  plotPP3(pswFNE,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
          label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,slices=slices,cbw=cbw,
          he0=he0,width=w,height=h,ptw=ptw,pngDir=pngDir,pngS="pswFNE_zoom")
  plotPP3(pswSAG,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
          label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,slices=slices,cbw=cbw,
          he0=he0,width=w,height=h,ptw=ptw,pngDir=pngDir,pngS="pswSAG_zoom")

def makeFNEVpvs(dw,pp):
  label = "fne"
  u = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  us = copy(u)
  ref.apply(us,us)
  cbw=182
  w = 1009
  h = 1024
  he0=320
  slices = [345,83,96]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(us)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsFNE: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,limits1=el,slices=slices,
          vInterval1=1.0,vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cbw=cbw,
          width=w,height=h,he0=he0,pngDir=pngDir,pngS="vpvsFNE")

def makeInterpolate(pp,ps1,dw):
  label = "sag_ng"
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,50,ng=13)
  u = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  umin = min(u)
  umax = max(u)
  print "u min=%g, max=%g:"%(umin,umax)
  ull = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                    Interp.LINEAR,Interp.LINEAR)
  ulm = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                    Interp.MONOTONIC,Interp.LINEAR)
  uls = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                    Interp.SPLINE,Interp.LINEAR)
  print "Checking shifts ull:",checkShifts3(ull)
  print "Checking shifts ulm:",checkShifts3(ulm)
  print "Checking shifts uls:",checkShifts3(uls)
  print "ull min=%9.3g, max=%5.3g:"%(min(ull),max(ull))
  print "ulm min=%9.3g, max=%5.3g:"%(min(ulm),max(ulm))
  print "uls min=%9.3g, max=%5.3g:"%(min(uls),max(uls))
  vll = DynamicWarpingC.vpvs(ull)
  vlm = DynamicWarpingC.vpvs(ulm)
  vls = DynamicWarpingC.vpvs(uls)
  print "vll min=%8.4g, max=%8.4g:"%(min(vll),max(vll))
  print "vlm min=%8.4g, max=%8.4g:"%(min(vlm),max(vlm))
  print "vls min=%8.4g, max=%8.4g:"%(min(vls),max(vls))
  ull = mul(d1,ull[k3][k2])
  ulm = mul(d1,ulm[k3][k2])
  uls = mul(d1,uls[k3][k2])
  vll = vll[k3][k2]
  vlm = vlm[k3][k2]
  vls = vls[k3][k2]
  gum = [ulm,"r--",3.0]
  hus = [uls,"b-.",3.0]
  gvm = [vlm,"r--",3.0]
  hvs = [vls,"b-.",3.0]
  s1 = Sampling(len(ull),d1,f1)
  vl = [s1.getFirst(),s1.getLast()]
  hl = [1.0,2.5]
  ptw = 222/2.0
  plot1(ull,g=gum,h=hus,lineWidth=2.0,s=s1,vLabel="PP time (s)",
        hLabel="Time shift (s)",o=x1down_x2right,width=480,height=900,
        vLimits=vl,pngDir=pngDir,png1="u_interp",ptw=ptw)
  plot1(vll,g=gvm,h=hvs,lineWidth=2.0,s=s1,vLabel="PP time (s)",
        hLabel="Interval Vp/Vs ratio",o=x1down_x2right,width=480,height=900,
        vLimits=vl,hLimits=hl,pngDir=pngDir,png1="vpvs_interp",ptw=ptw)

def getEnvelopeSum(ff):
  htf = HilbertTransformFilter()
  n3 = len(ff)
  n2 = len(ff[0])
  n1 = len(ff[0][0])
  es = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      x = ff[i3][i2]
      y = copy(x)
      htf.apply(n1,x,y)
      for i1 in range(n1):
        es[i1] = es[i1] + sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return es

def getEnvelopeImage(ff):
  htf = HilbertTransformFilter()
  n3 = len(ff)
  n2 = len(ff[0])
  n1 = len(ff[0][0])
  e = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      x = ff[i3][i2]
      y = copy(x)
      htf.apply(n1,x,y)
      for i1 in range(n1):
        e[i3][i2][i1] = sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return e

def getEnvelope(x,n1):
  htf = HilbertTransformFilter()
  y = copy(x)
  htf.apply(n1,x,y)
  e = zerofloat(n1)
  for i1 in range(n1):
    e[i1] = sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return e

def makeAutomaticGrid(x1m,ff,dw,dg1,ng=None):
  es = getEnvelopeSum(ff);
  # n1 = len(x1m[0])
  # g1t = Subsample.subsample(n1,dg1)
  # g1Flat = Subsample.subsample(es,dg1,len(g1t))
  if ng:
    g1Flat = Subsample.subsample(es,dg1,ng)
  else:
    g1Flat = Subsample.subsample(es,dg1)
  ng = len(g1Flat)
  g1 = zeroint(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        x = x1m[i3][i2][g1Flat[i1]]
        g1[i3][i2][i1] = int(x+0.5)
  g2 = Subsample.subsample(n2,dg2);
  g3 = Subsample.subsample(n3,dg3);
  print "g1Flat:"; dump(g1Flat)
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  return g1,g1Flat,g2,g3

#############################################################################
# Plotting

def setPlotVars(dw):
  global ne1,nel,se,su,el,ul,xl,yl
  ne1 = dw.getPPErrorLength()
  nel = dw.getNumberOfLags()
  se = Sampling(ne1,d1,f1) # error sampling
  su = Sampling(nel,d1,f1) # shift sampling
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  yl = [s3.getFirst(),s3.getLast()] # frame limits for plotting

#############################################################################
run(main)
