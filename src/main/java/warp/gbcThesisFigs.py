#############################################################################
# Figures for CWP report

from gbcUtils import *

#############################################################################
#pngDir = "/Users/scompton/Pictures/gbcThesis/"
pngDir = None
baseDir = getBaseDir()
threeDDir = baseDir+"3D_thesis/"
s1,s2,s3 = getSamplings()
n1 = s1.getCount()
n2 = s2.getCount()
n3 = s3.getCount()

# subset for PS data
n1ps = 1501
s1ps = Sampling(n1ps,s1.getDelta(),s1.getFirst())

vpvsAvg = 1.9
vpvsMin,vpvsMax = 1.5,2.2
dg1REG = 80
dg1SAG = 50
ng1SAG = 13
r2min,r2max,dg2 = -0.15,0.15,30
r3min,r3max,dg3 = -0.15,0.15,30
# Compute slope parameters from vpvsMin/Max.
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

vClips = [vpvsMin,vpvsMax]
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
iMap = bwr
#iMap = gray
#k1,k2,k3 = 450,75,72
k1,k2,k3 = 348,30,97
slices = [k1,k2,k3]
sigma = 16 # for smoothing of fine grid shifts
he0=320

#############################################################################
   
def main(args):
  pp,ps1 = getImages3DSmooth()
  dw = DynamicWarpingC.fromVpVs(s1,s1ps,vpvsAvg,0.0,s2,s3)
  setGlobals(dw)
  #printConstraints()
  makeREG(s1,pp,dg1REG) # make PP with REG grid image
  makeSAG(s1,pp,dg1SAG,ng1SAG) # make PP with SAG grid image
  goSmooth3D(s1,pp,s1ps,ps1,dw) # find shifts w/ smooth dw
  #go3D(s1,pp,s1ps,ps1,dw) # find shifts w/ orig dw
  #makePpPs1(pp,ps1) # make PP and PS1 image
  #makeFNEvsSAGErrors(dw,pp,ps1)
  #makeFNEvsSAGWarp(pp,ps1)
  #makeFNEVpvs(dw,pp)
  #makeSAGVpvs(dw,pp,ps1)
  #make1DErrors(dw,pp,ps1)

  #makeFlatPlots(dw,pp)
  #makeREGGridAndVpvs(dw,pp)
  #makeSAGGridAndVpvs(dw,pp)
  #makeInterpolate(pp,ps1,dw)

def makeREG(sf,f,dg1):
  g1 = getG1D(se,ne1,dg1)
  g1 = getG3D(g1)
  g2 = getG1D(s2,n2,dg2)
  g3 = getG1D(s3,n3,dg3)
  plotGrid(sf,s2,s3,f,g1,g2,g3)

def makeSAG(sf,f,dg1,ng):
  g1 = getG1SAG(sf.getDelta(),dg1,ng)
  g2 = getG1D(s2,n2,dg2)
  g3 = getG1D(s3,n3,dg3)
  plotGrid(sf,s2,s3,f,g1,g2,g3)

def go3D(sf,f,sg,g,dw):
  dw.setStrainLimits(0,1,-1,1,-1,1)
  g1 = getG1D(se,ne1,1)
  g1 = getG3D(g1)
  g2 = getG1D(s2,n2,1)
  g3 = getG1D(s3,n3,1)
  u = dw.findShifts(sf,f,sg,g,g1,g2,g3)
  checkShifts3(u)
  writeImage(threeDDir,"u_fne",u)

def goSmooth3D(sf,f,sg,g,dw):
  dw.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max)
  dw.setInterpolationMethod(Method.LINEAR)
  dw.setWorkTracker(WarperProgress())
  g2 = getG1D(s2,n2,dg2)
  g3 = getG1D(s3,n3,dg3)
  def goShifts(g1,label):
    u = dw.findShifts(sf,f,sg,g,g1,g2,g3)
    checkShifts3(u)
    writeImage(threeDDir,"u_"+label,u)
  g1r = getG1D(se,ne1,dg1REG)
  g1r = getG3D(g1r)
  g1s = getG1SAG(sf.getDelta(),dg1SAG,ng1SAG)
  goShifts(g1r,"reg_linear")
  goShifts(g1s,"sag_linear")
  dw.setInterpolationMethod(Method.MONOTONIC)
  goShifts(g1r,"reg_monotonic")
  goShifts(g1s,"sag_monotonic")

def getImages3DSmooth():
  pp = getGbcImage(baseDir,"pp_smooth")
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  # plotPP3(ps1,title="PS1",s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=slices,
  #         label1="PS time (s)",vInterval1=0.2,cmap1=iMap,cbw=100)
  return pp,ps1

def makePpPs1(pp,ps1):
  ppSlices = [348,75,72]
  psSlices = [474,75,72]
  ptw = 504.0/2
  cbw=100
  w = 842
  h = 1024
  he0=320
  plotPP3(pp,s1=s1,s2=s2,s3=s3,clips1=iClips,slices=ppSlices,width=w,height=h,
          label1="PP time (s)",limits1=el,limits2=xl,limits3=yl,vInterval1=1.0,
          vInterval0=2.0,hInterval0=2.0,hInterval1=2.0,cbar="Amplitude",
          pngDir=pngDir,png2="pp",ptw=ptw,he0=he0,cbw=cbw,
          lineColor=Color.YELLOW)
  plotPP3(ps1,s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=psSlices,width=w,
          height=h,label1="PS time (s)",vInterval1=1.0,vInterval0=2.0,
          hInterval0=2.0,hInterval1=2.0,cbar="Amplitude",pngDir=pngDir,
          png2="ps1",ptw=ptw,he0=he0,limits1=[0.0,2.8],limits2=xl,limits3=yl,
          cbw=cbw,lineColor=Color.YELLOW)

def makeFNEvsSAGErrors(dw,pp,ps1):
  k2,k3 = 30,78
  uFNE = readImage(threeDDir,"u_fne",n1,n2,n3)
  uSAG = readImage(threeDDir,"u_sag_monotonic",n1,n2,n3)
  uFNE = uFNE[k3][k2]
  uSAG = uSAG[k3][k2]
  ref = RecursiveExponentialFilter(sigma)
  uFNEs = copy(uFNE)
  ref.apply(uFNEs,uFNEs)
  pA = [[copy(ne1,uFNE),"FNE","w-",4.0],
        [copy(ne1,uFNEs),"FNE smooth","c--",4.0],
        [copy(ne1,uSAG),"SAG","y--",4.0]]
  e = dw.computeErrors(s1,pp[k3][k2],s1ps,ps1[k3][k2])
  e = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(e)
  plot2(e,pA=pA,s1=se,s2=su,hLabel="PP time (s)",vInterval=0.1,
        vLabel="Time shift (s)",cbar=None,clips1=[0,0.3],width=1109,height=536,
        hLimits=[0.0,1.0],vLimits=[0.0,0.4],o=x1right_x2up,pngDir=pngDir,
        png2="eFneVsSag",cwp=False)

def makeFNEvsSAGWarp(pp,ps1):
  slices = [190,30,78]
  uFNE = readImage(threeDDir,"u_fne",n1,n2,n3)
  uSAG = readImage(threeDDir,"u_sag_monotonic",n1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  uFNEs = copy(uFNE)
  ref.apply(uFNEs,uFNEs)
  pswFNE = WarpUtils.applyShifts(s1,uFNEs,s1ps,ps1)
  pswSAG = WarpUtils.applyShifts(s1, uSAG,s1ps,ps1)
  l1,l2,l3 = [0.0,1.0],[0.5,4.0],[0.5,4.0]
  ptw = 504/3.0
  def plot(f,pngName):
    plotPP3(f,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
            label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,
            slices=slices,width=778,height=1000,ptw=ptw,he0=he0,pngDir=pngDir,
            png1=pngName)
  plot(pp,"pp_zoom")
  plot(pswFNE,"pswFNE_zoom")
  plot(pswSAG,"pswSAG_zoom")

def makeFNEVpvs(dw,pp):
  u = readImage(threeDDir,"u_fne",n1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  us = copy(u)
  ref.apply(us,us)
  ptw = 504.0/2
  cbw = 100
  w = 830
  h = 1000
  slices = [490,30,78]
  zm = ZeroMask(pp)
  vpvs = WarpUtils.vpvs(s1,us)
  print "vpvsFNE: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,he0=he0,pngDir=pngDir,png1="vpvsFNE")

def makeSAGVpvs(dw,pp,ps1):
  u = readImage(threeDDir,"u_sag_monotonic",n1,n2,n3)
  #u = readImage(threeDDir,"u_sag_linear",n1,n2,n3)
  ptw = 504.0/2
  cbw = 100
  w = 830
  h = 1000
  slices = [490,30,78]
  zm = ZeroMask(pp)
  vpvs = WarpUtils.vpvs(s1,u)
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          ptw=ptw,he0=he0,limits1=el,slices=slices,pngDir=pngDir,png1="vpvsSAG")

def make1DErrors(dw,pp,ps1):
  k2,k3 = 28,72
  pp1,ps11 = pp[k3][k2],ps1[k3][k2]
  r1min = (1.4-1.0)/2
  dw.setStrainLimits(r1min,r1max)
  e = dw.computeErrors(s1,pp1,s1ps,ps11)
  e = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(e)
  g1REG = getG1D(se,ne1,dg1REG)
  ppSub = copy(ne1,pp1)
  env = getEnvelope(ppSub,ne1);
  g1SAG = Subsample.subsample(env,dg1SAG,ng1SAG)
  g1SAG = Subsample.indicesToSampling(s1,g1SAG)
  uREG = dw.findShifts(s1,pp1,s1ps,ps11,g1REG)
  uSAG = dw.findShifts(s1,pp1,s1ps,ps11,g1SAG)
  def getShiftCoords(g1,u):
    n = len(g1)
    x1,x2 = zerofloat(n),zerofloat(n)
    for i in range(n):
      x1[i] = g1[i]
      g1i = s1.indexOfNearest(g1[i])
      x2[i] = u[g1i]
    return x1,x2
  x1REG,x2REG = getShiftCoords(g1REG,uREG)
  x1SAG,x2SAG = getShiftCoords(g1SAG,uSAG)
  h = 600
  w = 1133
  pA = [[copy(ne1,uSAG),"u-sag","y-.",4.0],[copy(ne1,uREG),"u-reg","c--",4.0]]
  x12SingleA = [[x1SAG,x2SAG,"sag","yS",18.0],[x1REG,x2REG,"reg","cO",18.0]]
  fsa = [[ppSub,"k-"],[env,"r-"]]
  plot1(fsa,se,hLabel="PP time (s)",vLabel="Amplitude",hLimits=el,
        vFormat="%11f",width=w,height=400,pngDir=pngDir,png2="envelope",
        cwp=False,vInterval=1.0)
  def plotE(pngName,pA=None,x12SingleA=None):
    plot2(e,pA=pA,x12SingleA=x12SingleA,s1=se,s2=su,hLabel="PP time (s)",
          vLabel="Time shift (s)",clips1=[0,0.3],width=w,height=h,hLimits=el,
          vLimits=ul,o=x1right_x2up,pngDir=pngDir,png2=pngName,cwp=False,
          vInterval=0.2,hInterval=0.2,cbar=None)
  plotE("e_1D")
  plotE("e_1D_sag_reg",pA,x12SingleA)

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
  ptw = 469.0/2
  cbw = 100
  w = 830
  h = 1000
  slices = [348,83,96]
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          slices=slices,pngDir=pngDir,png1="pp",lineColor=Color.YELLOW,ptw=ptw,
          cbw=cbw,width=w,height=h)
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,pngDir=pngDir,png1="pp_sag",
          coordMap=cm,lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,
          height=h)
  plotPP3(ff,s1=se,s2=s2,s3=s3,label1="PP time (s)",slices=slices,pngDir=pngDir,
          png1="ppflat",lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,height=h)
  plotPP3(ff,s1=se,s2=s2,s3=s3,label1="PP time (s)",slices=slices,limits1=el,
          limits2=xl,limits3=yl,pngDir=pngDir,png1="ppflat_grid",
          coordMap=cmf,lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,
          height=h)
  plotPP3(fs,s1=se,s2=s2,s3=s3,label1="PP time (s)",clips1=[-50,50],
          cbar="Flattening shift (ms)",slices=slices,pngDir=pngDir,
          png1="ppflatshifts",cmap1=jet,ptw=ptw,cbw=cbw,width=w,height=h)

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
  u = DynamicWarpingC.interpolate(ne1,n2,n3,g11,g2,g3,u,
                                  Interp.MONOTONIC,Interp.LINEAR)
  print "checking shifts:",checkShifts3(u)
  print "u: min=%9.3g, max=%5.3g:"%(min(u),max(u))
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  cm = [coordMap,"rO",6.0]
  ptw = 469.0/2
  cbw = 100
  w = 830
  h = 1000
  slices = [345,83,48]
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsREG")
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,pngDir=pngDir,png1="pp_reg",
          lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,height=h,coordMap=cm)

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
  # uc = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,False)
  # ul = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,True)
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  cm = [coordMap,"rO",6.0]
  ptw = 469.0/2
  cbw = 100
  w = 830
  h = 1000
  slices = [345,83,48]
  zm = ZeroMask(pp)
  vpvsC = DynamicWarpingC.vpvs(uc)
  vpvsC = DynamicWarpingC.extrapolate(vpvsC,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvsC),max(vpvsC))
  zm.apply(0.00,vpvsC)
  plotPP3(vpvsC,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsSAGvREG")
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,pngDir=pngDir,png1="pp_sagVreg",
          lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,height=h,coordMap=cm)
  ptw = 222.0
  vpvsL = DynamicWarpingC.vpvs(ul)
  vpvsL = DynamicWarpingC.extrapolate(vpvsL,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvsL),max(vpvsL))
  zm.apply(0.00,vpvsL)
  plotPP3(vpvsL,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsSAGvREG_L")

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

def getEnvelope(x,n1):
  htf = HilbertTransformFilter()
  y = copy(x)
  htf.apply(n1,x,y)
  e = zerofloat(n1)
  for i1 in range(n1):
    e[i1] = sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return e

def getG1SAG(d1,dg1,ng):
  x1m = readImage(threeDDir,"x1",ne1,n2,n3) # flattening mappings
  ff = readImage(threeDDir,"ff",ne1,n2,n3) # flattened image
  g1Flat = Subsample.subsampleEnvelope(ff,dg1,ng)
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        x = x1m[i3][i2][g1Flat[i1]]
        g1[i3][i2][i1] = x*d1
  #print "g1Flat:"; dump(g1Flat)
  return g1

def getG1D(s,n,dg):
  g = Subsample.subsample(n,dg); # simple regular interval grid
  return Subsample.indicesToSampling(s,g)

def getG3D(g1D):
  ng1 = len(g1D)
  g = zerofloat(ng1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      g[i3][i2] = copy(g1D) # use the same grid g1D for all traces
  return g

#############################################################################
# Plotting

def setGlobals(dw):
  global ne1,nel,se,su,el,ul,xl,yl
  se = dw.getSampling1()
  su = dw.getSamplingU()
  ne1 = se.getCount()
  nel = su.getCount()
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  yl = [s3.getFirst(),s3.getLast()] # frame limits for plotting

def plotGrid(sf,s2,s3,f,g1,g2,g3):
  coordMap = Viewer3P.getSparseCoordsMap(s1,g1,s2,g2,s3,g3)
  cm = [coordMap,"rO",6.0]
  plotPP3(f,coordMap=cm,title="PP",s1=sf,s2=s2,s3=s3,clips1=iClips,
          label1="PP time (s)",vInterval1=0.2,cbw=100,slices=slices,
          limits1=el,limits2=xl,limits3=yl,he0=he0)


#############################################################################
run(main)
