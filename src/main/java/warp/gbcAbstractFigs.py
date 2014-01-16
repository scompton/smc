#############################################################################
# Figures for CWP report

from gbcUtils import *

#############################################################################
pngDir = "./pngAbstract_3D/"
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
  # makeTeaserFigures(dw,pp)
  makeFNEvsSAGErrors(dw,pp,ps1)
  # makeFNEvsSAGWarp(dw,pp,ps1)
  # makeFNEVpvs(dw,pp)
  # makeSAGVpvs(dw,pp)
  # make1DErrors(pp,ps1)
  # makeFlatPlots(dw,pp)
  # makeREGGridAndVpvs(dw,pp)
  # makeSAGGridAndVpvs(dw,pp)
  # makeInterpolate(pp,ps1,dw)

def getImages3DSmooth():
  pp = getGbcImage(baseDir,"pp_smooth")
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  # plotPP3(ps1,title="PS1",s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=slices,
  #         label1="PS time (s)",vInterval1=0.2,cmap1=iMap,cbw=100)
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

def makeTeaserFigures(dw,pp):
  label = "sag_ng"
  dg1 = 50
  ng = 13
  k3 = 78
  slices = [k1,k2,k3]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  # read sparse shifts from disk
  u = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  print "u: min=%g, max=%g:"%(min(u),max(u))
  # Linear interpolation
  # uLinear = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,True)
  uLinear = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                        Interp.LINEAR,Interp.LINEAR)
  print "checking linear shifts:",checkShifts3(uLinear)
  print "uLinear: min=%9.3g, max=%5.3g:"%(min(uLinear),max(uLinear))
  pswLinear = dw.applyShifts(uLinear)
  vpvsLinear = DynamicWarpingC.vpvs(uLinear)
  vpvsLinear = DynamicWarpingC.extrapolate(vpvsLinear,len(pp[0][0]))
  zm.apply(0.00,vpvsLinear)
  # Cubic interpolation
  uCubic = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,u,x1m,
                                       Interp.MONOTONIC,Interp.LINEAR)
  print "checking monotonic shifts:",checkShifts3(uCubic)
  print "uCubic: min=%9.3g, max=%5.3g:"%(min(uCubic),max(uCubic))
  pswCubic = dw.applyShifts(uCubic)
  vpvsCubic = DynamicWarpingC.vpvs(uCubic)
  vpvsCubic = DynamicWarpingC.extrapolate(vpvsCubic,len(pp[0][0]))
  zm.apply(0.00,vpvsCubic)
  ptw = 469.0/3
  cbw = 114
  plot2(pp[k3],s1=s1,s2=s2,vLabel="PP time (s)",hLimits=xl,vLimits=el,
        pngDir=pngDir,png1="teaser_pp_2D",ptw=ptw,cbw=cbw)
  plot2(pswLinear[k3],s1=s1,s2=s2,vLabel="PP time (s)",vLimits=el,
        pngDir=pngDir,png1="teaser_ps1w_2D_linear",ptw=ptw,cbw=cbw)
  plot2(pswCubic[k3],s1=s1,s2=s2,vLabel="PP time (s)",vLimits=el,
        pngDir=pngDir,png1="teaser_ps1w_2D_cubic",ptw=ptw,cbw=cbw)
  plot2(vpvsLinear[k3],s1=s1,s2=s2,vLabel="PP time (s)",
        cbar="Interval Vp/Vs ratio",cmap1=jet,clips1=vClips,vLimits=el,
        pngDir=pngDir,png1="teaser_vpvs_2D_linear",ptw=ptw,cbw=cbw)
  plot2(vpvsCubic[k3],s1=s1,s2=s2,vLabel="PP time (s)",
        cbar="Interval Vp/Vs ratio",cmap1=jet,clips1=vClips,vLimits=el,
        pngDir=pngDir,png1="teaser_vpvs_2D_cubic",ptw=ptw,cbw=cbw)

def makeFNEvsSAGErrors(dw,pp,ps1):
  label = "fne"
  uFNE = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  uFNEs = copy(uFNE)
  ref.apply(uFNEs,uFNEs)
  ptw = 222.0
  k2,k3 = 30,78
  label = "sag_ng"
  x1m = readImage(threeDDir,"x1",ne1,n2,n3)
  ff = readImage(threeDDir,"ff",ne1,n2,n3)
  dg1 = 50 # minimum interval (samples)
  ng = 13 # number of time sample intervals
  g1,g1Flat,g2,g3 = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng)
  uSAG = readImage(threeDDir,"u_"+label,len(g1Flat),len(g2),len(g3))
  uSAG = DynamicWarpingC.interpolate(ne1,n2,n3,g1Flat,g2,g3,uSAG,x1m,
                                     Interp.MONOTONIC,Interp.LINEAR)
  print "checking shifts:",checkShifts3(uSAG)
  print "uSAG: min=%9.3g, max=%5.3g:"%(min(uSAG),max(uSAG))
  uFNE = mul(uFNE[k3][k2],d1)
  uFNEs = mul(uFNEs[k3][k2],d1)
  uSAG = mul(uSAG[k3][k2],d1)
  # x1SAG = zerofloat(ng)
  # x2SAG = zerofloat(ng)
  # for i1 in range(ng):
  #   x1SAG[i1] = g1[k3][k2][i1]
  #   x2SAG[i1] = uSAG[g1[k3][k2][i1]]
  # x2SAG = mul(x2SAG,d1)
  # x1SAG = mul(x1SAG,d1)
  p = [uFNE,"FNE","w-",4.0]
  p2 = [uFNEs,"r--",3.0]
  p3 = [uSAG,"b-",2.0]
  # x1x2 = [x1SAG,x2SAG,"bO",8.0]
  plot1(uFNE,g=p2,h=p3,lineWidth=2.0,s=se,hLabel="PP time (s)",
        vLabel="Time shift (s)",width=1064,height=490,pngDir=pngDir,
        png1="FneVsSagShifts",ptw=ptw)
  # dw1 = DynamicWarpingC(pp[k3][k2],ps1[k3][k2],vpvsAvg)
  # e1 = dw1.computeErrors()
  # e1 = DynamicWarpingC.transposeLag(e1)
  # DynamicWarpingC.normalizeErrors(e1)
  # plot2(e1,s1=se,s2=su,hLabel="PP time (s)",vInterval=0.1,
  #       vLabel="Time shift (s)",cbar=None,clips1=[0,0.3],width=1064,height=536,
  #       hLimits=[0.0,1.0],vLimits=[0.0,0.4],o=x1rx2u,pngDir=pngDir,
  #       png1="e_plain",ptw=ptw)
  # plot2(e1,s1=se,s2=su,hLabel="PP time (s)",vInterval=0.1,
  #       vLabel="Time shift (s)",cbar=None,clips1=[0,0.3],width=1064,height=536,
  #       hLimits=[0.662,0.758],vLimits=[0.248,0.288],o=x1rx2u,pngDir=pngDir,
  #       png1="eFneVsSag",ptw=ptw)
  # plot2(e1,p=p,p2=p2,p3=p3,s1=se,s2=su,hLabel="PP time (s)",vInterval=0.1,
  #       vLabel="Time shift (s)",cbar=None,clips1=[0,0.3],width=1064,height=536,
  #       hLimits=[0.0,1.0],vLimits=[0.0,0.4],o=x1rx2u,pngDir=pngDir,
  #       png1="eFneVsSag",ptw=ptw)

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
  ppTSlice = zerofloat(n2,n3)
  pswFNETSlice = zerofloat(n2,n3)
  pswSAGTSlice = zerofloat(n2,n3)
  i1 = 190
  for i3 in range(n3):
    for i2 in range(n2):
      ppTSlice[i3][i2] = pp[i3][i2][i1]
      pswFNETSlice[i3][i2] = pswFNE[i3][i2][i1]
      pswSAGTSlice[i3][i2] = pswSAG[i3][i2][i1]
  lh,lv = [0.4,4.3],[1.0,3.0]
  w = 868
  h = 638
  ptw = 222.0
  plot2(ppTSlice,s1=s2,s2=s3,hLabel="Crossline (km)",vLabel="Inline (km)",
        clips1=iClips,cmap1=iMap,width=w,height=h,o=x1rx2u,hLimits=lh,
        vLimits=lv,pngDir=pngDir,png1="pp_tslice",vInterval=1.0,hInterval=1.0,
        cbar=None,ptw=ptw)
  plot2(pswFNETSlice,s1=s2,s2=s3,hLabel="Crossline (km)",
        vLabel="Inline (km)",clips1=iClips,cmap1=iMap,width=w,height=h,
        o=x1rx2u,pngDir=pngDir,png1="pswFNE_tslice",hLimits=lh,vLimits=lv,
        vInterval=1.0,hInterval=1.0,cbar=None,ptw=ptw)
  plot2(pswSAGTSlice,s1=s2,s2=s3,hLabel="Crossline (km)",
        vLabel="Inline (km)",clips1=iClips,cmap1=iMap,width=w,height=h,
        o=x1rx2u,pngDir=pngDir,png1="pswSAG_tslice",hLimits=lh,vLimits=lv,
        vInterval=1.0,hInterval=1.0,cbar=None,ptw=ptw)
  # l1,l2,l3 = [0.0,1.0],[0.5,4.0],[0.5,4.0]
  # w = 648
  # h = 1000
  # ptw = 469/3.0
  # plotPP3(pp,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
  #         label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,
  #         slices=slices,width=w,height=h,ptw=ptw,pngDir=pngDir,png1="pp_zoom")
  # plotPP3(pswFNE,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
  #         label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,
  #         slices=slices,width=w,height=h,ptw=ptw,pngDir=pngDir,
  #         png1="pswFNE_zoom")
  # plotPP3(pswSAG,s1=s1,s2=s2,s3=s3,clips1=iClips,cmap1=iMap,cbar=None,
  #         label1="PP time (s)",limits1=l1,limits2=l2,limits3=l3,
  #         slices=slices,width=w,height=h,ptw=ptw,pngDir=pngDir,
  #         png1="pswSAG_zoom")

def makeFNEVpvs(dw,pp):
  label = "fne"
  u = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  ref = RecursiveExponentialFilter(sigma)
  us = copy(u)
  ref.apply(us,us)
  ptw = 222.0
  cbw = 100
  w = 830
  h = 1000
  slices = [490,30,78]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(us)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsFNE: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsFNE")

def makeSAGVpvs(dw,pp):
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
  ptw = 222.0
  cbw = 100
  w = 830
  h = 1000
  slices = [490,30,78]
  vClips = [1.3,2.3]
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsSAG")
  # plotPP3(pp,g=vpvs,s1=s1,s2=s2,s3=s3,clips2=vClips,label1="PP time (s)",
  #         cbar="Interval Vp/Vs ratio",cmap2=ColorMap.setAlpha(jet,0.5),cbw=cbw,width=w,height=h,
  #         limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsSAG")

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
  ppSub = copy(ne1,pp1)
  env = getEnvelope(ppSub,ne1);
  g1SAG = Subsample.subsample(env,50,13)
  # g1SAG = Subsample.subsample(env,dg1)
  ngSAG = len(g1SAG)
  uREG = dw.findSparseShifts(r1min,r1max,g1REG)
  uREG = DynamicWarpingC.interpolate(ne1,g1REG,uREG,True)
  x1REG = zerofloat(ngREG)
  x2REG = zerofloat(ngREG)
  for i1 in range(ngREG):
    x1REG[i1] = g1REG[i1]
    x2REG[i1] = uREG[g1REG[i1]]
  x2REG = mul(x2REG,d1)
  x1REG = mul(x1REG,d1)
  uREG = mul(uREG,d1)
  # preg = [uREG,"u-reg","w-",2.0]
  # plot2(e,p=preg,x11=x1REG,x21=x2REG,s1=se,s2=su,hLabel="PP time (s)",
  #       vLabel="Time shift (s)",clips1=[0,0.3],width=1096,height=600,
  #       hLimits=el,vLimits=ul,o=x1rx2u,pngDir=pngDir,png1="ae_1D_reg",
  #       vInterval=0.2,hInterval=0.2,cbar=None)
  plot2(e,s1=se,s2=su,hLabel="PP time (s)",
        vLabel="Time shift (s)",clips1=[0,0.3],width=1096,height=600,
        hLimits=el,vLimits=ul,o=x1rx2u,pngDir=pngDir,png1="e_1D",
        vInterval=0.2,hInterval=0.2,cbar=None)
  uSAG = dw.findSparseShifts(r1min,r1max,g1SAG)
  uSAG = DynamicWarpingC.interpolate(ne1,g1SAG,uSAG,True)
  x1SAG = zerofloat(ngSAG)
  x2SAG = zerofloat(ngSAG)
  for i1 in range(ngSAG):
    x1SAG[i1] = g1SAG[i1]
    x2SAG[i1] = uSAG[g1SAG[i1]]
  x2SAG = mul(x2SAG,d1)
  x1SAG = mul(x1SAG,d1)
  uSAG = mul(uSAG,d1)
  psag = [uSAG,"u-sag","y-.",4.0]
  preg = [uREG,"u-reg","c--",4.0]
  x12SAG = [x1SAG,x2SAG,"yS",18.0]
  x12REG = [x1REG,x2REG,"cO",18.0]
  plot1(ppSub,g=env,s=se,hLabel="PP time (s)",vLabel="Amplitude",hLimits=el,
        vFormat="%11f",width=1096,height=400,pngDir=pngDir,png1="envelope")
  # plot2(e,p=psag,p2=preg,x11=x1SAG,x21=x2SAG,s1=se,s2=su,hLabel="PP time (s)",
  #       vLabel="Time shift (s)",clips1=[0,0.3],width=1096,height=600,
  #       hLimits=el,vLimits=ul,o=x1rx2u,pngDir=pngDir,png1="ae_1D_sag",
  #       vInterval=0.2,hInterval=0.2,cbar=None)
  plot2(e,p=psag,p2=preg,x12a=x12SAG,x12b=x12REG,s1=se,s2=su,hLabel="PP time (s)",
        vLabel="Time shift (s)",clips1=[0,0.3],width=1096,height=600,
        hLimits=el,vLimits=ul,o=x1rx2u,pngDir=pngDir,png1="e_1D_sag_reg",
        vInterval=0.2,hInterval=0.2,cbar=None)

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
  cm = [coordMap,"rO",5.0]
  ptw = 469.0/2
  cbw = 70
  w = 570
  h = 700
  slices = [345,83,48]
  vClips = [1.3,2.3]
  plotPP3(pp,s1=s1,s2=s2,s3=s3,label1="PP time (s)",limits1=el,limits2=xl,
          limits3=yl,slices=slices,pngDir=pngDir,png1="pp_reg",
          lineColor=Color.YELLOW,ptw=ptw,cbw=cbw,width=w,height=h,coordMap=cm)
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(vpvs,len(pp[0][0]))
  print "vpvsSAG: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  ptw = 222.0
  plotPP3(vpvs,s1=s1,s2=s2,s3=s3,clips1=vClips,label1="PP time (s)",
          cbar="Interval Vp/Vs ratio",cmap1=jet,cbw=cbw,width=w,height=h,
          limits1=el,slices=slices,ptw=ptw,pngDir=pngDir,png1="vpvsREG")
  
  
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
  cm = [coordMap,"rO",5.0]
  ptw = 469.0/2
  cbw = 70
  w = 570
  h = 700
  slices = [345,83,48]
  vClips = [1.3,2.3]
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
        hLabel="Time shift (s)",o=x1dx2r,width=480,height=900,
        vLimits=vl,pngDir=pngDir,png1="u_interp",ptw=ptw)
  plot1(vll,g=gvm,h=hvs,lineWidth=2.0,s=s1,vLabel="PP time (s)",
        hLabel="Interval Vp/Vs ratio",o=x1dx2r,width=480,height=900,
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
