##############################################################################
# Simultaneous warping figures for thesis

from imports import *
from plotUtils import *

##############################################################################
pngDir = "/Users/scompton/Pictures/simwThesis/"
#pngDir=None
nf,df,ff = 2000,0.002,0.0
# nf,df,ff = 1000,0.002,0.0
# nf,df,ff = 500,0.002,0.0
sf = Sampling(nf,df,ff)

n2 = 500
nrms = 1.0 # rms noise/signal ratio
fpeak = 0.05

m1 =  2.0 # centered Vp/Vs value
a1 =  0.8 # scale for Vp/Vs curve
b1 = 0.09 # frequency in Hz for Vp/Vs curve

# m2 =  2.2 # centered Vp/Vs value
# a2 =  0.8 # scale for Vp/Vs curve
# b2 = 0.09 # frequency in Hz for Vp/Vs curve

aS =  0.05 # scale for shift curve
bS = 0.04 # frequency in Hz for shift curve

# vpvs1Min,vpvs1Max,vpvs1Avg = 1.3,3.0,2.3
vpvs1Min,vpvs1Max,vpvs1Avg = 1.3,3.0,2.4
vpvs2Min,vpvs2Max,vpvs2Avg = 1.3,3.5,2.4
rSMin,rSMax,uSMin,uSMax = 0.0,0.02,0.0,0.08
#rSMin,rSMax,uSMin,uSMax = 0.0,0.02,0.0,0.2

sw = Stopwatch()

#bgColor = Color.BLACK
bgColor = None
pw = False # plot warped traces

##############################################################################

def main(args):
  #makeSyntheticTraces()
  #makeCompress()
  makeIndividual(200,200,200)
  #makeIndividual(1,1,1)
  #makeSimultaneous(200)

def makeSyntheticTraces():
  init(1.0)
  ptw = 504.0
  def plot(fA,s,vl,pngName):
    plot1(fA,s=s,hLabel=vl,vLabel="Amplitude",vLimits=[-2.0,2.0],width=1000,
          height=200,o=x1rx2u,pngDir=pngDir,png2=pngName,ptw=ptw)
  plot([[h,"k-"]],sh,"PS2 time (s)","simw_h")
  plot([[g,"k-"]],sg,"PS1 time (s)","simw_g")
  plot([[f,"k-"]],sf,"PP time (s)" ,"simw_f")

def makeCompress():
  init(0.0)
  e, u, dw = goPpPs1( 200)
  ec,u,dwc = goPpPs1C(200)
  c = dwc.getCompression()
  ku1c = WarpUtils.getScaledShifts(sf,ku1,1.0/c)
  ptw = 504.0/2
  def plot(dw,e,ku,pngName):
    se = dw.getSampling1()
    su = dw.getSamplingU()
    el = [se.getFirst(),se.getLast()]
    ul = [su.getFirst(),su.getLast()]
    e = Transpose.transpose12(e)
    WarpUtils.normalizeErrors(e)
    uA = [[copy(se.getCount(),ku),"ku","w--",2.0]]
    plot2(e,pA=uA,s1=se,s2=su,hLabel="PP time (s)",vLabel="Time shift (s)",
          cbar=None,clips1=[0,0.3],width=750,height=600,hLimits=el,vLimits=ul,
          o=x1rx2u,pngDir=pngDir,png2=pngName,ptw=ptw)
  plot( dw, e, ku1, "ae_ku1")
  plot(dwc,ec,ku1c,"aec_ku1")

def makeIndividual(dg1,dg2,dgS):
  for nrms in [0.0,1.0]:
    label = str(nrms)
    init(nrms)
    e1,u1c,dw1 = goPpPs1C(dg1)
    e2,u2c,dw2 = goPpPs2C(dg2)
    eS,uSg,dwS = goPs1Ps2(dgS)
    #plotVpVs(u1c,kvpvs1,dw1,sf,"PP","Vp/Vs 1")
    #plotVpVs(u2c,kvpvs2,dw2,sf,"PP","Vp/Vs 2")
    c1 = dw1.getCompression()
    c2 = dw2.getCompression()
    u1 = WarpUtils.getScaledShifts(sf,u1c,c1)
    u2 = WarpUtils.getScaledShifts(sf,u2c,c2)
    uS = WarpUtils.ps1ToPpTime(sf,u1,sg,uSg)
    ku1c = WarpUtils.getScaledShifts(sf,ku1,1.0/c1)
    ku2c = WarpUtils.getScaledShifts(sf,ku2,1.0/c2)
    kuSPP = WarpUtils.ps1ToPpTime(sf,ku1,sg,kuS)
    ptw = 504.0/2
    def plot(dw,e,u,ku,pngName):
      se = dw.getSampling1()
      su = dw.getSamplingU()
      el = [se.getFirst(),se.getLast()]
      ul = [su.getFirst(),su.getLast()]
      e = Transpose.transpose12(e)
      WarpUtils.normalizeErrors(e)
      ne = se.getCount()
      uA = [[copy(ne,u),"u","w-",2.0],[copy(ne,ku),"ku","w--",2.0]]
      plot2(e,pA=uA,s1=se,s2=su,hLabel="PP time (s)",vLabel="Time shift (s)",
            cbar=None,clips1=[0,0.3],width=750,height=400,hLimits=el,vLimits=ul,
            o=x1rx2u,png=pngDir+pngName,paint=paintGeo2,npng=2)
    #plot(dw1,e1, u1,ku1c,"ae_u1")
    #plot(dw2,e2, u2,ku2c,"ae_u2")
    #plot(dwS,eS,uSg, kuS,"ae_uS")
    def plotU(u,ku,c,pngName,w=1000,h=435,vLimits=[-0.15,0.2],vFormat=None):
      plot1([[u,c+"-"],[ku,"w--"]],s=sf,lineWidth=3.0,hLabel="PP time (s)",
            vLabel="Time shift (s)",vFormat=vFormat,vLimits=vLimits,width=w,
            height=h,o=x1rx2u,bgColor=bgColor,png=pngDir+pngName,
            paint=paintGeo2,npng=2)
    plotU(u1c, ku1c,"r","simw_u1_"+label,vFormat="%10f")
    #plotU( u1,  ku1,"r","simw_u1_"+label,vFormat="%10f",vLimits=None)
    plotU(u2c, ku2c,"b","simw_u2_"+label,vFormat="%10f")
    #plotU( u2,  ku2,"b","simw_u2_"+label,vFormat="%10f",vLimits=None)
    plotU( uS,kuSPP,"g","simw_uS_"+label,h=375,vLimits=[0.0,0.1])
    def plotGS(u1,u2,method,c,pngName):
      gs = method(sf,u1,u2)
      plot1([[gs,c+"-"],[kgs,"k--"]],s=sf,lineWidth=2.0,hLabel="PP time (s)",
            vLabel="GammaS",width=1000,height=435,vLimits=[-0.05,0.08],o=x1rx2u,
            png=pngDir+pngName,paint=paintGeo2,npng=2)
    #plotGS(u1,u2,WarpUtils.gammaSu12,"r","simw_gs12_"+label)
    #plotGS(u1,uS,WarpUtils.gammaSu1S,"b","simw_gs1S_"+label)
    #plotGS(u2,uS,WarpUtils.gammaSu2S,"g","simw_gs2S_"+label)
    gs12 = WarpUtils.gammaSu12(sf,u1,u2)
    gs1S = WarpUtils.gammaSu1S(sf,u1,uS)
    gs2S = WarpUtils.gammaSu1S(sf,u2,uS)
    #gsdA = [[sub(gs12,kgs),"m-"],[sub(gs1S,kgs),"y-."],[sub(gs2S,kgs),"c--."]]
    gsdA = [[sub(gs12,kgs),"m-"],[sub(gs1S,kgs),"y-"],[sub(gs2S,kgs),"c-"]]
    plot1(gsdA,s=sf,lineWidth=4.0,hLabel="PP time (s)",vLabel="Error",
          vLimits=[-0.1,0.1],vFormat="%10f",width=1000,height=435,o=x1rx2u,
          bgColor=bgColor,png=pngDir+"simw_gsr_"+label,paint=paintGeo2,npng=2)

def makeSimultaneous(dg1):
  for nrms in [0.0,1.0]:
    label = str(nrms)
    init(nrms)
    e,u1c,uSPPc,dw,sgc,gc,shc,hc = warp3Compress(dg1)
    c  = dw.getCompression()
    ci = 1.0/c
    u1 = WarpUtils.getScaledShifts(sf,u1c,c)
    uSPP = mul(c,uSPPc) # stretch
    #uSPP = WarpUtils.interpUS(sf,u1,uSPP)
    u2 = add(u1,uSPP)
    u2c  = WarpUtils.getScaledShifts(sf, u2,ci)
    ku1c = WarpUtils.getScaledShifts(sf,ku1,ci)
    ku2c = WarpUtils.getScaledShifts(sf,ku2,ci)
    kuSPP = WarpUtils.ps1ToPpTime(sf,ku1,sg,kuS)
    def plotU(u,ku,c,pngName,w=1000,h=435,vLimits=[-0.15,0.2],vFormat=None):
      plot1([[u,c+"-"],[ku,"w--"]],s=sf,lineWidth=2.0,hLabel="PP time (s)",
            vLabel="Time shift (s)",vFormat=vFormat,vLimits=vLimits,width=w,
            height=h,o=x1rx2u,bgColor=bgColor,pngDir=pngDir,png2=pngName,
            cwp=False)
    plotU( u1c, ku1c,"r","simw_u1Sim_"+label,vFormat="%10f")
    plotU( u2c, ku2c,"b","simw_u2Sim_"+label,vFormat="%10f")
    plotU(uSPP,kuSPP,"g","simw_uSSim_"+label,h=375,vLimits=[0.0,0.1])
    gs12 = WarpUtils.gammaSu12(sf,u1,u2)
    gs1S = WarpUtils.gammaSu1S(sf,u1,uSPP)
    gs2S = WarpUtils.gammaSu1S(sf,u2,uSPP)
    gsdA = [[sub(gs12,kgs),"m-"],[sub(gs1S,kgs),"y-"],[sub(gs2S,kgs),"c-"]]
    plot1(gsdA,s=sf,lineWidth=4.0,hLabel="PP time (s)",vLabel="Error",
          vLimits=[-0.1,0.1],vFormat="%10f",width=1000,height=435,o=x1rx2u,
          bgColor=bgColor,pngDir=pngDir,png2="simw_gsrSim_"+label,ptw=ptw,
          cwp=False)
    s1 = dw.getSampling1()
    n1 = s1.getCount()
    plotError3(e,dw,copy(n1,u1c),copy(n1,uSPPc),"simw_e3d_"+label)

def init(nrms):
  global sg,sh,f,g,h,ku1,ku2,kuS,kvpvs1,kvpvs2,kgs
  f,g,h,sg,sh,ku1,ku2,kuS,kvpvs1,kvpvs2,kgs = makeSequences(nrms)

def goPpPs1(dg):
  #dw = DynamicWarpingC.fromVpVs(sf,sg,vpvs1Avg,0.0)
  dw = DynamicWarpingC(0.0,1361*sf.getDelta(),Sampling(1950,0.002,0.0))
  rmin,rmax = (vpvs1Min-1.0)/2,(vpvs1Max-1.0)/2
  dw.setStrainLimits(rmin,rmax)
  e,u = warp(sf,f,sg,g,dw,dg)
  return e,u,dw

def goPpPs2(dg):
  # dw = DynamicWarpingC.fromVpVs(sf,sh,vpvs2Avg,0.0)
  dw = DynamicWarpingC(0.0,1361*sf.getDelta(),Sampling(1950,0.002,0.0))
  rmin,rmax = (vpvs2Min-1.0)/2,(vpvs2Max-1.0)/2
  dw.setStrainLimits(rmin,rmax)
  e,u = warp(sf,f,sh,h,dw,dg)
  return e,u,dw

def goPpPs1C(dg):
  # MAJOR HACK - need to have the same sampling grid as goPpPs2C, so we need
  # to force the same warping samplings. The second 'dw' uses values close
  # to those computed from the factory constructor 'fromVpVs' and the sampling
  # s1 from goPpPs2C
  dw = DynamicWarpingC.fromVpVs(sf,sg,vpvs1Avg-0.08,vpvs1Avg+0.1,vpvs1Avg,0.0)
  c = dw.getCompression()
  dw = DynamicWarpingC(-0.15,0.2,Sampling(1947,0.002,0.0))
  dw.setCompression(c)
  dw.setInterpolationMethod(Method.LINEAR)
  #dw.setInterpolationMethod(Method.MONOTONIC)
  gc = copy(g)
  WarpUtils.compress(c,gc)
  se = dw.getSampling1()
  sgc = Sampling(se.getCount(),se.getDelta(),se.getFirst())
  if dg==1:
    rmin,rmax = -1,1
  else:
    rmin,rmax = getSlopes(vpvs1Min,vpvs1Max,c)
  print "PP-PS1: rmin=%g, rmax=%g"%(rmin,rmax)
  dw.setStrainLimits(rmin,rmax)
  e,u = warp(sf,f,sgc,gc,dw,dg)
  return e,u,dw

def goPpPs2C(dg):
  # MAJOR HACK - need to have the same sampling grid as goPpPs1C, so we need
  # to force the same warping samplings. The second 'dw' uses values close
  # to those computed from the factory constructor 'fromVpVs'
  dw = DynamicWarpingC.fromVpVs(sf,sh,vpvs2Avg-0.1,vpvs2Avg+0.1,vpvs2Avg,0.0)
  c = dw.getCompression()
  dw = DynamicWarpingC(-0.15,0.2,Sampling(1947,0.002,0.0))
  dw.setCompression(c)
  dw.setInterpolationMethod(Method.LINEAR)
  #dw.setInterpolationMethod(Method.MONOTONIC)
  hc = copy(h)
  WarpUtils.compress(c,hc)
  se = dw.getSampling1()
  shc = Sampling(se.getCount(),se.getDelta(),se.getFirst())
  if dg==1:
    rmin,rmax = -1,1
  else:
    rmin,rmax = getSlopes(vpvs2Min,vpvs2Max,c)
  print "PP-PS2: rmin=%g, rmax=%g"%(rmin,rmax)
  dw.setStrainLimits(rmin,rmax)
  e,u = warp(sf,f,shc,hc,dw,dg)
  return e,u,dw

def goPs1Ps2(dg):
  dw = DynamicWarpingC(uSMin,uSMax,sg)
  if dg==1:
    dw.setStrainLimits(0,1)
  else:
    dw.setStrainLimits(rSMin,rSMax)
  dw.setInterpolationMethod(Method.LINEAR)
  #dw.setInterpolationMethod(Method.MONOTONIC)
  e,u = warp(sg,g,sh,h,dw,dg)
  return e,u,dw

def warp3Compress(dg1):
  dw = DynamicWarpingC3.fromVpVs(sf,sg,vpvs1Avg-0.1,vpvs1Avg+0.1,vpvs1Avg,
                                 0.0,uSMin,uSMax)
  c = dw.getCompression()
  dw.setWeight1(1.0)
  dw.setWeight2(1.0)
  dw.setWeight3(1.0)
  gc = copy(g)
  hc = copy(h)
  WarpUtils.compress(c,gc)
  WarpUtils.compress(c,hc)
  se = dw.getSampling1()
  sgc = Sampling(se.getCount(),se.getDelta(),se.getFirst())
  nhc = int(sh.getLast()/c/se.getDelta()+0.5)
  shc = Sampling(nhc,se.getDelta(),se.getFirst())
  print "se: first=%g, last=%g, n=%g, d=%g"%(se.getFirst(),se.getLast(),
                                             se.getCount(),se.getDelta())
  print "sgc: first=%g, last=%g, n=%g, d=%g"%(sgc.getFirst(),sgc.getLast(),
                                              sgc.getCount(),sgc.getDelta())
  print "shc: first=%g, last=%g, n=%g, d=%g"%(shc.getFirst(),shc.getLast(),
                                              shc.getCount(),shc.getDelta())
  gc = copy(sgc.getCount(),gc)
  hc = copy(shc.getCount(),hc)
  # plot1([[gc,"k-"]],sgc,"g")
  # plot1([[hc,"k-"]],shc,"h")
  rmin,rmax = getSlopes(vpvs1Min,vpvs1Max,c)
  e = computeErrors3(dw,sf,f,sgc,gc,shc,hc)
  e,u1,uS = warp3(dw,e,rmin,rmax,dg1)
  return e,u1,uS,dw,sgc,gc,shc,hc

def warp3(dw,e,rmin,rmax,dg1):
  dw.setStrainLimits(rmin,rmax,rSMin,rSMax)
  se = dw.getSampling1()
  g1 = Subsample.subsampleEnvelope(se.getCount(),f,dg1)
  g1 = Subsample.indicesToSampling(sf,g1)
  sw.restart()
  u1,uS = dw.findShifts(sf,e,g1)
  sw.stop()
  print "Find shifts time: %g seconds"%sw.time()
  return e,u1,uS
  
def warp(sf,f,sg,g,dw,dg):
  se = dw.getSampling1()
  g1 = Subsample.subsample(se.getCount(),dg)
  g1 = Subsample.indicesToSampling(sf,g1)
  #print "g1: "; dump(g1)
  e = dw.computeErrors(sf,f,sg,g)
  u = dw.findShifts(sf,e,g1)
  return e,u

def computeErrors3(dw,sf,f,sg,g,sh,h):
  sw.restart()
  e = dw.computeErrors(sf,f,sg,g,sh,h)
  sw.stop()
  print "Compute error time: %g seconds"%sw.time()
  return e

def getSlopes(vpvsMin,vpvsMax,c):
  rmin = (vpvsMin-2*c+1)/(2*c)
  rmax = (vpvsMax-2*c+1)/(2*c)
  return rmin,rmax

def makeSequences(nrms):
  seed = 1000
  vpvs1,u1,sg = get1()
  uS = getS(sg)
  uSPP = WarpUtils.ps1ToPpTime(sf,u1,sg,uS)
  u2 = add(u1,uSPP)
  vpvs2 = WarpUtils.vpvs(sf,u2)
  gammaSu12 = WarpUtils.gammaSu12(sf,u1,u2)
  gammaSu1S = WarpUtils.gammaSu1S(sf,u1,uSPP)
  gammaSu2S = WarpUtils.gammaSu2S(sf,u2,uSPP)
  if nf!=len(u2):
    print "nf!=len(u2)"
    sys.exit(0)
  nh = int((sf.getLast()+u2[nf-1])/df+0.5)
  sh = Sampling(nh,df,ff)
  print "nh=%g"%nh
  h = Synthetic.makeRandomEvents(nh,seed)
  g = VpVsWarp.warp(h,uS,sh,sg)
  f = VpVsWarp.warp(g,u1,sg,sf)
  h = Synthetic.addRickerWavelet(fpeak,h)
  g = Synthetic.addRickerWavelet(fpeak,g)
  f = Synthetic.addRickerWavelet(fpeak,f)
  h = Synthetic.addNoise(nrms,fpeak,h,seed*11+3)
  g = Synthetic.addNoise(nrms,fpeak,g,seed*11+2)
  f = Synthetic.addNoise(nrms,fpeak,f,seed*11+1)
#  if plot:
#    ulabel = "Time shit (s)"
#    flabel = "PP time (s)"
#    glabel = "PS1 time (s)"
#    hlabel = "PS2 time (s)"
#    def plot(fA,s,vl,hl,pngName,hLimits=None):
#      plot1(fA,s=s,vLabel=vl,hLabel=hl,hLimits=hLimits,width=400,height=900,
#            o=x1dx2r,pngDir=pngDir,png2=pngName)
#    plot([[h,"k-"]],sh,hlabel,"Amplitude","simw_h",hLimits=[-2.0,2.0])
#    plot([[g,"k-"]],sg,glabel,"Amplitude","simw_g",hLimits=[-2.0,2.0])
#    plot([[f,"k-"]],sf,flabel,"Amplitude","simw_f",hLimits=[-2.0,2.0])
#    plot([[vpvs1,"k-"],[vpvs2,"r-"]],sf,flabel,"Interval Vp/Vs","simw_kv")
#    plot([[gammaSu12,"k-"],[gammaSu1S,"r-"],[gammaSu2S,"b-"]],sf,flabel,
#          "GammaS","simw_kgs")
#    plot([[uS,"k-"]],sg,glabel,ulabel,"simw_kuS")
#    plot([[u1,"k-"],[u2,"b-"],[uSPP,"r-"]],sf,flabel,ulabel,"sim_ku",
#         hLimits=[0.0,3.0])
  return f,g,h,sg,sh,u1,u2,uS,vpvs1,vpvs2,gammaSu12

def get1():
  vpvs1 = VpVsWarp.getVpVs(m1,a1,b1,sf)
  vpvsAvg1 = VpVsWarp.getAverage(vpvs1)
  print "vpvsMin1=%g, vpvsMax1=%g, vpvsAvg1=%g"%(min(vpvs1),max(vpvs1),vpvsAvg1)
  ng = int(nf*(vpvsAvg1+1.0)*0.5)
  print "ng=%g"%ng
  sg = Sampling(ng,df,ff)
  u1 = VpVsWarp.getU(vpvs1,df)
  return vpvs1,u1,sg

def get2():
  vpvs2 = VpVsWarp.getVpVs(m2,a2,b2,sf)
  vpvsAvg2 = VpVsWarp.getAverage(vpvs2)
  print "vpvsMin2=%g, vpvsMax2=%g, vpvsAvg2=%g"%(min(vpvs2),max(vpvs2),vpvsAvg2)
  nh = int(nf*(vpvsAvg2+1.0)*0.5)
  print "nh=%g"%nh
  sh = Sampling(nh,df,ff)
  u2 = VpVsWarp.getU(vpvs2,df)

def getS(sg):
  uS = VpVsWarp.getU(aS,bS,sg)
  uS = mul(uS,VpVsWarp.getWeights(0.0,0.2,sg))
  return uS

def getShiftCoords(g1,u):
  ng1 = len(g1)
  x1 = zerofloat(ng1)
  x2 = zerofloat(ng1)
  for i1 in range(ng1):
    x1[i1] = g1[i1];
    x2[i1] = u[g1[i1]];
  return x1,x2

def normalize1(f,nMin,nMax):
  fMin = min(f)
  fMax = max(f)
  fr = fMax-fMin
  nr = nMax-nMin
  for i1 in range(len(f)):
    f[i1] = nr * (f[i1]-fMin) / fr + nMin

def toFloats1(f):
  n1 = len(f)
  ff = zerofloat(n1)
  for i1 in range(n1):
    ff[i1] = float(f[i1])
  return ff

################################################################################
# Plotting

def plotError2(e,se,su,uA,hname,desc):
  el = [se.getFirst(),se.getLast()]
  ul = [su.getFirst(),su.getLast()]
  e = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(e)
  plot2(e,pA=uA,title="AE "+desc,s1=se,s2=su,hLabel=hname+" time (s)",
        vLabel="Time shift (s)",cbar="Error",clips1=[0,0.2],width=900,
        height=600,hLimits=el,vLimits=ul,o=x1rx2u)

def plotWarped(pA,s,fname,desc):
  plot1(pA,s,title="Warp "+desc,vLabel=fname+" time (s)",width=400,height=900,
        vLimits=None,hLimits=[-1.5,1.5])

def plotVpVs(u,kvpvs,dw,s,vname,desc):
  vpvs = WarpUtils.vpvs(s,u,dw.getCompression())
  #vpvs = WarpUtils.vpvs(s,u,1.0)
  plot1([[vpvs,"k-"],[kvpvs,"r--"]],title=desc,s=s,vLabel=vname+" time (s)",
        hLabel="Vp/Vs C",width=400,height=900,hLimits=[vpvs1Min,vpvs1Max],
        o=x1dx2r)

def plotGammaS(u1,u2,s,method,vname,desc):
  #gs = WarpUtils.gammaS(s,u1,u2)
  gs = method(s,u1,u2)
  plot1([[gs,"k-"],[kgs,"r--"]],title=desc,s=s,vLabel=vname+" time (s)",
        hLabel="Gamma S",width=400,height=900,hLimits=[-0.05,0.08],o=x1dx2r)

def plotError3(e,dw,u1,uS,pngName):
  se = dw.getSampling1()
  su1 = dw.getSamplingU1()
  suS = dw.getSamplingUS()
  e = Transpose.transpose312(e)
  WarpUtils.normalizeErrors(e)
  xyz = getXYZ(se,u1,uS)
  #rgb = getRGB(se.getCount(),Color.YELLOW)
  rgb = getRGB(se.getCount(),Color.WHITE)
  plot3D(e,su1,suS,se,xyz=xyz,rgb=rgb,clips=[0.0,0.3],pngName=pngName)

def getXYZ(se,u1,uS):
  n1 = se.getCount()
  xyz = zerofloat(n1*3)
  j = 0
  for i1 in range(n1):
    #xyz[j  ] = uS[i1]
    #xyz[j+1] = u1[i1]
    #xyz[j+2] = se.getValue(i1)
    xyz[j  ] = se.getValue(i1)
    xyz[j+1] = uS[i1]
    xyz[j+2] = u1[i1]
    j += 3
  return xyz

def getRGB(n1,color):
  rgb = zerofloat(n1*3)
  r,g,b = color.getRed(),color.getGreen(),color.getBlue()
  j = 0
  for i in range(n1):
    rgb[j  ] = r
    rgb[j+1] = g
    rgb[j+2] = b
    j = j+3
  return rgb
    
def plot3D(f,s1,s2,s3,xyz=None,size=0.004,rgb=None,clips=None,cmap=None,
           width=800,height=1000,pngName=None):
  sf = SimpleFrame()
  ipg = sf.addImagePanels(s1,s2,s3,f)
  #ipg.setSlices(int(s1.getCount()*0.6),s2.getCount()-1,s3.getCount()-1)
  ipg.setSlices(int(s1.getCount()*0.6),s2.getCount()-1,890)
  if xyz:
    if rgb==None:
      rgb = getRGB(s1.getCount(),Color.WHITE)
    pg = PointGroup(size,xyz,rgb)
    world = sf.getWorld()
    world.addChild(pg)
  if clips:
    ipg.setClips(clips[0],clips[1])
  if cmap:
    ipg.setColorModel(cmap)
  sf.setSize(width,height)
  ov = sf.getOrbitView()
  ov.setScale(10.0)
  ov.setAxesScale(0.2,1.5,1.0)
  #ov.setAzimuth(45.0)
  ov.setElevation(230.0)
  if pngDir and pngName:
    sf.paintToFile(pngDir+pngName+".png")
  sf.setVisible(True)
 
##############################################################################
# Run the function main on the Swing thread
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
