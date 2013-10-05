#############################################################################
# 2D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
twoDDir = baseDir+"2D/" # directory for writing intermediate files
horizonDir = "/data/seis/gbc/horizons/"
s1pp,s2,s3 = getSamplings() # Samplings defined in gbcUtils
n1pp = s1pp.getCount() # n1pp is the fastest dimension
n2 = s2.getCount()
n3 = s3.getCount() # n3 is the slowest

# subset for PS data
n1ps = 1501
s1ps = Sampling(n1ps,d1,f1)

# default clips for seismic images(i), time shifts(u), interval vp/vp(v)
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.5,2.1]
# iMap = bwr
iMap = gray
i3 = 72 # index of 2D slice to extract from 3D volume

# just a guess... This controls the maximum possible time shift and the maximum
# PP time sample that could correspond to the last PS time sample. A larger
# value increases the maximum possible shift and includes less PP samples for
# warping.
vpvsAvg = 1.8

# Contstraints for time shift slopes, which are physically related to interval
# Vp/Vs ratios. Compute slope parameters from vpvsMin/Max.
vpvsMin,vpvsMax = 1.5,2.5
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

r2min,r2max,dg2 = -0.10,0.10,25

ppHorizons = [
  "a_pp_time",
  "m_pp_time",
  "b_pp_time",
  "c_pp_time",
  "d_pp_time",
  "n_pp_time",
  "e_pp_time",
  # "f_pp_time",
  "l_pp_time",
  "g_pp_time",
  # "h_pp_time",
  "i_pp_time",
  "j_pp_time",
  # "k_pp_time"
]

ps1Horizons = [
  "a_ps1_time",
  "b_ps1_time",
  "c_ps1_time",
  "d_ps1_time",
  "e_ps1_time",
  "f_ps1_time",
  "g_ps1_time",
  "h_ps1_time",
  "i_ps1_time"
]
sigma = 8.0 # sigma for smoothing horizons

#############################################################################

def main(args):
  goPpPs1(True)
  goPpPs1(False)
  # goPpPs2(True)
  # goPpPs2(False)
  # goPs1Ps2()
  # go3Images()

def goPpPs1(tbt):
  fname = "PP"
  gname = "PS1"
  f,sf = getPp(i3)
  g,sg = getPs1(i3)
  if tbt:
    dw = DynamicWarpingC.fromVpVs(sf,sg,vpvsAvg,0.0,s2)
    dw.setStrainLimits(r1min,r1max,r2min,r2max)
    goTraceByTrace(sf,f,sg,g,dw,ppHorizons,fname,gname)
  else:
    goResidual(sf,f,sg,g,80,ppHorizons,fname,gname)

def goPpPs2(tbt):
  fname = "PP"
  gname = "PS2"
  f,sf = getPp(i3)
  g,sg = getPs2(i3)
  if tbt:
    dw = DynamicWarpingC.fromVpVs(sf,sg,vpvsAvg,0.0,s2)
    dw.setStrainLimits(r1min,r1max,r2min,r2max)
    goTraceByTrace(sf,f,sg,g,dw,ppHorizons,fname,gname)
  else:
    goResidual(sf,f,sg,g,150,ppHorizons,fname,gname)

def goPs1Ps2():
  fname = "PS1"
  gname = "PS2"
  f,sf = getPs1(i3)
  g,sg = getPs2(i3)
  dw = DynamicWarpingC(0.0,0.05,sg)
  dw.setStrainLimits(0.0,0.08,r2min,r2max)
  goTraceByTrace(sf,f,sg,g,dw,ppHorizons,fname,gname)

def go3Images():
  fname = "PP"
  gname = "PS1"
  hname = "PS2"
  f,sf = getPp(i3)
  g,sg = getPs1(i3)
  h,sh = getPs2(i3)
  dwa = DynamicWarpingC3.fromVpVs(sf,sg,vpvsAvg-0.1,vpvsAvg+0.1,vpvsAvg,
                                  0.0,0.0,0.03)
  dwa.setInterpolationMethod(Method.LINEAR)
  c = dwa.getCompression()
  rmin = WarpUtils.getSlope(vpvsMin,c)
  rmax = WarpUtils.getSlope(vpvsMax,c)
  dwa.setStrainLimits(rmin,rmax,0.0,0.2)
  g1 = getG1Envelope(sf,f,100,dwa)
  sw = Stopwatch()
  sw.restart()
  e = dwa.computeErrorsSum(sf,f,sg,g,sh,h)
  print "Computed errors in %g"%sw.time()
  sw.restart()
  ua1,uaS = dwa.findShifts(sf,sg,e,g1)
  print "Computed average shifts in %g"%sw.time()
  n1 = len(ua1)
  ua12 = zerofloat(n1,n2)
  uaS2 = zerofloat(n1,n2)
  for i2 in range(n2):
    ua12[i2] = copy(ua1)
  for i2 in range(n2):
    uaS2[i2] = copy(uaS)

  # Plot average results
  setLimits(dwa.getSampling1,dwa.getSamplingU1)
  sgc = getCompressedSampling(sg,c)
  shc = getCompressedSampling(sh,c)
  wa1 = WarpUtils.applyShifts(sf,ua12,sgc,g)
  waS = WarpUtils.applyShifts(sf,uaS2,shc,h)
  plotWarped(sf,f,wa1,fname,gname+" average")
  plotWarped(sf,f,waS,fname,hname+" average")
  plotVpvs(sf,f,ua12,c,fname)

def goTraceByTrace(sf,f,sg,g,dw,horizons,fname,gname):
  g1 = getG1Horizons(horizons,dw)
  g2 = getG2(dg2)
  u = getShifts(sf,f,sg,g,dw,g1,g2)
  checkShifts2(u)

  # Plot results
  setLimits(dw.getSampling1,dw.getSamplingU)
  w = WarpUtils.applyShifts(sf,u,sg,g)
  x1i,x2i = getSparseGridCoords(g1,s2,g2)
  x1u,x2u = getShiftCoords(se,g1,s2,g2,u)
  x12SliceA = [[x1u,x2u,"coarse grid","rO",10.0]]
  uA = [[copy(se.getCount(),n2,u),"u","w-",2.0]]
  plotErrors(sf,f,sg,g,dw,uA,fname,x12SliceA=x12SliceA)
  plotGrid(sf,f,x1i,x2i,fname)
  plotWarped(sf,f,w,fname,gname)
  plotVpvs(sf,f,u,1.0,fname)
  plotShifts(sf,f,u,fname)

def goResidual(sf,f,sg,g,dg1,horizons,fname,gname):
  # Initial compression of PS1, then find average shifts
  dwa = DynamicWarpingC.fromVpVs(sf,sg,vpvsAvg-0.1,vpvsAvg+0.1,vpvsAvg,0.0)
  dwa.setInterpolationMethod(Method.LINEAR)
  c = dwa.getCompression()
  rmin = WarpUtils.getSlope(vpvsMin,c)
  rmax = WarpUtils.getSlope(vpvsMax,c)
  dwa.setStrainLimits(rmin,rmax)
  # g1 = getG1Regular(sf,200,dwa)
  g1 = getG1Envelope(sf,f,dg1,dwa)
  ua,e = getShiftsAvg(sf,f,sg,g,dwa,g1)

  # Plot average results
  setLimits(dwa.getSampling1,dwa.getSamplingU)
  sgc = getCompressedSampling(sg,c)
  wa = WarpUtils.applyShifts(sf,ua,sgc,g)
  x1u,x2u = getShiftCoords1(se,g1,ua[0])
  x12SingleA = [[x1u,x2u,"coarse grid","rO",10.0]]
  uA = [[copy(se.getCount(),ua[0]),"ua","w-",2.0]]
  plotErrorsAvg(e,dwa,uA,fname,x12SingleA=x12SingleA)
  plotWarped(sf,f,wa,fname,gname+" average")
  plotVpvs(sf,f,ua,c,fname)

  # Find residual shifts
  dw = DynamicWarpingC(-0.05,0.05,dwa.getSampling1(),s2)
  dw.setInterpolationMethod(Method.LINEAR)
  r2min,r2max,dg2 = -0.10,0.10,25
  dw.setStrainLimits(rmin,rmax,r2min,r2max)
  g1 = getG1Horizons(horizons,dw)
  g2 = getG2(dg2)
  ur = getShifts(sf,f,sf,wa,dw,g1,g2)

  # Plot residual results
  setLimits(dw.getSampling1,dw.getSamplingU)
  w = WarpUtils.applyShifts(sf,ur,sf,wa)
  x1i,x2i = getSparseGridCoords(g1,s2,g2)
  x1u,x2u = getShiftCoords(se,g1,s2,g2,ur)
  x12SliceA = [[x1u,x2u,"coarse grid","rO",10.0]]
  uA = [[copy(se.getCount(),n2,ur),"ur","w-",2.0]]
  plotErrors(sf,f,sf,wa,dw,uA,fname,x12SliceA=x12SliceA)
  plotGrid(sf,f,x1i,x2i,fname)
  plotWarped(sf,f,w,fname,gname)
  uc = WarpUtils.compositeShifts(sf,ua,ur)
  plotVpvs(sf,f,uc,c,fname)
  fsA = [[ua[72],"r-"],[ur[72],"b-"],[uc[72],"k-"]]
  plot1(fsA,sf,title="Shifts",hLabel="Time shift (s)",vLabel=fname+" time (s)",
        width=400,height=800,o=x1down_x2right)

def getPp(i3):
  pp = getGbcImage(baseDir,"pp_smooth")
  spp = Sampling(len(pp[0][0]),d1,f1)
  plot2(pp[i3],title="PP",s1=spp,s2=s2,clips1=iClips,vLabel="PP time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  return pp[i3],spp

def getPs1(i3):
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  # ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth")
  sps1 = Sampling(len(ps1[0][0]),d1,f1)
  # plot2(ps1[i3],title="PS1",s1=sps1,s2=s2,clips1=iClips,vLabel="PS1 time (s)",
  #       vInterval=0.2,cmap1=iMap,cbw=100)
  return ps1[i3],sps1

def getPs2(i3):
  # makeSubset(baseDir,"ps2",0,n1ps-1,0,n2-1,0,n3-1)
  # ps2 = getGbcImage(baseDir,"ps2_1501_150_145",1501)
  ps2 = getGbcImage(baseDir,"ps2_fk_smooth")
  sps2 = Sampling(len(ps2[0][0]),d1,f1)
  plot2(ps2[i3],title="PS2",s1=sps2,s2=s2,clips1=iClips,vLabel="PS2 time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  return ps2[i3],sps2
  # makeSubset2(baseDir,"ps2_fk72",0,n1ps-1,0,n2-1)
  # return readImage(baseDir,"ps2_fk72_1501_150",1501,150)

def getShifts(sf,f,sg,g,dw,g1,g2):
  return dw.findShifts(sf,f,sg,g,g1,g2)

def getShiftsAvg(sf,f,sg,g,dw,g1):
  e = dw.computeErrorsSum(sf,f,sg,g)
  u = dw.findShifts(sf,e,g1)
  n1 = len(u)
  u2 = zerofloat(n1,n2)
  for i2 in range(n2):
    u2[i2] = copy(u)
  return u2,e

def getG1Regular(sf,dg1,dw):
  se = dw.getSampling1()
  g1 = Subsample.subsample(se.getCount(),dg1)
  return Subsample.indicesToSampling(sf,g1)

def getG1Envelope(sf,f,dg1,dw):
  se = dw.getSampling1()
  g1 = Subsample.subsampleEnvelope(se.getCount(),f,dg1)
  return Subsample.indicesToSampling(sf,g1)

def getG1Horizons(horizons,dw):
  se = dw.getSampling1()
  g1 = makeHorizonGrid(se.getCount(),horizons) # g1 samples
  return mul(d1,g1)

def getG2(dg2):
  g2 = Subsample.subsample(n2,dg2) # simple regular interval grid
  g2 = Subsample.indicesToSampling(s2,g2)
  print "g2:"; dump(g2)
  return g2

def makeHorizonGrid(ne,horizons):
  n1 = len(horizons)+2
  n1m = n1-1
  g1 = zerofloat(n1,n2)
  last = ne-1
  for i2 in range(n2):
    g1[i2][0  ] = 0.0
    g1[i2][n1m] = last
  ref = RecursiveExponentialFilter(sigma)
  z = zerofloat(n2)
  d1i = 1.0/d1; # convert from seconds to samples
  for i1 in range(1,n1-1):
    h = horizonRead(horizons[i1-1])
    xyz = h.getI3X321(s2,s3,i3)
    for i2 in range(n2): # get horizon row for smoothing
      z[i2] = xyz[3*i2+2]*d1i
    ref.apply(z,z)
    for i2 in range(n2):
      g1[i2][i1] = z[i2]
  return g1

def getShiftCoords(s1,g1,s2,g2,u):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = fillfloat(-10,ng1,n2)
  x2 = fillfloat(-10,ng1,n2)
  for i2 in range(ng2):
    g2i = s2.indexOfNearest(g2[i2])
    for i1 in range(ng1):
      g1i = s1.indexOfNearest(g1[g2i][i1])
      x1[g2i][i1] = g1[g2i][i1]
      x2[g2i][i1] = u[g2i][g1i]
  return x1,x2

def getShiftCoords1(s1,g1,u):
  ng1 = len(g1)
  x1 = fillfloat(-10,ng1)
  x2 = fillfloat(-10,ng1)
  for i1 in range(ng1):
    g1i = s1.indexOfNearest(g1[i1])
    x1[i1] = g1[i1]
    x2[i1] = u[g1i]
  return x1,x2

def getSparseGridCoords(g1,s2,g2):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = zerofloat(ng1,ng2)
  x2 = zerofloat(ng1,ng2)
  for i2 in range(ng2):
    g2i = s2.indexOfNearest(g2[i2])
    for i1 in range(ng1):
      x1[i2][i1] = g1[g2i][i1];
      x2[i2][i1] = g2[i2];
  return x1,x2

def horizonFile(name):
  return horizonDir+name+".dat"

def horizonRead(name):
  return Horizon.readBinary(horizonFile(name))

def getCompressedSampling(sg,c):
  fg = sg.getFirst()
  lg = sg.getLast()
  dg = sg.getDelta()
  return Sampling(int(lg/c/dg)+1,dg,fg)

#############################################################################

# def setLimits(dw):
#   global se,el,su,ul,xl
#   se = dw.getSampling1()
#   su = dw.getSamplingU()
#   el = [se.getFirst(),se.getLast()]
#   ul = [su.getFirst(),su.getLast()]
#   xl = [s2.getFirst(),s2.getLast()]

def setLimits(mse,msu):
  global se,el,su,ul,xl
  se = mse()
  su = msu()
  el = [se.getFirst(),se.getLast()]
  ul = [su.getFirst(),su.getLast()]
  xl = [s2.getFirst(),s2.getLast()]

def plotGrid(sf,f,x1i,x2i,fname):
  # Plot f with grid points
  plot2(f,x12=x1i,x22=x2i,title=fname,s1=sf,s2=s2,vLabel=fname+" time (s)",
        hLimits=xl,vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotShifts(sf,f,u,fname):
  # Plot 2D shifts
  zm = ZeroMask(f)
  uc = copy(u)
  zm.apply(0.00,uc)
  plot2(uc,title="Time shifts",s1=sf,s2=s2,vLabel=fname+" time (s)",
        cbar="Time shift (s)",vLimits=el,cmap1=jet,clips1=uClips,cbw=100)

def plotVpvs(sf,f,u,c,fname):
  # Plot Vp/Vs ratios
  zm = ZeroMask(f)
  vpvs = WarpUtils.vpvs(sf,u,c)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plot2(f,g=vpvs,title="Vp/Vs",s1=s1,s2=s2,vLabel=fname+" time (s)",
        cbar="Vp/Vs ratio",clips1=iClips,cmap2=jet,clips2=vClips,hLimits=xl,
        vLimits=el,cmap1=iMap,cbw=100)

def plotWarped(sf,f,w,fname,gname):
  # Plot Warped g image and the difference between the warped and input.
  plot2(w,title=gname+" warped",s1=sf,s2=s2,vLabel=fname+" time (s)",
        vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)
  print "NRMS 1: %g"%WarpUtils.computeNrms(se.getCount(),f,w)
  # d = sub(h,f)
  # plot2(d,title=fname+"-"+gname+" warped"+s1=s1,s2=s2,
  #       vLabel=fname+" time (s)",vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotErrors(sf,f,sg,g,dw,uA,fname,x12SliceA=None):
  # Plot Alignment Errors with coarse grid points and shifts
  e = dw.computeErrors2(sf,f,sg,g)
  e = Transpose.transpose312(e)
  WarpUtils.normalizeErrors(e)
  plot3(e,pA=uA,x12SliceA=x12SliceA,title="AE",s1=se,s2=su,
        hLabel=fname+" time (s)",vLabel="Time shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1right_x2up)

def plotErrorsAvg(e,dw,uA,fname,x12SingleA=None):
  # Plot Alignment Errors with coarse grid points and shifts
  e = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(e)
  plot2(e,pA=uA,title="AE",s1=se,s2=su,x12SingleA=x12SingleA,
        hLabel=fname+" time (s)",vLabel="Time shift (s)",cbar="Error",
        clips1=[0,0.5],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1right_x2up)

###############################################################################
run(main)
