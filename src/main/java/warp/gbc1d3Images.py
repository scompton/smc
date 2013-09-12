#############################################################################
# 2D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
s1pp,s2,s3 = getSamplings() # Samplings defined in gbcUtils
n1pp = s1pp.getCount() # n1pp is the fastest dimension
n2 = s2.getCount()
n3 = s3.getCount() # n3 is the slowest

# subset for PS data
n1ps = 1501
s1ps = Sampling(n1ps,s1pp.getDelta(),s1pp.getFirst())

# default clips for seismic images(i), time shifts(u), interval vp/vp(v)
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.5,2.1]
# iMap = bwr
iMap = gray
i2,i3 = 75,72 # central trace

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

#############################################################################

def main(args):
  goPpPs1()
  # goPpPs2()
  # goPs1Ps2()

def goPpPs1():
  pp,ps1 = getPp(i3),getPs1(i3)
  # plot2(ps1,title="PS1",s1=s1ps,s2=s2,clips1=iClips,vLabel="PS time (s)",
  #       vInterval=0.2,cmap1=iMap,cbw=100)
  lMax = DynamicWarpingC.computeMaxLag(n1ps,vpvsAvg)
  n1w = DynamicWarpingC.computeMaxLength(n1pp,n1ps,vpvsAvg)
  ppw = copy(n1w,n2,pp)
  dw = DynamicWarpingC(0,lMax) # DynamicWarpingC.java
  nl = dw.getNumberOfLags()
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(n1w,nl,n1pp,s1pp)
  # us = goShiftsSAG(pp,ppw,ps1,dw,50,13,"PP","PS1")
  uh = goShiftsHOR(pp,ppw,ps1,dw,ppHorizons,"PP","PS1")

def goPpPs2():
  pp,ps2 = getPp(i3),getPs2(i3)
  lMax = DynamicWarpingC.computeMaxLag(n1ps,vpvsAvg)
  n1w = DynamicWarpingC.computeMaxLength(n1pp,n1ps,vpvsAvg)
  ppw = copy(n1w,n2,pp)
  dw = DynamicWarpingC(0,lMax) # DynamicWarpingC.java
  nl = dw.getNumberOfLags()
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(n1w,nl,n1pp,s1pp)
  # us = goShiftsSAG(pp,ppw,ps2,dw,50,13,"PP","PS2")
  uh = goShiftsHOR(pp,ppw,ps2,dw,ppHorizons,"PP","PS2")
  
def goPs1Ps2():
  ps1,ps2 = getPs1(i3),getPs2(i3)
  plot2(ps1,title="PS1",s1=s1ps,s2=s2,clips1=iClips,vLabel="PS1 time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  plot2(ps2,title="PS2",s1=s1ps,s2=s2,clips1=iClips,vLabel="PS2 time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  lMax = DynamicWarpingC.computeMaxLag(n1ps,1.1)
  n1w = DynamicWarpingC.computeMaxLength(n1ps,n1ps,1.1)
  ps1w = copy(n1w,n2,ps1)
  dw = DynamicWarpingC(0,lMax) # DynamicWarpingC.java
  nl = dw.getNumberOfLags()
  dw.setStrainLimits(0.0,0.4,r2min,r2max)
  setPlotVars(n1w,nl,n1ps,s1ps)
  # us = goShiftsSAG(ps1,ps1w,ps2,dw,50,13,"PS1","PS2")
  uh = goShiftsHOR(ps1,ps1w,ps2,dw,ps1Horizons,"PS1","PS2")

def getPp(i3):
  pp = getGbcImage(baseDir,"pp_smooth")
  return pp[i3]

def getPs1(i3):
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  return ps1[i3]

def getPs2(i3):
  # makeSubset(baseDir,"ps2",0,n1ps-1,0,n2-1,0,n3-1)
  ps2 = getGbcImage(baseDir,"ps2_1501_150_145",1501)
  return ps2[i3]
  # makeSubset2(baseDir,"ps2_fk72",0,n1ps-1,0,n2-1)
  # return readImage(baseDir,"ps2_fk72_1501_150",1501,150)

def printConstraints(dg1):
  k1min = int(math.ceil( r1min*dg1))
  k1max = int(math.floor(r1max*dg1))
  k2min = int(math.ceil( r2min*dg2))
  k2max = int(math.floor(r2max*dg2))
  info = """Constraints:
  r1min=%g, r1max=%g, dg1=%g, k1min=%g, k1max=%g
  r2min=%g, r2max=%g, dg2=%g, k2min=%g, k2max=%g"""%(r1min,r1max,dg1,k1min,
                                                     k1max,r2min,r2max,dg2,
                                                     k2min,k2max)
  print info

def goShiftsSAG(f,fw,g,dw,dg1,ng,fname,gname):
  printConstraints(dg1)
  goSlopes(f) # compute slopes for flattening
  goFlat(f) # flatten
  x1m = readImage(twoDDir,"x1",ne1,n2) # flattening mappings
  ff = readImage(twoDDir,"ff",ne1,n2) # flattened image
  g1,g1Flat = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng) # structure aligned grid
  g2 = Subsample.subsample(n2,dg2) # simple regular interval grid
  print "g2:"; dump(g2)
  # Compute time shifts u with values only at g1 and g2 coordinates and
  # interpolate to get time shifts everywhere. Use Interp.MONOTONIC for smoother
  # shifts, but Interp.LINEAR is really the answer most consistent with the
  # method of finding time shifts. Interp.LINEAR for 2nd dimension should not be
  # changed.
  u = dw.findShifts(fw,g,getG1Floats(x1m,g1Flat),g2,Interp.MONOTONIC,
                    Interp.LINEAR)
  checkShifts2(u)
  g1f = zeroint(ne1,n2)
  for i2 in range(n2):
    g1f[i2] = copy(g1Flat)
  x1u,x2u = getShiftCoords(g1,g2,u)
  x1u,x2u = mul(x1u,d1),mul(x2u,d1)
  x12SliceA = [[x1u,x2u,"coarse grid","rO",10.0]]
  uA = [[mul(u,d1),"u","w-",2.0]]
  plotErrors(fw,g,dw,uA,fname,x12SliceA=x12SliceA)
  x1i,x2i = getSparseGridCoords(g1,g2); plotGrid(f,x1i,x2i,fname)
  plotShifts(f,u,fname)
  plotVpvs(f,u,fname)
  plotWarped(f,g,dw,u,fname,gname)
  return u

def getEnvelopeSum(ff):
  htf = HilbertTransformFilter()
  n2 = len(ff)
  n1 = len(ff[0])
  es = zerofloat(n1)
  for i2 in range(n2):
    x = ff[i2]
    y = copy(x)
    htf.apply(n1,x,y)
    for i1 in range(n1):
      es[i1] = es[i1] + sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return es

def makeAutomaticGrid(x1m,ff,dw,dg1,ng=None):
  es = getEnvelopeSum(ff)
  if ng:
    g1Flat = Subsample.subsample(es,dg1,ng)
  else:
    g1Flat = Subsample.subsample(es,dg1)
  ng = len(g1Flat)
  g1 = zeroint(ng,n2)
  for i2 in range(n2):
    for i1 in range(ng):
      x = x1m[i2][g1Flat[i1]]
      g1[i2][i1] = int(x+0.5)
  print "g1Flat:"; dump(g1Flat)
  return g1,g1Flat

def makeHorizonGrid(horizons):
  n1 = len(horizons)+2
  n1m = n1-1
  g1 = zerofloat(n1,n2)
  last = ne1-1
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

def getG1Floats(x1m,g1Flat):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2)
  for i2 in range(n2):
    for i1 in range(ng):
      g1[i2][i1] = x1m[i2][g1Flat[i1]]
  return g1

def getG1Ints(g1):
  ng = len(g1)
  g1i = zeroint(ng)
  for i1 in range(ng):
    g1i[i1] = int(g1[i1]+0.5)
  return g1i

def getG1Ints2(g1):
  ng2 = len(g1)
  ng1 = len(g1[0])
  g1i = zeroint(ng1,ng2)
  for i2 in range(ng2):
    for i1 in range(ng1):
      g1i[i2][i1] = int(g1[i2][i1]+0.5)
  return g1i

def getShiftCoords(g1,g2,u):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = fillfloat(-10,ng1,n2)
  x2 = fillfloat(-10,ng1,n2)
  for i2 in range(ng2):
    g2i = g2[i2]
    for i1 in range(ng1):
      x1[g2i][i1] = g1[g2i][i1]
      x2[g2i][i1] = u[g2i][g1[g2i][i1]]
  return x1,x2

def getShiftCoordsI2(g1,i2,u):
  ng1 = len(g1[0])
  x1 = fillfloat(-10,ng1)
  x2 = fillfloat(-10,ng1)
  for i1 in range(ng1):
    x1[i1] = g1[i2][i1]
    x2[i1] = u[int(g1[i2][i1]+0.5)]
  return x1,x2

def getShiftCoords1(g1,u):
  ng1 = len(g1)
  x1 = fillfloat(-10,ng1)
  x2 = fillfloat(-10,ng1)
  for i1 in range(ng1):
    x1[i1] = g1[i1]
    x2[i1] = u[int(g1[i1]+0.5)]
  return x1,x2

def getSparseGridCoords(g1,g2):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = zerofloat(ng1,ng2)
  x2 = zerofloat(ng1,ng2)
  for i2 in range(ng2):
    for i1 in range(ng1):
      x1[i2][i1] = g1[g2[i2]][i1];
      x2[i2][i1] = g2[i2];
  return x1,x2

#############################################################################

def setPlotVars(n1w,nl,nf1,sf1):
  global ne1,nel,se,su,el,ul,xl,n1,s1
  n1,s1 = nf1,sf1
  ne1,nel = n1w,nl
  se = Sampling(ne1,d1,f1) # error sampling
  su = Sampling(nel,d1,f1) # shift sampling
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  # el = [0.0,1.0] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting

def plotGrid(f,x1i,x2i,fname):
  # Plot f with grid points
  x1i = mul(x1i,d1)
  # x2i = add(25,x2i)
  x2i = mul(x2i,d2)
  plot2(f,x12=x1i,x22=x2i,title=fname,s1=s1,s2=s2,vLabel=fname+" time (s)",
        hLimits=xl,vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotShifts(f,u,fname):
  # Plot 2D shifts
  us = mul(u,d1)
  zm = ZeroMask(f)
  us = DynamicWarpingC.extrapolate(len(f[0]),us)
  zm.apply(0.00,us)
  plot2(us,title="Vertical shifts",s1=s1,s2=s2,vLabel=fname+" time (s)",
        cbar="Vertical shift (s)",cmap1=jet,clips1=uClips,vLimits=el,
        cbw=100)

def plotVpvs(f,u,fname):
  # Plot Vp/Vs ratios
  zm = ZeroMask(f)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(len(f[0]),vpvs)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plot2(f,g=vpvs,title="Vp/Vs",s1=s1,s2=s2,vLabel=fname+" time (s)",
        cbar="Vp/Vs ratio",clips1=iClips,cmap2=jet,clips2=vClips,hLimits=xl,
        vLimits=el,cmap1=iMap,cbw=100)

def plotWarped(f,g,dw,u,fname,gname):
  # Plot Warped g image and the difference between the warped and input.
  h = dw.applyShifts(n1,g,u)
  plot2(h,title=gname+" warped",s1=s1,s2=s2,vLabel=fname+" time (s)",
        vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)
  nrms = DynamicWarpingC.computeNrms(ne1,f,h)
  print "nrms=%g"%nrms
  # d = sub(h,f)
  # plot2(d,title=fname+"-"+gname+" warped"+s1=s1,s2=s2,
  #       vLabel=fname+" time (s)",vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotErrors(fw,g,dw,uA,fname,x12SliceA=None):
  # Plot Alignment Errors with coarse grid points and shifts
  e = dw.computeErrors2(fw,g)
  e = DynamicWarpingC.transposeLag12(e)
  DynamicWarpingC.normalizeErrors(e)
  plot3(e,pA=uA,x12SliceA=x12SliceA,title="AE",s1=se,s2=su,
        hLabel=fname+" time (s)",vLabel="Vertical shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1right_x2up)

###############################################################################
run(main)
