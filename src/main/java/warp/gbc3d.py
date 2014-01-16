#############################################################################
# 3D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
threeDDir = baseDir+"3D/" # directory for writing intermediate files
s1,s2,s3 = getSamplings() # Samplings defined in gbcUtils
n1 = s1.getCount() # n1 is the fastest dimension
n2 = s2.getCount()
n3 = s3.getCount() # n3 is the slowest

# subset for PS data
n1ps = 1501
s1ps = Sampling(n1ps,s1.getDelta(),s1.getFirst())

# default clips for seismic images(i), time shifts(u), interval vp/vp(v)
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.5,2.2]
# iMap = bwr
iMap = gray
k1,k2,k3 = 345,83,48 # default indices for plotting
slices = [k1,k2,k3]
he0 = 332 # height of time slice panel in plots

# just a guess... This controls the maximum possible time shift and the maximum
# PP time sample that could correspond to the last PS time sample. A larger
# value increases the maximum possible shift and includes less PP samples for
# warping.
vpvsAvg = 1.9

# Contstraints for time shift slopes, which are physically related to interval
# Vp/Vs ratios. Compute slope parameters from vpvsMin/Max.
vpvsMin,vpvsMax = 1.5,2.3
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

# Constraints for horizontal directions. A better match between the PS and PP
# images can be achieved by loosening these constraints, but I'm not sure that
# we want to correct for all lateral variations seen in the PS image but not
# in the PP image. These differences may not be related to Vp/Vs ratios.
r2min,r2max,dg2 = -0.15,0.15,25
r3min,r3max,dg3 = -0.15,0.15,25

#############################################################################

def main(args):
  goSagVsReg()
  # goSagFixedU()

def goSagVsReg():
  pp,ps1,dw = getPpPs1Data()
  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  # Fasle to read shifts from disk, True to find shifts.
  us = goShiftsSAG(s1,pp,s1ps,ps1,dw,50,13,False)

  # Warping using a regular grid.
  # dg1=80 - use a grid with layers spaced at a minimum of 80 samples apart,
  # including the first and last time samples. This is not optimal, but is
  # simple. The structure aligned grid should be more accurate.
  # Fasle to read shifts from disk, True to find shifts.
  ur = goShiftsREG(s1,pp,s1ps,ps1,dw,80,False)

  # uA = [[mul(us[k3],d1),"us","w-",2.0],[mul(ur[k3],d1),"ur","w--",2.0]]
  # plotErrors(ppw,ps1,dw,uA," REG vs SAG") # Compare time shifts on error plot

def goSagFixedU():
  pp,ppw,ps1,dw = getPpPs1Data()
  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  # Fasle to read shifts from disk, True to find shifts.
  us = goShiftsSAGFixedU(pp,ppw,ps1,dw,50,13,True)

def getPpPs1Data():
  # pp,ps1 = getImages3D() # FKK and AGC
  pp,ps1 = getImages3DSmooth() # FKK, AGC, and Structure Oriented Smoothing
  dw = DynamicWarpingC.fromVpVs(s1,s1ps,vpvsAvg,0.0,s2,s3)
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(dw)
  return pp,ps1,dw

def getImages3D():
  pp = getGbcImage(baseDir,"pp")
  # makeSubset(baseDir,"ps1_fkk",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_1501_150_145",1501)
  plotPP3(ps1,title="PS1",s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=slices,
          label1="PS time (s)",vInterval1=0.2,cmap1=iMap,cbw=100,he0=he0)
  return pp,ps1

def getImages3DSmooth():
  pp = getGbcImage(baseDir,"pp_smooth")
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  plotPP3(ps1,title="PS1",s1=s1ps,s2=s2,s3=s3,clips1=iClips,slices=slices,
          label1="PS time (s)",vInterval1=0.2,cmap1=iMap,cbw=100,he0=he0)
  return pp,ps1

def printConstraints(desc,dg1):
  k1min = int(math.ceil( r1min*dg1))
  k1max = int(math.floor(r1max*dg1))
  k2min = int(math.ceil( r2min*dg2))
  k2max = int(math.floor(r2max*dg2))
  k3min = int(math.ceil( r3min*dg3))
  k3max = int(math.floor(r3max*dg3))
  info = """Constraints: %s
  r1min=%g, r1max=%g, dg1=%g, k1min=%g, k1max=%g
  r2min=%g, r2max=%g, dg2=%g, k2min=%g, k2max=%g
  r3min=%g, r3max=%g, dg3=%g, k3min=%g, k3max=%g"""%(desc,r1min,r1max,dg1,k1min,
                                                     k1max,r2min,r2max,dg2,
                                                     k2min,k2max,r3min,r3max,
                                                     dg3,k3min,k3max)
  print info

def goShiftsREG(sf,f,sg,g,dw,dg1,find):
  desc = " regular grid"
  label = "reg"
  printConstraints(desc,dg1)
  g11 = Subsample.subsample(ne1,dg1) # simple regular interval grid
  ng1 = len(g11)
  g11f = zerofloat(ng1)
  for i1 in range(ng1):
    g11f[i1] = se.getValue(g11[i1])
  g1 = zerofloat(ne1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      g1[i3][i2] = copy(g11f) # use the same grid g11 for all traces
  g2 = getG2(dg2)
  g3 = getG3(dg3)
  print "g1:"; dump(g1[0][0])
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  return goShifts(sf,f,sg,g,g1,g2,g3,dw,desc,label,find)

def goShiftsSAG(sf,f,sg,g,dw,dg1,ng,find):
  desc = " structure aligned grid"
  label = "sag"
  printConstraints(desc,dg1)
  # Can do flattening separately, comment these steps out, and read the results
  # from disk.
  # goSlopes(f,dw) # compute slopes for flattening
  # goFlat(f,dw) # flatten
  x1m = readImage(threeDDir,"x1",ne1,n2,n3) # flattening mappings
  ff = readImage(threeDDir,"ff",ne1,n2,n3) # flattened image
  g1,g1Flat = makeAutomaticGrid(x1m,ff,dw,dg1,ng) # structure aligned grid
  g1 = getG1Floats(x1m,g1Flat,sf.getDelta())
  g2 = getG2(dg2)
  g3 = getG3(dg3)
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  return goShifts(sf,f,sg,g,g1,g2,g3,dw,desc,label,find)

def goShifts(sf,f,sg,g,g1,g2,g3,dw,desc,label,find):
  if find:
    u = dw.findShifts(sf,f,sg,g,g1,g2,g3)
    writeImage(threeDDir,"u_"+label,u)
  u = readImage(threeDDir,"u_"+label,n1,n2,n3)
  checkShifts3(u)
  coordMap = Viewer3P.getSparseCoordsMap(s1,g1,s2,g2,s3,g3)
  plotGrid(f,coordMap,desc)
  plotWarped(sf,f,sg,g,u,desc)
  plotVpvs(f,u,desc)
  plotShifts(f,u,desc)
  return u

def goShiftsSAGFixedU(pp,ppw,ps1,dw,dg1,ng,find):
  desc = " structure aligned grid"
  label = "sag"
  printConstraints(desc,dg1)
  # Can do flattening separately, comment these steps out, and read the results
  # from disk.
  # goSlopes(pp,dw) # compute slopes for flattening
  # goFlat(pp,dw) # flatten
  x1m = readImage(threeDDir,"x1",ne1,n2,n3) # flattening mappings
  ff = readImage(threeDDir,"ff",ne1,n2,n3) # flattened image
  g1,g1Flat = makeAutomaticGrid(x1m,ff,dw,dg1,ng) # structure aligned grid
  g2 = Subsample.subsample(n2,dg2) # simple regular interval grid
  g3 = Subsample.subsample(n3,dg3) # simple regular interval grid
  # Set known shifts
  ksu = KnownShiftUtil(su,se,s2,s3)
  ksu.setStrainMin(r1min)
  # p1s1tA = [[0.68,0.937,0.833,2.414],[0.706,0.955,3.613,2.414],
  #           [1.916,2.616,2.399,2.414]]
  # p1s1tA = [[0.68,0.850,0.833,2.414]] # shift violates Vp/Vs constraints
  p1s1tA = [[0.68,1.037,0.833,2.414]] # incorrect shift
  nks = len(p1s1tA)
  for i in range(nks):
    p1s1t = p1s1tA[i]
    s = p1s1t[1]-p1s1t[0]
    ksu.add(s,p1s1t[0],p1s1t[2],p1s1t[3])
  map = ksu.getMap()
  # dw.setStrainLimits(r1min,r1max,-0.2,0.2)
  g2ks,g3ks = ksu.getG23(g2,g3,filldouble(r2min,n2),filldouble(r2max,n2),
                         filldouble(r3min,n3),filldouble(r3max,n3))
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  print "g2ks:",dump(g2ks)
  print "g3ks:",dump(g3ks)
  if find:
    # Computes time shifts u with values only at g1, g2, and g3 coordinates and
    # interpolate to get time shifts everywhere. Use Interp.MONOTONIC for
    # smoother shifts, but Interp.LINEAR is really the answer most consistent
    # with the method of finding time shifts. With MONOTONIC the min/max Vp/Vs
    # constraints may be violated. Interp.LINEAR for 2nd/3rd dimension should
    # not be changed.
    g1f = getG1Floats(x1m,g1Flat)
    # u = dw.findShifts(ppw,ps1,g1f,g2,g3,Interp.MONOTONIC,Interp.LINEAR)
    # writeImage(threeDDir,"u_"+label,u)

    g1ks = g1f
    uks = dw.findShifts(ppw,ps1,g1ks,g2ks,g3ks,Interp.MONOTONIC,Interp.LINEAR,
                        map)
    writeImage(threeDDir,"uks_"+label,uks)
  u = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  uks = readImage(threeDDir,"uks_"+label,ne1,n2,n3)
  for i in range(nks):
    p1s1t = p1s1tA[i]
    s = p1s1t[1]-p1s1t[0]
    i3 = s3.indexOfNearest(p1s1t[3])
    i2 = s2.indexOfNearest(p1s1t[2])
    i1 = se.indexOfNearest(p1s1t[0])
    csu,csuks = u[i3][i2][i1]*d1,uks[i3][i2][i1]*d1
    print "i3=%d, i2=%d, i1=%d: s=%g, csu=%g, csuks=%g"%(i3,i2,i1,s,csu,csuks)
  checkShifts3(u)
  checkShifts3(uks)
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  plotGrid(pp,coordMap,desc)
  plotShifts(pp,u,desc)
  plotVpvs(pp,u,desc)
  plotWarped(pp,ps1,dw,u,desc)
  uA = [[mul(u[k3],d1),"u","w-",2.0],[mul(uks[k3],d1),"uks","w--",2.0]]
  plotErrors(ppw,ps1,dw,uA," SAG vs SAG KS") # Compare time shifts on error plot

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

def makeAutomaticGrid(x1m,ff,dw,dg1,ng):
  es = getEnvelopeSum(ff)
  g1Flat = Subsample.subsample(es,dg1,ng)
  ng = len(g1Flat)
  g1 = zeroint(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        x = x1m[i3][i2][g1Flat[i1]]
        g1[i3][i2][i1] = int(x+0.5)
  print "g1Flat:"; dump(g1Flat)
  return g1,g1Flat

def getG1Floats(x1m,g1Flat,d1):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        g1[i3][i2][i1] = x1m[i3][i2][g1Flat[i1]]*d1
  return g1

def getG2(dg2):
  g2 = Subsample.subsample(n2,dg2); # simple regular interval grid
  ng2 = len(g2)  
  g2f = zerofloat(ng2)
  for i2 in range(ng2):
    g2f[i2] = s2.getValue(g2[i2])
  return g2f

def getG3(dg3):
  g3 = Subsample.subsample(n3,dg3); # simple regular interval grid
  ng3 = len(g3)  
  g3f = zerofloat(ng3)
  for i3 in range(ng3):
    g3f[i3] = s3.getValue(g3[i3])
  return g3f

def toFloats(f):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  ff = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        ff[i3][i2][i1] = float(f[i3][i2][i1])
  return ff
  
def getShiftCoords(g1,g2,u):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = fillfloat(-10,ng1,n2)
  x2 = fillfloat(-10,ng1,n2)
  for i2 in range(ng2):
    g2i = g2[i2]
    for i1 in range(ng1):
      x1[g2i][i1] = g1[g2i][i1];
      x2[g2i][i1] = u[g2i][g1[g2i][i1]];
  return x1,x2

def goSlopes(f,dw):
  fs = copy(ne1,n2,n3,f) # f with only ne1 samples
  zm = ZeroMask(0.1,1.0,1.0,1.0,fs)
  sigma1 = 8.0
  sigma2 = 8.0
  sigma3 = 8.0
  print "LSF: sigma1=%g, sigma2=%g, sigma3=%g"%(sigma1,sigma2,sigma3)
  pmax = 2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma3,pmax)
  p2 = zerofloat(ne1,n2,n3)
  p3 = zerofloat(ne1,n2,n3)
  ep = zerofloat(ne1,n2,n3)
  lsf.findSlopes(fs,p2,p3,ep)
  zero = 0.00;
  tiny = 0.01;
  zm.apply(zero,p2);
  zm.apply(zero,p3);
  zm.apply(tiny,ep);
  writeImage(threeDDir,"p2",p2)
  writeImage(threeDDir,"p3",p3)
  writeImage(threeDDir,"ep",ep)
  s1f = Sampling(ne1,d1,f1)
  plotPP3(p2,title="P2",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",cbar="Slope",
          cmap1=bwr,slices=slices,he0=he0)
  plotPP3(p3,title="P3",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",cbar="Slope",
        cmap1=bwr,slices=slices,he0=he0)
  plotPP3(ep,title="Planarity",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",
        cbar="Planarity",cmap1=jet,clips1=[0.0,1.0],slices=slices,he0=he0)

def goFlat(f,dw):
  s1f = Sampling(ne1)
  s2f = Sampling(n2)
  s3f = Sampling(n3)
  p2 = readImage(threeDDir,"p2",ne1,n2,n3)
  p3 = readImage(threeDDir,"p3",ne1,n2,n3)
  ep = readImage(threeDDir,"ep",ne1,n2,n3)
  ep = pow(ep,2)
  fl = Flattener3()
  fl.setWeight1(1.0)
  fl.setIterations(0.1,1000)
  fm = fl.getMappingsFromSlopes(s1f,s2f,s3f,p2,p3,ep)
  ff = fm.flatten(f)
  h = fm.unflatten(ff)
  s = fm.getShiftsS()
  x1 = fm.x1
  y1 = fm.u1
  writeImage(threeDDir,"x1",x1)
  writeImage(threeDDir,"y1",y1)
  writeImage(threeDDir,"ff",ff)
  writeImage(threeDDir,"fs",s)
  sfs = Sampling(ne1,d1,f1)
  plotPP3(ff,title="PP flat",s2=s2,s3=s3,label1="Tau",slices=slices,he0=he0)
  plotPP3(h,title="PP unflattened",s1=sfs,s2=s2,s3=s3,label1="PP time (s)",
          slices=slices,he0=he0)
  plotPP3(s,title="PP flattening shifts",s1=sfs,s2=s2,s3=s3,slices=slices,
          label1="PP time (s)",cbar="Shift (samples)",cmap1=jet,clips1=None,
          he0=he0)
  print "average shift =",sum(s)/(n1*n2),"samples"

#############################################################################

def setPlotVars(dw):
  global ne1,nel,se,su,el,ul,xl,yl
  se = dw.getSampling1() # error sampling
  su = dw.getSamplingU() # shift sampling
  ne1 = se.getCount()
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  yl = [s3.getFirst(),s3.getLast()] # frame limits for plotting
  
def plotGrid(pp,coordMap,desc):
  # Plot PP with grid points
  s1 = Sampling(len(pp[0][0]),d1,f1)
  cm = [coordMap,"rO",6.0]
  plotPP3(pp,coordMap=cm,title="PP",s1=s1,s2=s2,s3=s3,clips1=iClips,
          label1="PP time (s)",vInterval1=0.2,cmap1=iMap,cbw=100,slices=slices,
          limits1=el,limits2=xl,limits3=yl,he0=he0)

def plotShifts(f,u,desc):
  # Plot 2D shifts
  zm = ZeroMask(f)
  zm.apply(0.00,copy(u))
  plotPP3(u,title="Vertical shifts"+desc,s1=s1,s2=s2,s3=s3,clips1=uClips,
          label1="PS time (s)",cbar="Vertical shift (s)",cmap1=jet,cbw=100,
          limits1=el,slices=slices,he0=he0)

def plotVpvs(f,u,desc,coordMap=None):
  # Plot Vp/Vs ratios
  zm = ZeroMask(f)
  vpvs = WarpUtils.vpvs(su,u)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  cm = None
  if coordMap:
    cm = [coordMap,"rO",6.0]
  plotPP3(vpvs,coordMap=cm,title="Vp/Vs"+desc,s1=s1,s2=s2,s3=s3,
          clips1=vClips,label1="PP time (s)",cbar="Vp/Vs ratio",cmap1=jet,
          cbw=100,limits1=el,limits2=xl,limits3=yl,slices=slices,he0=he0)

def plotWarped(sf,f,sg,g,u,desc):
  # Plot Warped PS1 image and the difference between the warped and input.
  h = WarpUtils.applyShifts(sf,u,sg,g)
  plotPP3(h,title="PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
          cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices,
          he0=he0)
  nrms = WarpUtils.computeNrms(ne1,f,h)
  print desc+": nrms=%g"%nrms
  # d = sub(h,f)
  # plotPP3(d,title="PP-PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
  #         cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices)

def plotErrors(ppw,ps1,dw,uA,desc,x12SliceA=None):
  # Plot Alignment Errors with coarse grid points and shifts
  # u = mul(u[k3],d1)
  # p = [u,"shifts","w-",2.0]
  # p2,p3 = None,None
  # if u2:
  #   u2 = mul(u2[k3],d1)
  #   p2 = [u2,"shifts 2","c--",2.0]
  # if u3:
  #   u3 = mul(u3[k3],d1)
  #   p3 = [u3,"shifts 3","r-.",2.0]
  # if x1u and x2u: 
  #   x1u = mul(x1u,d1)
  #   x2u = mul(x2u,d1)
  e = dw.computeErrors2(ppw[k3],ps1[k3])
  e = DynamicWarpingC.transposeLag12(e)
  DynamicWarpingC.normalizeErrors(e)
  plot3(e,pA=uA,x12SliceA=x12SliceA,title="AE"+desc,s1=se,s2=su,
        hLabel="PP time (s)",vLabel="Vertical shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1rx2u)

###############################################################################
run(main)
