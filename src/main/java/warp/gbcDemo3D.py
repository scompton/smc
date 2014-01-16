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
vpvsMin,vpvsMax = 1.5,2.5
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

# Constraints for horizontal directions. A better match between the PS and PP
# images can be achieved by loosening these constraints, but I'm not sure that
# we want to correct for all lateral variations seen in the PS image but not
# in the PP image. These differences may not be related to Vp/Vs ratios.
r2min,r2max,dg2 = -0.15,0.15,15
r3min,r3max,dg3 = -0.15,0.15,15

#############################################################################

def main(args):
  # pp,ps1 = getImages3D() # FKK and AGC
  pp,ps1 = getImages3DSmooth() # FKK, AGC, and Structure Oriented Smoothing
  lMax = DynamicWarpingC.computeMaxLag(n1ps,vpvsAvg)
  n1w = DynamicWarpingC.computeMaxLength(n1,n1ps,vpvsAvg)
  ppw = copy(n1w,n2,n3,pp)
  dw = DynamicWarpingC(0,lMax) # DynamicWarpingC.java
  nl = dw.getNumberOfLags()
  dw.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max)
  setPlotVars(n1w,nl)

  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  # Fasle to read shifts from disk, True to find shifts.
  us = goShiftsSAG(pp,ppw,ps1,dw,50,13,False)

  # Warping using a regular grid.
  # dg1=80 - use a grid with layers spaced at a minimum of 80 samples apart,
  # including the first and last time samples. This is not optimal, but is
  # simple. The structure aligned grid should be more accurate.
  # Fasle to read shifts from disk, True to find shifts.
  ur = goShiftsREG(pp,ppw,ps1,dw,80,False)
  
  plotErrors(ppw,ps1,dw,ur," REG vs SAG",u2=us) # Compare time shifts on error plot
             
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

def goShiftsREG(pp,ppw,ps1,dw,dg1,find):
  desc = " regular grid"
  label = "reg"
  printConstraints(desc,dg1)
  g1 = zeroint(ne1,n2,n3)
  g11 = Subsample.subsample(ne1,dg1) # simple regular interval grid
  g2 = Subsample.subsample(n2,dg2) # simple regular interval grid
  g3 = Subsample.subsample(n3,dg3) # simple regular interval grid
  for i3 in range(n3):
    for i2 in range(n2):
      g1[i3][i2] = copy(g11) # use the same grid g11 for all traces
  print "g1:"; dump(g11)
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  if find:
    # Computes time shifts u with values only at g1, g2, and g3 coordinates and
    # interpolate to get time shifts everywhere. Use Interp.MONOTONIC for
    # smoother shifts, but Interp.LINEAR is really the answer most consistent
    # with the method of finding time shifts. With MONOTONIC the min/max Vp/Vs
    # constraints may be violated. Interp.LINEAR for 2nd/3rd dimension should
    # not be changed.
    u = dw.findShifts(ppw,ps1,toFloats(g1),g2,g3,Interp.MONOTONIC,Interp.LINEAR)
    writeImage(threeDDir,"u_"+label,u)
  u = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  checkShifts3(u)
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  plotGrid(pp,coordMap,desc)
  plotShifts(pp,u,desc)
  plotVpvs(pp,u,desc)
  plotWarped(pp,ps1,dw,u,desc)
  return u
    
def goShiftsSAG(pp,ppw,ps1,dw,dg1,ng,find):
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
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  if find:
    # Computes time shifts u with values only at g1, g2, and g3 coordinates and
    # interpolate to get time shifts everywhere. Use Interp.MONOTONIC for
    # smoother shifts, but Interp.LINEAR is really the answer most consistent
    # with the method of finding time shifts. With MONOTONIC the min/max Vp/Vs
    # constraints may be violated. Interp.LINEAR for 2nd/3rd dimension should
    # not be changed.
    u = dw.findShifts(ppw,ps1,getG1Floats(x1m,g1Flat),g2,g3,Interp.MONOTONIC,
                      Interp.LINEAR)
    writeImage(threeDDir,"u_"+label,u)
  u = readImage(threeDDir,"u_"+label,ne1,n2,n3)
  checkShifts3(u)
  coordMap = Viewer3P.getSparseCoordsMap(g1,g2,g3,d1,d2,d3)
  plotGrid(pp,coordMap,desc)
  plotShifts(pp,u,desc)
  plotVpvs(pp,u,desc)
  plotWarped(pp,ps1,dw,u,desc)
  return u

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

def getG1Floats(x1m,g1Flat):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        g1[i3][i2][i1] = x1m[i3][i2][g1Flat[i1]]
  return g1

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

def setPlotVars(n1,nl):
  global ne1,nel,se,su,el,ul,xl,yl
  ne1 = n1
  nel = nl
  se = Sampling(ne1,d1,f1) # error sampling
  su = Sampling(nel,d1,f1) # shift sampling
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

def plotShifts(pp,u,desc):
  # Plot 2D shifts
  us = mul(u,d1)
  zm = ZeroMask(pp)
  us = DynamicWarpingC.extrapolate(len(pp[0][0]),us)
  zm.apply(0.00,us)
  plotPP3(us,title="Vertical shifts"+desc,s1=s1,s2=s2,s3=s3,clips1=uClips,
          label1="PS time (s)",cbar="Vertical shift (s)",cmap1=jet,cbw=100,
          limits1=el,slices=slices,he0=he0)

def plotVpvs(pp,u,desc,coordMap=None):
  # Plot Vp/Vs ratios
  zm = ZeroMask(pp)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(len(pp[0][0]),vpvs)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  cm = None
  if coordMap:
    cm = [coordMap,"rO",6.0]
  plotPP3(vpvs,coordMap=cm,title="Vp/Vs"+desc,s1=s1,s2=s2,s3=s3,
          clips1=vClips,label1="PP time (s)",cbar="Vp/Vs ratio",cmap1=jet,
          cbw=100,limits1=el,limits2=xl,limits3=yl,slices=slices,he0=he0)

def plotWarped(pp,ps1,dw,u,desc):
  # Plot Warped PS1 image and the difference between the warped and input.
  psw = dw.applyShifts(n1,ps1,u)
  plotPP3(psw,title="PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
          cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices,
          he0=he0)
  nrms = DynamicWarpingC.computeNrms(ne1,pp,psw)
  print desc+": nrms=%g"%nrms
  # d = sub(psw,pp)
  # plotPP3(d,title="PP-PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
  #         cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices)

def plotErrors(ppw,ps1,dw,u,desc,u2=None,u3=None,x1u=None,x2u=None):
  # Plot Alignment Errors with coarse grid points and shifts
  u = mul(u[k3],d1)
  p = [u,"shifts","w-",2.0]
  p2,p3 = None,None
  if u2:
    u2 = mul(u2[k3],d1)
    p2 = [u2,"shifts 2","c--",2.0]
  if u3:
    u3 = mul(u3[k3],d1)
    p3 = [u3,"shifts 3","r-.",2.0]
  if x1u and x2u: 
    x1u = mul(x1u,d1)
    x2u = mul(x2u,d1)
  e = dw.computeErrors2(ppw[k3],ps1[k3])
  e = DynamicWarpingC.transposeLag12(e)
  DynamicWarpingC.normalizeErrors(e)
  plot3(e,p=p,p2=p2,p3=p3,x11=x1u,x21=x2u,title="AE"+desc,s1=se,s2=su,
        hLabel="PP time (s)",vLabel="Vertical shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1rx2u)

###############################################################################
run(main)
