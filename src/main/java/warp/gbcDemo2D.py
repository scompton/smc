#############################################################################
# 2D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
twoDDir = baseDir+"2D/" # directory for writing intermediate files
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
vClips = [1.5,2.1]
# iMap = bwr
iMap = gray
i3 = 72 # index of 2D slice to extract from 3D volume

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

# Constraints for horizontal direction. A better match between the PS and PP
# images can be achieved by loosening these constraints, but I'm not sure that
# we want to correct for all lateral variations seen in the PS image but not
# in the PP image. These differences may not be related to Vp/Vs ratios.
r2min,r2max,dg2 = -0.15,0.15,15

#############################################################################

def main(args):
  # pp,ps1 = getImages2D() # FKK and AGC
  pp,ps1 = getImages2DSmooth() # FKK, AGC, and Structure Oriented Smoothing
  lMax = DynamicWarpingC.computeMaxLag(n1ps,vpvsAvg)
  n1w = DynamicWarpingC.computeMaxLength(n1,n1ps,vpvsAvg)
  ppw = copy(n1w,n2,pp)
  dw = DynamicWarpingC(0,lMax) # DynamicWarpingC.java
  nl = dw.getNumberOfLags()
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(n1w,nl)

  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  us = goShiftsSAG(pp,ppw,ps1,dw,50,13)

  # Warping using a regular grid.
  # dg1=80 - use a grid with layers spaced at a minimum of 80 samples apart,
  # including the first and last time samples. This is not optimal, but is
  # simple. The structure aligned grid should be more accurate.
  ur = goShiftsREG(pp,ppw,ps1,dw,80)

  # Compare time shifts on error plot
  plotErrors(ppw,ps1,dw,ur," REG vs SAG",u2=us) 
             
def getImages2D():
  pp = getGbcImage(baseDir,"pp")
  pp = pp[i3]
  # makeSubset(baseDir,"ps1_fkk",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_1501_150_145",1501)
  ps1 = ps1[i3]
  plot2(ps1,title="PS1",s1=s1ps,s2=s2,clips1=iClips,vLabel="PS time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  return pp,ps1

def getImages2DSmooth():
  pp = getGbcImage(baseDir,"pp_smooth")
  pp = pp[i3]
  # makeSubset(baseDir,"ps1_fkk_smooth",0,n1ps-1,0,n2-1,0,n3-1)
  ps1 = getGbcImage(baseDir,"ps1_fkk_smooth_1501_150_145",1501)
  ps1 = ps1[i3]
  plot2(ps1,title="PS1",s1=s1ps,s2=s2,clips1=iClips,vLabel="PS time (s)",
        vInterval=0.2,cmap1=iMap,cbw=100)
  return pp,ps1

def printConstraints(desc,dg1):
  k1min = int(math.ceil( r1min*dg1))
  k1max = int(math.floor(r1max*dg1))
  k2min = int(math.ceil( r2min*dg2))
  k2max = int(math.floor(r2max*dg2))
  info = """Constraints: %s
  r1min=%g, r1max=%g, dg1=%g, k1min=%g, k1max=%g
  r2min=%g, r2max=%g, dg2=%g, k2min=%g, k2max=%g"""%(desc,r1min,r1max,dg1,k1min,
                                                     k1max,r2min,r2max,dg2,
                                                     k2min,k2max)
  print info

def goShiftsREG(pp,ppw,ps1,dw,dg1):
  desc = " regular grid"
  printConstraints(desc,dg1)
  g1 = zeroint(ne1,n2)
  g11 = Subsample.subsample(ne1,dg1); # simple regular interval grid
  g2 = Subsample.subsample(n2,dg2); # simple regular interval grid
  for i2 in range(n2):
    g1[i2] = copy(g11) # use the same grid g11 for all traces
  print "g1:"; dump(g11)
  print "g2:"; dump(g2)
  # Compute time shifts u with values only at g1 and g2 coordinates and 
  # interpolate to get time shifts everywhere. Use Interp.MONOTONIC for smoother
  # shifts, but Interp.LINEAR is really the answer most consistent with the
  # method of finding time shifts. With MONOTONIC the min/max Vp/Vs constraints
  # may be violated. Interp.LINEAR for 2nd dimension should not be changed.
  u = dw.findShifts(ppw,ps1,toFloats2(g1),g2,Interp.MONOTONIC,Interp.LINEAR)
  checkShifts2(u)
  x1u,x2u = getShiftCoords(g1,g2,u)
  x1i,x2i = getSparseGridCoords(g1,g2);
  plotGrid(pp,x1i,x2i,desc)
  plotShifts(pp,u,desc)
  plotVpvs(pp,u,desc)
  plotWarped(pp,ps1,dw,u,desc)
  return u
    
def goShiftsSAG(pp,ppw,ps1,dw,dg1,ng):
  desc = " structure aligned grid"
  printConstraints(desc,dg1)
  goSlopes(pp) # compute slopes for flattening
  goFlat(pp) # flatten
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
  u = dw.findShifts(ppw,ps1,getG1Floats(x1m,g1Flat),g2,Interp.MONOTONIC,
                    Interp.LINEAR)
  checkShifts2(u)
  g1f = zeroint(ne1,n2)
  for i2 in range(n2):
    g1f[i2] = copy(g1Flat)
  x1u,x2u = getShiftCoords(g1,g2,u)
  plotErrors(ppw,ps1,dw,u,desc,x1u=x1u,x2u=x2u)
  x1i,x2i = getSparseGridCoords(g1,g2);
  plotGrid(pp,x1i,x2i,desc)
  plotShifts(pp,u,desc)
  plotVpvs(pp,u,desc)
  plotWarped(pp,ps1,dw,u,desc)
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

def getG1Floats(x1m,g1Flat):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2)
  for i2 in range(n2):
    for i1 in range(ng):
      g1[i2][i1] = x1m[i2][g1Flat[i1]]
  return g1
  
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
    
def goSlopes(f):
  fs = zerofloat(ne1,n2) # f with only ne1 samples
  for i2 in range(n2):
    for i1 in range(ne1):
      fs[i2][i1] = f[i2][i1]
  s1f = Sampling(len(fs[0]))
  s2f = Sampling(len(fs))
  sigma1 = 8
  sigma2 = 12
  print "LSF: sigma1=%g, sigma2=%g"%(sigma1,sigma2)
  pmax = 2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(s1f.getCount(),s2f.getCount())
  el = zerofloat(s1f.getCount(),s2f.getCount())
  lsf.findSlopes(fs,p2,el)
  el = pow(el,4)
  writeImage(twoDDir,"p2",p2)
  writeImage(twoDDir,"el",el)
  # plot2(el,title="Linearity",s1=s1f,hLabel="crossline",vLabel="time (samples)",
  #       cbar="linearity",cmap1=jet,clips1=[0.0,1.0],width=600,height=900,
  #       cbw=100)
  # plot2(p2,title="Slopes",s1=s1f,hLabel="crossline",vLabel="time (samples)",
  #       cbar="slope",cmap1=bwr,width=600,height=900,cbw=100)
  
def goFlat(f):
  fs = zerofloat(ne1,n2) # f with only ne1 samples
  for i2 in range(n2):
    for i1 in range(ne1):
      fs[i2][i1] = f[i2][i1]
  s1f = Sampling(len(fs[0]))
  s2f = Sampling(len(fs))
  p2 = readImage(twoDDir,"p2",ne1,n2)
  el = readImage(twoDDir,"el",ne1,n2)
  fl = Flattener2()
  fl.setWeight1(0.2)
  fm = fl.getMappingsFromSlopes(s1f,s2f,p2,el)
  ff = fm.flatten(fs)
  h = fm.unflatten(ff)
  s = fm.getShiftsS()
  r = fm.getShiftsR()
  x1 = fm.x1
  y1 = fm.u1
  writeImage(twoDDir,"x1",x1)
  writeImage(twoDDir,"y1",y1)
  writeImage(twoDDir,"ff",ff)
  writeImage(twoDDir,"s",s)
  writeImage(twoDDir,"r",r)
  # plot2(ff,title="PP Flat",hLabel="crossline",vLabel="Tau",cbw=100)
  # plot2(h,title="PP Unflattened",s1=s1f,hLabel="crossline",
  #       vLabel="time (samples)",cbw=100)
  # plot2(s,title="PP Flat Shifts",s1=s1f,hLabel="crossline",
  #       vLabel="time (samples)",cbar="shift (samples)",cmap1=jet,clips1=None,
  #       cbw=100)
  print "average shift =",sum(s)/(n1*n2),"samples"

#############################################################################

def setPlotVars(n1,nl):
  global ne1,nel,se,su,el,ul,xl
  ne1 = n1
  nel = nl
  se = Sampling(ne1,d1,f1) # error sampling
  su = Sampling(nel,d1,f1) # shift sampling
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  # el = [0.0,1.0] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  
def plotGrid(pp,x1i,x2i,desc):
  # Plot PP with grid points
  x1i = mul(x1i,d1)
  # x2i = add(25,x2i)
  x2i = mul(x2i,d2)
  s1 = Sampling(len(pp[0]),d1,f1)
  plot2(pp,x12=x1i,x22=x2i,title="PP",s1=s1,s2=s2,vLabel="PP time (s)",
        hLimits=xl,vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotShifts(pp,u,desc):
  # Plot 2D shifts
  us = mul(u,d1)
  zm = ZeroMask(pp)
  us = DynamicWarpingC.extrapolate(len(pp[0]),us)
  zm.apply(0.00,us)
  plot2(us,title="Vertical shifts"+desc,s1=s1,s2=s2,vLabel="PP time (s)",
        cbar="Vertical shift (s)",cmap1=jet,clips1=uClips,vLimits=el,
        cbw=100)

def plotVpvs(pp,u,desc,x1i=None,x2i=None):
  # Plot Vp/Vs ratios
  zm = ZeroMask(pp)
  if x1i and x2i:
    x1i = mul(x1i,d1)
    x2i = mul(x2i,d2)
  vpvs = DynamicWarpingC.vpvs(u)
  vpvs = DynamicWarpingC.extrapolate(len(pp[0]),vpvs)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plot2(pp,g=vpvs,x12=x1i,x22=x2i,title="Vp/Vs"+desc,s1=s1,s2=s2,
        vLabel="PP time (s)",cbar="Vp/Vs ratio",clips1=iClips,cmap2=jet,
        clips2=vClips,hLimits=xl,vLimits=el,cmap1=iMap,cbw=100)

def plotWarped(pp,ps1,dw,u,desc):
  # Plot Warped PS1 image and the difference between the warped and input.
  psw = dw.applyShifts(n1,ps1,u)
  plot2(psw,title="PS1 warped"+desc,s1=s1,s2=s2,vLabel="PP time (s)",
        vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)
  nrms = DynamicWarpingC.computeNrms(ne1,pp,psw)
  print desc+": nrms=%g"%nrms
  # d = sub(psw,pp)  
  # plot2(d,title="PP-PS1 warped"+desc,s1=s1,s2=s2,vLabel="PP time (s)",
  #       vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotErrors(ppw,ps,dw,u,desc,u2=None,u3=None,x1u=None,x2u=None):
  # Plot Alignment Errors with coarse grid points and shifts
  u = mul(u,d1)
  p = [u,"shifts","w-",2.0]
  p2,p3 = None,None
  if u2:
    u2 = mul(u2,d1)
    p2 = [u2,"shifts 2","c--",2.0]
  if u3:
    u3 = mul(u3,d1)
    p3 = [u3,"shifts 3","r-.",2.0]
  if x1u and x2u:
    x1u = mul(x1u,d1)
    x2u = mul(x2u,d1)
  e = dw.computeErrors2(ppw,ps)
  e = DynamicWarpingC.transposeLag12(e)
  DynamicWarpingC.normalizeErrors(e)
  plot3(e,p=p,p2=p2,p3=p3,x11=x1u,x21=x2u,title="AE"+desc,s1=se,s2=su,
        hLabel="PP time (s)",vLabel="Vertical shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1rx2u)

###############################################################################
run(main)
