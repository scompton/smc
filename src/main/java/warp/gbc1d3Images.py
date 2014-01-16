#############################################################################
# 2D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
s1pp,s2,s3 = getSamplings() # Samplings defined in gbcUtils
n1pp = s1pp.getCount() # n1pp is the fastest dimension
n2 = s2.getCount()
n3 = s3.getCount() # n3 is the slowest

# default clips for seismic images(i), time shifts(u), interval vp/vp(v)
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.5,2.1]
w,h = 1200,600
# iMap = bwr
iMap = gray
i2,i3 = 80,82 # central trace
hw2 = 10
hw3 = 10

# just a guess... This controls the maximum possible time shift and the maximum
# PP time sample that could correspond to the last PS time sample. A larger
# value increases the maximum possible shift and includes less PP samples for
# warping.
vpvsAvg = 2.0

# Contstraints for time shift slopes, which are physically related to interval
# Vp/Vs ratios. Compute slope parameters from vpvsMin
vpvsMin,vpvsMax = 1.5,2.0
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2
rSmin =  -0.0
rSmax =  0.1
uSmin = 0.0
uSmax = 0.04

# ppName = "pp_smooth"
# ps1Name = "ps1_fkk_smooth"
# ps2Name = "ps2_fk_smooth"
ppName = "pp"
ps1Name = "ps1"
ps2Name = "ps2"
#############################################################################

def main(args):
  displayGBC()
  # displayGBC(ppName,ps1Name,ps2Name)
  #go3(ppName,ps1Name,ps2Name,75,75,100,hw2,hw3)
  # go3C(ppName,ps1Name,ps2Name,75,hw2,hw3)

def go3(ppName,ps1Name,ps2Name,dg1,dg2,dgS,hw2,hw3):
  # sf,u1e = goPpPs1Envelope(ppName,ps1Name,dg1,hw2,hw3)
  sf,u1 = goPpPs1(ppName,ps1Name,dg1,hw2,hw3)
  sg,uS = goPs1Ps2(ps1Name,ps2Name,dgS,hw2,hw3)
  uSPP = WarpUtils.ps1ToPpTime(sf,u1,sg,uS)
  u2c = WarpUtils.compositeShifts(sf,u1,uSPP)
  sf,u2 = goPpPs2(ppName,ps2Name,dg2,hw2,hw3,uc=u2c)
  plot1([[u1,"k-"],[u2,"r-"],[u2c,"b-"]],sf,title="shifts",width=400,height=800,
        o=x1dx2r)
  gammaS = WarpUtils.gammaS(sf,u1,u2)
  gammaSc = WarpUtils.gammaS(sf,u1,u2c)
  plotGammaS([[gammaS,"k-"]],sf,"PP",desc="")
  plotGammaS([[gammaSc,"k-"]],sf,"PP",desc="composite")

def go3C(ppName,ps1Name,ps2Name,dg1,hw2,hw3):
  desc1 = "PP & PS1 C3"
  desc2 = "PP & PS2 C3"
  descS = "PS1 & PS2 C3"
  f,sf =  getPp(i3,hw3,i2,hw2, ppName),s1pp
  g,sg = getPs1(i3,hw3,i2,hw2,ps1Name) # returns data and sampling
  h,sh = getPs2(i3,hw3,i2,hw2,ps2Name) # returns data and sampling
  dw = DynamicWarpingC3.fromVpVs(sf,sg,vpvsAvg,0.0,uSmin,uSmax)
  dw.setStrainLimits(r1min,r1max,rSmin,rSmax)
  dw.setInterpolationMethod(Method.LINEAR)
  se = dw.getSampling1()
  su1 = dw.getSamplingU1()
  suS = dw.getSamplingUS()
  g1 = Subsample.subsampleEnvelope(se.getCount(),f,dg1)
  g1 = Subsample.indicesToSampling(sf,g1)
  sw = Stopwatch()
  sw.restart()
  e = dw.computeErrorsSum(sf,f,sg,g,sh,h)
  sw.stop()
  print "Compute error time: %g seconds"%sw.time()
  sw.restart()
  u1,uS = dw.findShifts(sf,e,g1)
  sw.stop()
  print "Find shift time: %g seconds"%sw.time()
  # u2 = WarpUtils.computeU2(sf,u1,sg,uS)
  u2 = WarpUtils.compositeShifts(sf,u1,uS)
  gpp = WarpUtils.ps1ToPpTime(sf,u1,sg,g[hw3][hw2])
  # hpp = WarpUtils.ps1ToPpTime(sf,u1,sh,h[hw3][hw2])
  # print "gpp: min=%g, max=%g"%(min(gpp),max(gpp))
  # print "hpp: min=%g, max=%g"%(min(hpp),max(hpp))
  w1 = WarpUtils.applyShifts(sf,u1,sg,g[hw3][hw2])
  w2 = WarpUtils.applyShifts(sf,u2,sh,h[hw3][hw2])
  wh = WarpUtils.applyShifts(sf,u1,sh,h[hw3][hw2])
  # wS = WarpUtils.applyShifts(sg,uS,sh,h[hw3][hw2])
  wS = WarpUtils.applyShifts(sf,uS,sf,wh)
  print "NRMS 1: %g"%WarpUtils.computeNrms(se.getCount(),f[hw3][hw2],w1)
  print "NRMS 2: %g"%WarpUtils.computeNrms(se.getCount(),f[hw3][hw2],w2)
  # print "NRMS S: %g"%WarpUtils.computeNrms(se.getCount(),g[hw3][hw2],wS)
  print "NRMS S: %g"%WarpUtils.computeNrms(se.getCount(),w1,wS)
  # Plot warped traces
  plotWarped([[f[hw3][hw2],"k-"],[w1,"r-"]],sf,se,"PP",desc1)
  plotWarped([[f[hw3][hw2],"k-"],[w2,"b-"]],sf,se,"PP",desc2)
  # plotWarped([[g[hw3][hw2],"b-"],[wS,"m-"]],sf,se,"PS1",descS)
  plotWarped([[w1,"b-"],[wS,"m-"]],sf,se,"PP",descS)
  plotVpVs(u1,sf,"PP",desc1)
  plotVpVs(u2,sf,"PP",desc2)
  gammaS = WarpUtils.gammaS(sf,u1,u2)
  plotGammaS([[gammaS,"k-"]],sf,"PP")
  plotError3(e,se,su1,suS,u1,uS)

def goPpPs1(ppName,ps1Name,dg1,hw2,hw3):
  desc = "PP & PS1"
  f,sf = getPp(i3,hw3,i2,hw2,ppName),s1pp
  g,sg = getPs1(i3,hw3,i2,hw2,ps1Name) # returns data and sampling
  dw = DynamicWarpingC.fromVpVs(sf,sg,vpvsAvg,0.0)
  dw.setStrainLimits(r1min,r1max)
  dw.setInterpolationMethod(Method.LINEAR)
  e,u = goWarp(sf,f,sg,g,dg1,dw)
  se = dw.getSampling1()
  su = dw.getSamplingU()
  uA = [[copy(se.getCount(),u),"u","w-",2.0]]
  w = WarpUtils.applyShifts(sf,u,sg,g[hw3][hw2])
  print "NRMS 1: %g"%WarpUtils.computeNrms(se.getCount(),f[hw3][hw2],w)
  wA = [[f[hw3][hw2],"k-"],[w,"r-"]]
  plotError2(e,se,su,uA,"PP",desc)
  plotWarped(wA,sf,se,"PP",desc)
  plotVpVs(u,sf,"PP",desc)
  return sf,u

def goPpPs1Envelope(ppName,ps1Name,dg1,hw2,hw3):
  desc = "PP & PS1 Envelope"
  f,sf = getPp(i3,hw3,i2,hw2,ppName),s1pp
  g,sg = getPs1(i3,hw3,i2,hw2,ps1Name) # returns data and sampling
  fe = Envelope.getEnvelope(f)
  ge = Envelope.getEnvelope(g)
  plot3(fe,title="f envelope")
  plot3(ge,title="g envelope")
  dw = DynamicWarpingC.fromVpVs(sf,sg,vpvsAvg,0.0)
  dw.setStrainLimits(r1min,r1max)
  dw.setInterpolationMethod(Method.LINEAR)
  e,u = goWarp(sf,fe,sg,ge,dg1,dw)
  se = dw.getSampling1()
  su = dw.getSamplingU()
  uA = [[copy(se.getCount(),u),"u","w-",2.0]]
  w = WarpUtils.applyShifts(sf,u,sg,g[hw3][hw2])
  wA = [[f[hw3][hw2],"k-"],[w,"r-"]]
  plotError2(e,se,su,uA,"PP",desc)
  plotWarped(wA,sf,se,"PP",desc)
  plotVpVs(u,sf,"PP",desc)
  return sf,u

def goPpPs2(ppName,ps2Name,dg1,hw2,hw3,uc=None):
  desc = "PP & PS2"
  f,sf = getPp(i3,hw3,i2,hw2,ppName),s1pp
  h,sh = getPs2(i3,hw3,i2,hw2,ps2Name) # returns data and sampling
  dw = DynamicWarpingC.fromVpVs(sf,sh,vpvsAvg,0.0)
  dw.setStrainLimits(r1min,r1max)
  dw.setInterpolationMethod(Method.LINEAR)
  e,u = goWarp(sf,f,sh,h,dg1,dw)
  se = dw.getSampling1()
  su = dw.getSamplingU()
  if uc:
    uA = [[copy(se.getCount(),u),"u","w-",2.0],
          [copy(se.getCount(),uc),"uc","y-",2.0]]
  else:
    uA = [[copy(se.getCount(),u),"u","w-",2.0]]
  w = WarpUtils.applyShifts(sf,u,sh,h[hw3][hw2])
  print "NRMS 2: %g"%WarpUtils.computeNrms(se.getCount(),f[hw3][hw2],w)
  wA = [[f[hw3][hw2],"k-"],[w,"b-"]]
  plotError2(e,se,su,uA,"PP",desc)
  plotWarped(wA,sf,se,"PP",desc)
  plotVpVs(u,sf,"PP",desc)
  return sf,u

def goPs1Ps2(ps1Name,ps2Name,dg1,hw2,hw3):
  desc = "PS1 & PS2"
  g,sg = getPs1(i3,hw3,i2,hw2,ps1Name) # returns data and sampling
  h,sh = getPs2(i3,hw3,i2,hw2,ps2Name) # returns data and sampling
  dw = DynamicWarpingC(uSmin,uSmax,sg)
  dw.setStrainLimits(rSmin,rSmax)
  dw.setInterpolationMethod(Method.LINEAR)
  e,u = goWarp(sg,g,sh,h,dg1,dw)
  se = dw.getSampling1()
  su = dw.getSamplingU()
  uA = [[copy(se.getCount(),u),"u","w-",2.0]]
  w = WarpUtils.applyShifts(sg,u,sh,h[hw3][hw2])
  print "NRMS S: %g"%WarpUtils.computeNrms(se.getCount(),g[hw3][hw2],w)
  wA = [[g[hw3][hw2],"b-"],[w,"m-"]]
  plotError2(e,se,su,uA,"PS1",desc)
  plotWarped(wA,sg,se,"PS1",desc)
  return sg,u

def getPp(i3,hw3,i2,hw2,ppName):
  pp = getGbcImage(baseDir,ppName)
  ppc = getCube(pp,i3,hw3,i2,hw2)
  # plot3(ppc,title="PP Subset")
  return ppc

def getPs1(i3,hw3,i2,hw2,ps1Name):
  ps1 = getGbcImage(baseDir,ps1Name)
  n1ps = len(ps1[0][0])
  sps1 = Sampling(n1ps,s1pp.getDelta(),s1pp.getFirst())
  ps1c = getCube(ps1,i3,hw3,i2,hw2)
  # plot3(ps1c,title="PS1 Subset")
  return ps1c,sps1

def getPs2(i3,hw3,i2,hw2,ps2Name):
  ps2 = getGbcImage(baseDir,ps2Name)
  n1ps = len(ps2[0][0])
  sps2 = Sampling(n1ps,s1pp.getDelta(),s1pp.getFirst())
  ps2c = getCube(ps2,i3,hw3,i2,hw2)
  # plot3(ps2c,title="PS2 Subset")
  return ps2c,sps2

def goWarp(sf,f,sg,g,dg1,dw):
  se = dw.getSampling1()
  e = dw.computeErrorsSum(sf,f,sg,g)
  g1 = Subsample.subsampleEnvelope(se.getCount(),f,dg1)
  g1 = Subsample.indicesToSampling(sf,g1)
  u = dw.findShifts(sf,e,g1)
  return e,u

def getCube(f,i3,hw3,i2,hw2):
  n2c = 2*hw2+1
  n3c = 2*hw3+1
  i2s,i2e = i2-hw2,i2+hw2
  i3s,i3e = i3-hw3,i3+hw3
  # print "i2s=%d, i2e=%d, i3s=%d, i3e=%d"%(i2s,i2e,i3s,i3e)
  if i2s<0:
    is2 = 0
  if i2e>n2-1:
    i2e = n2-1
  if i3s<0:
    i3s = 0
  if i3e>n3-1:
    i3e = n3-1
  fc = zerofloat(n1pp,n2c,n3c)
  ii = 0
  for i in range(i3s,i3e+1):
    jj = 0
    for j in range(i2s,i2e+1):
      fc[ii][jj] = f[i][j]
      jj += 1
    ii += 1
  return fc

# Copy a 1D array f into a 3D array with dimensions in g
def copy1To3(f,g):
  n1 = len(g[0][0])
  n2 = len(g[0])
  n3 = len(g)
  f3 = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      f3[i3][i12] = copy(f)
  return f3

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

#############################################################################
# Plotting (most in gbcUtils.py)

def plotError2(e,se,su,uA,hname,desc):
  e = Transpose.transpose12(e)
  WarpUtils.normalizeErrors(e)
  el = [se.getFirst(),se.getLast()]
  ul = [su.getFirst(),su.getLast()]
  plot2(e,pA=uA,title="AE "+desc,s1=se,s2=su,hLabel=hname+" time (s)",
        vLabel="Vertical shift (s)",cbar="Error",clips1=[0,0.2],width=1200,
        height=600,hLimits=el,vLimits=ul,cbw=100,o=x1rx2u)

def plotError3(e,se,su1,suS,u1,uS):
  e = Transpose.transpose132(e)
  e = Transpose.transpose312(e)
  WarpUtils.normalizeErrors(e)
  xyz = getXYZ(se,u1,uS)
  rgb = getRGB(se.getCount(),Color.YELLOW)
  plot3D(e,se,su1,suS,xyz=xyz,rgb=rgb,clips=[0.0,0.3])

def plotWarped(wA,s,se,vname,desc):
  el = [se.getFirst(),se.getLast()]
  plot1(wA,title=desc+" warped",s=s,vLabel=vname+" time (s)",
        width=400,height=900,vLimits=el,o=x1dx2r)

def plotVpVs(u,s,vname,desc):
  vpvs = WarpUtils.vpvs(s,u)
  plot1([[vpvs,"k-"]],title=desc,s=s,vLabel=vname+" time (s)",hLabel="Vp/Vs",
        width=400,height=900,hLimits=[vpvsMin,vpvsMax],o=x1dx2r)

def plotGammaS(gA,s,vname,desc=""):
  title="gammaS"
  if desc:
    title = title+" "+desc
  plot1(gA,title="gammaS "+desc,s=s,vLabel=vname+" time (s)",
      hLabel="",width=400,height=900,hLimits=[-0.1,0.1],o=x1dx2r)

def getXYZ(se,u1,uS):
  n1 = se.getCount()
  xyz = zerofloat(n1*3)
  j = 0
  for i1 in range(n1):
    xyz[j  ] = uS[i1]
    xyz[j+1] = u1[i1]
    xyz[j+2] = se.getValue(i1)
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
    
def plot3D(f,s1,s2,s3,xyz=None,size=0.004,rgb=None,clips=None,cmap=None,width=800,
           height=1000):
  sf = SimpleFrame()
  ipg = sf.addImagePanels(s1,s2,s3,f)
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
  ov.setScale(2.0)
  sf.setVisible(True)

###############################################################################
run(main)
