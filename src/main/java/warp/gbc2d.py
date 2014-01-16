#############################################################################
# 2D Dynamic warping for GBC images

from gbcUtils import *

#############################################################################
baseDir = getBaseDir() # defined in gbcUtils
twoDDir = baseDir+"2D/" # directory for writing intermediate files
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
i3 = 72 # index of 2D slice to extract from 3D volume

# just a guess... This controls the maximum possible time shift and the maximum
# PP time sample that could correspond to the last PS time sample. A larger
# value increases the maximum possible shift and includes less PP samples for
# warping.
vpvsAvg = 1.9

# Contstraints for time shift slopes, which are physically related to interval
# Vp/Vs ratios. Compute slope parameters from vpvsMin/Max.
vpvsMin,vpvsMax = 1.5,2.5
# vpvsMin,vpvsMax = 1.0,2.5
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

# Constraints for horizontal direction. A better match between the PS and PP
# images can be achieved by loosening these constraints, but I'm not sure that
# we want to correct for all lateral variations seen in the PS image but not
# in the PP image. These differences may not be related to Vp/Vs ratios.
r2min,r2max,dg2 = -0.15,0.15,15

#############################################################################

def main(args):
  goSagVsReg()
  # goSagFixedU()
  # goErrorTest()

def goSagVsReg():
  pp,ps1,dw = getPpPs1Data()

  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  us = goShiftsSAG(s1pp,pp,s1ps,ps1,dw,50,13,"PP","PS1")

  # Warping using a regular grid.
  # dg1=80 - use a grid with layers spaced at a minimum of 80 samples apart,
  # including the first and last time samples. This is not optimal, but is
  # simple. The structure aligned grid should be more accurate.
  ur = goShiftsREG(s1pp,pp,s1ps,ps1,dw,80,"PP","PS1")

  # Compare time shifts on error plot
  # uA = [[us,"us","w-",2.0],[ur,"ur","w--",2.0]]
  # plotErrors(s1pp,pp,s1ps,ps1,dw,uA," REG vs SAG","PP")

def goSagFixedU():
  pp,ps1,dw = getPpPs1Data()

  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  us = goShiftsSAGFixedU(pp,ps1,dw,50,13,"PP","PS1")

def goErrorTest():
  pp,ps1,dw = getPpPs1Data()
  e = dw.computeErrors2(pp,ps1)
  e = WarpUtils.transposeLag12(e)
  WarpUtils.normalizeErrors(e)
  plot3(e,title="Alignment Errors",s1=se,s2=su,hLabel="PP time (s)",
        vLabel="Vertical shift (s)",cbar="Error",clips1=[0,0.15],width=900,
        height=600,hLimits=el,vLimits=ul,o=x1rx2u)
  es = dw.computeErrorsSum2(pp,ps1)
  es = WarpUtils.transposeLag(es)
  WarpUtils.normalizeErrors(es)
  plot2(es,title="Alignment Errors Stack",s1=se,s2=su,hLabel="PP time (s)",
        vLabel="Vertical shift (s)",cbar="Error",clips1=[0,0.6],width=900,
        height=600,hLimits=el,vLimits=ul,o=x1rx2u)

def getPpPs1Data():
  # pp,ps1 = getImages2D() # FKK and AGC
  pp,ps1 = getImages2DSmooth() # FKK, AGC, and Structure Oriented Smoothing
  dw = DynamicWarpingC.fromVpVs(s1pp,s1ps,vpvsAvg,0.0,s2)
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(dw,n1pp,s1pp)
  return pp,ps1,dw

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

def goShiftsREG(sf,f,sg,g,dw,dg1,fname,gname):
  desc = " regular grid"
  printConstraints(desc,dg1)
  g11 = Subsample.subsample(ne1,dg1); # simple regular interval grid
  ng1 = len(g11)
  g11f = zerofloat(ng1)
  for i1 in range(ng1):
    g11f[i1] = se.getValue(g11[i1])
  g1 = zerofloat(ng1,n2)
  for i2 in range(n2):
    g1[i2] = copy(g11f)
  g2 = getG2(dg2)
  print "g1[0]:"; dump(g1[0])
  print "g2:"; dump(g2)
  return goShifts(sf,f,sg,g,g1,g2,dw,fname,gname,desc)

def goShiftsSAG(sf,f,sg,g,dw,dg1,ng,fname,gname):
  desc = " structure aligned grid"
  printConstraints(desc,dg1)
  goSlopes(f) # compute slopes for flattening
  goFlat(f) # flatten
  x1m = readImage(twoDDir,"x1",ne1,n2) # flattening mappings
  ff = readImage(twoDDir,"ff",ne1,n2) # flattened image
  g1,g1Flat = makeAutomaticGrid(x1m,ff,dw,dg1,ng=ng) # structure aligned grid
  g1 = getG1Floats(x1m,g1Flat,sf.getDelta())
  g2 = getG2(dg2)
  print "g2:"; dump(g2)
  return goShifts(sf,f,sg,g,g1,g2,dw,fname,gname,desc)

def goShifts(sf,f,sg,g,g1,g2,dw,fname,gname,desc):
  dw.setInterpolationMethod(Method.MONOTONIC)
  dw.setWorkTracker(WarperProgress())
  u = dw.findShifts(sf,f,sg,g,g1,g2)
  checkShifts2(u)
  une1 = copy(ne1,n2,u)
  x1u,x2u = getShiftCoords(sf,g1,g2,u)
  x12SliceA = [[x1u,x2u,"coarse grid","rO",10.0]]
  uA = [[une1,"u","w-",2.0]]
  plotErrors(sf,f,sg,g,dw,uA,desc,fname,x12SliceA=x12SliceA)
  x1i,x2i = getSparseGridCoords(g1,g2)
  dump(x1i)
  dump(x2i)
  plotGrid(f,x1i,x2i,fname)
  plotWarped(sf,f,sg,g,u,desc,fname,gname)
  plotVpvs(f,u,desc,fname)
  plotShifts(f,u,desc,fname)
  return u

def goShiftsSAGFixedU(f,g,dw,dg1,ng,fname,gname):
  printConstraints("goShiftsSAGFixedU",dg1)
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
  g1f = getG1Floats(x1m,g1Flat)
  u = dw.findShifts(f,g,g1f,g2,Interp.MONOTONIC,Interp.LINEAR)

  # Set known shifts
  ksu = KnownShiftUtil(su,se,s2,filldouble(r1min,ne1))
  p1s1tA = [[0.68,0.937,0.833],[0.706,0.955,3.613],[1.916,2.616,2.399]]
  # p1s1tA = [[0.68,0.850,0.833]] # should fail, shift violates Vp/Vs constraints
  # test1DKS(f,g,dw,getG1Floats(x1m,g1Flat),p1s1tA,su,se,s2)
  nks = len(p1s1tA)
  print "nks=",nks
  x1ks = fillfloat(-1000,nks,n2)
  x2ks = fillfloat(-1000,nks,n2)
  for i in range(nks):
    p1s1t = p1s1tA[i]
    s = p1s1t[1]-p1s1t[0]
    ksu.add(s,p1s1t[0],p1s1t[2])
    x1ks[s2.indexOfNearest(p1s1t[2])][i] = p1s1t[0]
    x2ks[s2.indexOfNearest(p1s1t[2])][i] = s
  map = ksu.getMap()
  # dw.setStrainLimits(r1min,r1max,-0.2,0.2)
  print "g2:",dump(g2)
  g2ks = ksu.getG2(g2,filldouble(r2min,n2),filldouble(r2max,n2))
  print "g2ks:",dump(g2ks)
  # g1ks = ksu.getG1(g1f,filldouble(r1min,ne1),filldouble(r1max,ne1))
  g1ks = g1f
  uks = dw.findShifts(f,g,g1ks,g2ks,Interp.MONOTONIC,Interp.LINEAR,map)
  for i in range(nks):
    p1s1t = p1s1tA[i]
    s = p1s1t[1]-p1s1t[0]
    i2 = s2.indexOfNearest(p1s1t[2])
    i1 = se.indexOfNearest(p1s1t[0])
    csu,csuks = u[i2][i1]*d1,uks[i2][i1]*d1
    print "i2=%d, i1=%d: s=%g, csu=%g, csuks=%g"%(i2,i1,s,csu,csuks)
  checkShifts2(uks)
  checkShifts2(u)
  x1u,x2u = getShiftCoords(g1,g2,u)
  x1u,x2u = mul(x1u,d1),mul(x2u,d1)
  x1uks,x2uks= getShiftCoords(getG1Ints2(g1ks),g2,uks)
  x1uks,x2uks= mul(x1uks,d1),mul(x2uks,d1)
  x12SliceA = [[x1u,x2u,"coarse grid","rO",10.0],
               [x1uks,x2uks,"coarse grid ks","cO",10.0],
               [x1ks,x2ks,"known shifts","yO",10.0]]
  uA = [[mul(u,d1),"u","w-",2.0],[mul(uks,d1),"u with ks","w--",2.0]]
  plotErrors(f,g,dw,uA," SAG vs SAG KS",fname,x12SliceA=x12SliceA)
  desc = " SAG"
  x1i,x2i = getSparseGridCoords(g1,g2)
  plotGrid(f,x1i,x2i,fname)
  plotShifts(f,u,desc,fname)
  plotVpvs(f,u,desc,fname)
  plotWarped(f,g,u,desc,fname,gname)
  desc = " SAG KS"
  x1i,x2i = getSparseGridCoords(g1,g2)
  plotGrid(f,x1i,x2i,fname)
  plotShifts(f,u,desc,fname)
  plotVpvs(f,u,desc,fname)
  plotWarped(f,g,u,desc,fname,gname)
  return u

def test1DKS(ppw,ps1,dw,g1,p1s1tA,su,se,s2):
  nks = len(p1s1tA)
  r1minA = filldouble(r1min,ne1)
  for i in range(nks):
    p1s1t = p1s1tA[i]
    s = p1s1t[1]-p1s1t[0]
    ksu = KnownShiftUtil(su,se,r1minA)
    ksu.add(s,p1s1t[0])
    i2 = s2.indexOfNearest(p1s1t[2])
    x1ks = [p1s1t[0]]
    x2ks = [s]
    x1ks = jarray.array(x1ks,"f")
    x2ks = jarray.array(x2ks,"f")
    # g1ks = ksu.getG1(g1[i2],filldouble(r1min,ne1),filldouble(r1max,ne1))
    u = dw.findShifts(ppw[i2],ps1[i2],g1[i2],Interp.MONOTONIC)
    uks = dw.findShifts(ppw[i2],ps1[i2],g1[i2],Interp.MONOTONIC,ksu.getCoords())
    # plot = False
    # if i2==25:
    #   plot = True
    # ugks = dw.findShifts(ppw[i2],ps1[i2],g1ks,Interp.MONOTONIC,ksu.getCoords(),plot)
    i1 = se.indexOfNearest(p1s1t[0])
    # csu,csuks,csugks = u[i1]*d1,uks[i1]*d1,ugks[i1]*d1
    # print "i2 %d, i1 %d: s=%g, csu=%g, csuks=%g, csugks=%g"%(i2,i1,s,csu,csuks,
    #                                                          csugks)
    csu,csuks = u[i1]*d1,uks[i1]*d1
    print "i2 %d, i1 %d: s=%g, csu=%g, csuks=%g"%(i2,i1,s,csu,csuks)

    # print "g1:",dump(g1[i2])
    # print "1D g1ks:",dump(g1ks)
    e = dw.computeErrors(ppw[i2],ps1[i2])
    et = WarpUtils.transposeLag(e)
    WarpUtils.normalizeErrors(et)
    # eks = copy(e)
    # dw.fixShifts(eks,ksu.getCoords())
    # efrs = dw.getSmoothErrors(e,getG1Ints(g1[i2]),len(ppw[0]))
    # efrsks = dw.getSmoothErrors(eks,getG1Ints(g1ks),len(ppw[0]))
    # ef = WarpUtils.transposeLag(efrs[0])
    # WarpUtils.normalizeErrors(ef,-Float.MAX_VALUE,Float.MAX_VALUE)
    # plot2(ef,title="Accumulate forward i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.8],width=900,height=600,o=x1rx2u)
    # efks = WarpUtils.transposeLag(efrsks[0])
    # # print "efks before normalization: min=%g, max=%f"%(min(efks),max(efks))
    # # emax = -Float.MAX_VALUE
    # # for ie2 in range(len(efks)):
    # #   for ie1 in range(len(efks[0])):
    # #     ev = efks[ie2][ie1]
    # #     # if ev!=Float.MAX_VALUE:
    # #     if ev<1.0E20:
    # #       if ev>emax:
    # #         emax = ev
    # # print "emax besides Float.MAX_VALUE: %g"%emax
    # WarpUtils.normalizeErrors(efks,-Float.MAX_VALUE,1E20)
    # # print "efks after normalization: min=%g, max=%g"%(min(efks),max(efks))
    # plot2(efks,title="Accumulate forward KS i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.8],width=900,height=600,o=x1rx2u)
    # er = WarpUtils.transposeLag(efrs[1])
    # WarpUtils.normalizeErrors(er,-Float.MAX_VALUE,Float.MAX_VALUE)
    # plot2(er,title="Accumulate reverse i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.8],width=900,height=600,o=x1rx2u)
    # erks = WarpUtils.transposeLag(efrsks[1])
    # WarpUtils.normalizeErrors(erks,-Float.MAX_VALUE,1E20)
    # plot2(erks,title="Accumulate reverse KS i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.8],width=900,height=600,o=x1rx2u)
    # es = WarpUtils.transposeLag(efrs[2])
    # WarpUtils.normalizeErrors(es,-Float.MAX_VALUE,1E10)
    # plot2(es,title="Smoothed i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.5],width=900,height=600,o=x1rx2u)
    # esks = WarpUtils.transposeLag(efrsks[2])
    # WarpUtils.normalizeErrors(esks,-Float.MAX_VALUE,1E10)
    # plot2(esks,title="Smoothed i2=%d"%i2,s2=su,
    #       hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Distance",
    #       clips1=[0,0.5],width=900,height=600,o=x1rx2u)
    x1u,x2u = getShiftCoordsI2(g1,i2,u)
    x1u,x2u = mul(x1u,d1),mul(x2u,d1)
    x1uks,x2uks= getShiftCoordsI2(g1,i2,uks)
    x1uks,x2uks= mul(x1uks,d1),mul(x2uks,d1)
    # x1ugks,x2ugks= getShiftCoords1(g1ks,ugks)
    # x1ugks,x2ugks= mul(x1ugks,d1),mul(x2ugks,d1)
    x12SingleA = [[x1u,x2u,"coarse grid","wO",10.0],
                 [x1uks,x2uks,"coarse grid ks","cO",10.0],
                 # [x1ugks,x2ugks,"coarse grid gks","mO",10.0],
                 [x1ks,x2ks,"known shifts","yO",10.0]]
    uA = [[mul(u,d1),"u","w--.",2.0],[mul(uks,d1),"u with ks","c--",2.0]]
          # [mul(ugks,d1),"u with gks","m-.",2.0]]
    plot2(et,pA=uA,x12SingleA=x12SingleA,title="AE i2=%d"%i2,s1=se,s2=su,
          hLabel="PP time (s)",vLabel="Time shift (s)",cbar="Error",
          clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
          o=x1rx2u)

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

def toFloats1(f):
  n = len(f)
  g = zerofloat(n)
  for i in range(n):
    g[i] = float(f[i])
  return g

def getG1Floats(x1m,g1Flat,d1):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2)
  for i2 in range(n2):
    for i1 in range(ng):
      g1[i2][i1] = x1m[i2][g1Flat[i1]]*d1
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

def getG2(dg2):
  g2 = Subsample.subsample(n2,dg2); # simple regular interval grid
  ng2 = len(g2)  
  g2f = zerofloat(ng2)
  for i2 in range(ng2):
    g2f[i2] = s2.getValue(g2[i2])
  return g2f
  
def getShiftCoords(sf,g1,g2,u):
  ng1 = len(g1[0])
  ng2 = len(g2)
  x1 = fillfloat(-10,ng1,n2)
  x2 = fillfloat(-10,ng1,n2)
  for i2 in range(ng2):
    g2i = s2.indexOfNearest(g2[i2])
    for i1 in range(ng1):
      x1[g2i][i1] = g1[g2i][i1]
      x2[g2i][i1] = u[g2i][sf.indexOfNearest(g1[g2i][i1])]
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
    g2i = s2.indexOfNearest(g2[i2])
    for i1 in range(ng1):
      x1[i2][i1] = g1[g2i][i1];
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

def setPlotVars(dw,nf1,sf1):
  global ne1,nu,se,su,el,ul,xl,n1,s1
  n1,s1 = nf1,sf1
  se = dw.getSampling1() # error sampling
  su = dw.getSamplingU() # shift sampling
  ne1 = se.getCount()
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting

def plotGrid(f,x1i,x2i,fname):
  # Plot f with grid points
  # x1i = mul(x1i,d1)
  plot2(f,x12=x1i,x22=x2i,title=fname,s1=s1,s2=s2,vLabel=fname+" time (s)",
        hLimits=xl,vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotShifts(f,u,desc,fname):
  # Plot 2D shifts
  zm = ZeroMask(f)
  zm.apply(0.00,u)
  plot2(u,title="Vertical shifts"+desc,s1=s1,s2=s2,vLabel=fname+" time (s)",
        cbar="Vertical shift (s)",cmap1=jet,clips1=uClips,vLimits=el,
        cbw=100)

def plotVpvs(f,u,desc,fname):
  # Plot Vp/Vs ratios
  zm = ZeroMask(f)
  vpvs = WarpUtils.vpvs(su,u)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  plot2(f,g=vpvs,title="Vp/Vs"+desc,s1=s1,s2=s2,vLabel=fname+" time (s)",
        cbar="Vp/Vs ratio",clips1=iClips,cmap2=jet,clips2=vClips,hLimits=xl,
        vLimits=el,cmap1=iMap,cbw=100)

def plotWarped(sf,f,sg,g,u,desc,fname,gname):
  # Plot Warped g image and the difference between the warped and input.
  h = WarpUtils.applyShifts(sf,u,sg,g)
  plot2(h,title=gname+" warped"+desc,s1=s1,s2=s2,vLabel=fname+" time (s)",
        vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)
  nrms = WarpUtils.computeNrms(ne1,f,h)
  print desc+": nrms=%g"%nrms
  # d = sub(h,f)
  # plot2(d,title=fname+"-"+gname+" warped"+desc,s1=s1,s2=s2,
  #       vLabel=fname+" time (s)",vLimits=el,clips1=iClips,cmap1=iMap,cbw=100)

def plotErrors(sf,f,sg,g,dw,uA,desc,fname,x12SliceA=None):
  # Plot Alignment Errors with coarse grid points and shifts
  e = dw.computeErrors2(sf,f,sg,g)
  e = WarpUtils.transposeLag12(e)
  WarpUtils.normalizeErrors(e)
  plot3(e,pA=uA,x12SliceA=x12SliceA,title="AE"+desc,s1=se,s2=su,
        hLabel=fname+" time (s)",vLabel="Vertical shift (s)",cbar="Error",
        clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
        o=x1rx2u)
  # plot3(e,title="AE"+desc,s1=se,s2=su,
  #       hLabel=fname+" time (s)",vLabel="Vertical shift (s)",cbar="Error",
  #       clips1=[0,0.15],width=900,height=600,hLimits=el,vLimits=ul,
  #       o=x1rx2u)

###############################################################################
run(main)
