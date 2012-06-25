"""
Converts sgy image files (SEG-Y format) to dat files (3D arrays of floats)
Removes all SEG-Y headers and, if necessary, converts the data format to
IEEE floats. The byte order for floats in dat files is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2012.06.18
"""
from imports import *
global n1,n2,n3

#############################################################################
def main(args):
  #goMbs()
  #goNorne()
  #goSino()
  goF3d()

def goMbs():
  """
  ***************************************************************************
  PstmSmall (pstm_fraw.sgy):
  number of bytes = 959341200
  number of traces = 153740
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  n1 =  1500 (number of samples per trace)
  n2 =   448 (number of traces in inline direction)
  n3 =   422 (number of traces in crossline direction)
  d1 = 0.002000 (time sampling interval, in s)
  d2 = 0.016764 (inline sampling interval, in km)
  d3 = 0.016764 (crossline sampling interval, in km)
  i2min =   601, i2max =  1048 (inline index bounds)
  i3min =  1001, i3max =  1422 (crossline index bounds)
  xmin = 823.725048, xmax = 831.574867 (x coordinate bounds, in km)
  ymin = 245.803217, ymax = 255.588516 (y coordinate bounds, in km)
  grid azimuth =  58.12 degrees
  ***************************************************************************
  PstmLarge (pstm_raw_cut.sgy):
  number of bytes = 4647376320
  number of traces = 795783
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  n1 =  1400 (number of samples per trace)
  n2 =  1119 (number of traces in inline direction)
  n3 =  1189 (number of traces in crossline direction)
  d1 = 0.002000 (time sampling interval, in s)
  d2 = 0.016764 (inline sampling interval, in km)
  d3 = 0.016764 (crossline sampling interval, in km)
  i2min =   350, i2max =  1468 (inline index bounds)
  i3min =   234, i3max =  1422 (crossline index bounds)
  xmin = 822.376308, xmax = 842.769562 (x coordinate bounds, in km)
  ymin = 235.175146, ymax = 255.588516 (y coordinate bounds, in km)
  grid azimuth =  58.12 degrees
  good subset:
    i1min,i1max = 150, 650,  m1 = 501
    i2min,i2max = 490,1258,  m2 = 769
    i3min,i3max = 358, 917,  m3 = 560
  """
  imageType = "PstmLarge" # which image to process
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/seis/mbs/"
  if imageType=="PstmSmall":
    sgyfile = basedir+"PstmSmall/Marathon20070228/pstm_fraw.sgy"
    datfile = basedir+"dat/pstm_fraw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,601,1048,1001,1422
  elif imageType=="PstmLarge":
    sgyfile = basedir+"PstmLarge/pstm_raw_cut.sgy"
    datfile = basedir+"dat/pstm_raw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,490,1258,358,917
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    f = readImage(datfile,n1,n2,n3)
    show3d(f,clip=10000.0)

def goNorne():
  """
  ***************************************************************************
  norne4d_2006-full.sgy:
  number of bytes = 1363689924
  number of traces = 321321
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  n1 =  1001 (number of samples per trace)
  n2 =  1001 (number of traces in inline direction)
  n3 =   321 (number of traces in crossline direction)
  d1 = 0.004000 (time sampling interval, in s)
  d2 = 0.012501 (inline sampling interval, in km)
  d3 = 0.012500 (crossline sampling interval, in km)
  i2min =  1300, i2max =  2300 (inline index bounds)
  i3min =   970, i3max =  1290 (crossline index bounds)
  xmin =  453.320000, xmax =  464.634000 (x coordinate bounds, in km)
  ymin = 7317.354880, ymax = 7329.340160 (y coordinate bounds, in km)
  grid corner points:
    i2min =  1300, i3min =   970, x =  453.320000, y = 7320.021120
    i2max =  2300, i3min =   970, x =  461.652000, y = 7329.340160
    i2min =  1300, i3max =  1290, x =  456.302000, y = 7317.354880
    i2max =  2300, i3max =  1290, x =  464.634000, y = 7326.673920
  grid azimuth: 41.80 degrees
  ***************************************************************************
  Full Norne corner coordinates (Knut, 16/7-09)
  line	trace	X	Y
  970	2300	461652	7329340.16
  970	1300	453320	7320021.12
  1290	1300	456302	7317355
  1290	2300	464634	7326674
  """
  firstLook = True # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/seis/norne/"
  sgyfile = basedir+"sgy/norne4d_2006-full.sgy"
  datfile = basedir+"dat/full2006.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,1000,1300,2300,970,1290
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1000.0)

def goSino():
  """
  Two 2D (x- and z-component) images
  nt = 2001, dt = 0.004 s
  nx =  721, dx = 0.030 km?
  samples [0:1250] in z component ~ samples [0:2000] in x component
  """
  basedir = "/data/seis/sino/"
  for component in ["x","z"]:
    sgyfile = basedir+"260_"+component+"_201-921_stack.segy"
    datfile = basedir+component+"260.dat"
    i1min,i1max,i2min,i2max = 0,2000,201,921
    n1,n2 = 1+i1max-i1min,1+i2max-i2min
    si = SegyImage(sgyfile,ByteOrder.LITTLE_ENDIAN) # non-standard byte order!
    #si.printSummaryInfo()
    #si.printBinaryHeader()
    #si.printTraceHeader(0)
    #si.printTraceHeader(1)
    #plotIbmIeeeFloats(si)
    si.setD2(0.015) # a guess, from looking at group coordinates in headers
    if component=="x":
      si.setFormat(5) # formats appear to be IEEE for x and IBM for z!???
    si.printAllInfo()
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,0,0)
    si.close()
    f = readImage(datfile,n1,n2)
    gain(500,f)
    if component=="z":
      stretch(1.6,f)
    show2d(f,title=component+" component")

def goF3d():
  """
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/f3d/f3draw.sgy
  byte order = BIG_ENDIAN
  number of bytes = 699003060
  number of traces = 600515
  format = 3 (2-byte two's complement integer)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   300, i2max =  1250 (inline indices)
    i3min =   100, i3max =   750 (crossline indices)
    xmin =  605.416700, xmax =  629.576300 (x coordinates, in km)
    ymin = 6073.556400, ymax = 6090.463200 (y coordinates, in km)
  grid sampling:
    n1 =   462 (number of samples per trace)
    n2 =   951 (number of traces in inline direction)
    n3 =   651 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.025000 (inline sampling interval, in km)
    d3 = 0.024999 (crossline sampling interval, in km)
  grid corner points:
    i2min =   300, i3min =   100, x =  605.835500, y = 6073.556400
    i2max =  1250, i3min =   100, x =  629.576300, y = 6074.219900
    i2min =   300, i3max =   750, x =  605.381800, y = 6089.799700
    i2max =  1250, i3max =   750, x =  629.122600, y = 6090.463200
  grid azimuth: 88.40 degrees
  ****** end of SEG-Y file info ******
  good subset
  i1min,i1max,i2min,i2max,i3min,i3max = 0,461,300,1250,100,690
  n1,n2,n3 = 462,951,591
  """
  firstLook = True # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/Users/scompton/opendtect/F3_Demo/Seismics/"
  sgyfile = basedir+"f3d.sgy"
  datfile = basedir+"f3draw.dat"
  datfileWrong = basedir+"f3dieee.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,461,300,1250,100,690
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    #si.writeFloats(datfile,i1min,i1max,i2min,i2max,i3min,i3max)
    si.setFormat(5)
    si.writeFloats(datfileWrong,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    y = readImage(datfileWrong,n1,n2,n3)
    show3d(x,clip=10000.0)
    show3d(y,clip=1000.0)

def show2d(f,clip=None,title=None):
  print "show2d: f min =",min(f)," max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if clip:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(2,98)
  if title:
    sp.setTitle(title)
  sp.setSize(600,1100)
  
def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)

def plot23(si):
  i2 = si.getI2sAsFloats()
  i3 = si.getI3sAsFloats()
  sp = SimplePlot()
  sp.setHLabel("inline sample index i2")
  sp.setVLabel("crossline sample index i3")
  pv = sp.addPoints(i2,i3)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(i2,i3)
  sp.setSize(w,h)

def plotXY(si):
  x = si.getXs()
  y = si.getYs()
  sp = SimplePlot()
  sp.setHLabel("x coordinate (km)")
  sp.setVLabel("y coordinate (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(x,y)
  sp.setSize(w,h)

def goodWidthHeight(x,y):
  xmin,xmax = min(x),max(x)
  ymin,ymax = min(y),max(y)
  w,h = 1000,1000
  if (xmax-xmin)>(ymax-ymin):
    h = int(h*(ymax-ymin)/(xmax-xmin))
  else:
    w = int(w*(xmax-xmin)/(ymax-ymin))
  return w,h

def readImage(datfile,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  aos = ArrayOutputStream(datfile)
  aos.writeFloats(x)
  aos.close()

def plotIbmIeeeFloats(si):
  ntrace = si.countTraces()
  itrace = ntrace/2
  fmt = si.getFormat()
  si.setFormat(1) # IBM floats
  fibm = si.getTrace(itrace)
  si.setFormat(5) # IEEE floats
  fieee = si.getTrace(itrace)
  si.setFormat(fmt)
  pp = PlotPanel(2,1)
  pp.setTitle("IBM (top) versus IEEE (bottom)")
  pp.setHLabel(0,"Sample index")
  pp.setVLabel(0,"IBM amplitude")
  pp.setVLabel(1,"IEEE amplitude")
  pp.addPoints(0,0,fibm)
  pp.addPoints(1,0,fieee)
  pf = PlotFrame(pp)
  pf.setSize(1000,800)
  pf.setVisible(True)

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0.0,1.0/c,n1)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  g = zerofloat(n1)
  for i2 in range(n2):
    si.setUniformSamples(f[i2])
    si.interpolate(n1,1.0/c,0.0,g)
    copy(g,f[i2])

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

#############################################################################
# Parihaka functions below are outdated; kept here for updating later

def goParihaka():
  """
  INLINES  : 1665 - 5599 (INC 1)   CROSSLINES : 2050 - 14026 (INC 2)
  TIME: 1501 samples  (6 sec, interval 4ms)
  nbytes = 69,479,807,644
  ntrace = (nbytes-nhead-nbhed)/(240+4*n1)
  """
  fmt = 1
  n1 = 1501
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"Parihaka3d_raw.sgy"
  datfile = datdir+"par11.dat"
  #displayParihaka(datfile)
  #makeSubsetsParihaka(datdir)
  #displaySubsetParihaka(datfile,0,500,500,751,501,501)
  #bigSubsetParihaka(n1,sgyfile,datfile)
  #convert(n1,n2,n3,sgyfile,datfile)

def displayParihaka(datfile):
  n1,n2,n3 = 751,1001,1001
  x = readImage(datfile,n1,n2,n3)
  display3d(x,1.0e5)

def makeSubsetsParihaka(datdir):
  m1,m2,m3 = 751,1001,1001
  x = readParihaka(datdir+"parBig.dat",0,0,0,m1,m2,m3)
  writeImage(datdir+"par00.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,0,m1,m2,m3)
  writeImage(datdir+"par01.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,0,1000,m1,m2,m3)
  writeImage(datdir+"par10.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,1000,m1,m2,m3)
  writeImage(datdir+"par11.dat",x)
  display3d(x,1.0e5)

def displaySubsetParihaka(datfile,j1,j2,j3,m1,m2,m3):
  x = readSubsetParihaka(datfile,j1,j2,j3,m1,m2,m3)
  display3d(x,1.0e5)

def readParihaka(datfile,j1,j2,j3,m1,m2,m3):
  n1,n2,n3 = 1501,2001,2001
  m1 = min(m1,n1)
  m2 = min(m2,n2)
  m3 = min(m3,n3)
  j1 = min(j1,n1-m1)
  j2 = min(j2,n2-m2)
  j3 = min(j3,n3-m3)
  x = zerofloat(m1,m2,m3)
  ais = ArrayInputStream(datfile)
  for i3 in range(j3):
    ais.skipBytes(4*n1*n2)
  ais.skipBytes(4*(j1+n1*j2))
  for i3 in range(m3):
    if i3%10==0: 
      print "i3 =",i3
    for i2 in range(m2):
      ais.readFloats(x[i3][i2])
      ais.skipBytes(4*(n1-m1))
    ais.skipBytes(4*n1*(n2-m2))
  return x

def bigSubsetParihaka(n1,sgyfile,datfile):
  """ big subset 1501x2001x2001 ~ 24 GB
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  """
  """ A: small subset 501x501x501 ~ 500 MB
  i1min,i1max =    0, 500 # time samples
  i2min,i2max = 7500,8500 # increment by 2
  i3min,i3max = 2500,3000 # increment by 1
  """
  """ B: subset 751x1001x1001 ~ 500 MB
  i1min,i1max =    0, 750 # time samples
  i2min,i2max = 7000,9000 # increment by 2
  i3min,i3max = 3000,4000 # increment by 1
  """
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  m1 = 1+i1max-i1min
  m2 = 1+(i2max-i2min)/2
  m3 = 1+i3max-i3min
  m23 = m2*m3
  ais = ArrayInputStream(sgyfile,bo)
  aos = ArrayOutputStream(datfile)
  ais.skipBytes(nhead)
  ais.skipBytes(nbhed)
  h = zeroint(nthed/4)
  x = zeroint(n1)
  y = zeroint(m1)
  z = zerofloat(m1)
  nread = 0
  i23 = 0
  while i23<m23:
    nread += 1
    ais.readInts(h)
    i3,i2 = h[49],h[50]
    if i2min<=i2 and i2<=i2max and i3min<=i3 and i3<=i3max:
      i23 += 1
      if i2==i2min:
        print "nread =",nread," i3min =",i3min," i3 =",i3," i3max =",i3max
      ais.readInts(x) # read trace samples
      copy(m1,x,y) # keep only m1 samples
      IbmIeee.ibmToFloat(y,z) # ibm to ieee
      aos.writeFloats(z) # write trace samples
    else:
      ais.skipBytes(4*n1)
  ais.close()
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 