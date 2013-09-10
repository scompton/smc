from imports import *

################################################################################

def main(args):
  basedir = "/disk2/code/git/smc/"
  sgyfile = basedir+"f3d11lines.sgy"
  si = SegyImage(sgyfile)
  fibm,fieee = getTraces(si)
  sampling = IbmIeeeAmplitudes.getSampling(fibm)
  n1 = len(fibm)
  #plotIbmIeeeFloats(fibm,fieee)
  p=3.8
  n=0.2
  #testWeights(p,n)
  db = True;
  aibm = IbmIeeeAmplitudes.computeAmplitudes(fibm,db)
  aieee = IbmIeeeAmplitudes.computeAmplitudes(fieee,db)

  hw = zerofloat(n1)
  rhw = Real1(sampling, hw)
  lw = zerofloat(n1)
  rlw = Real1(sampling,
  ribm = IbmIeeeAmplitudes.computeLinearRatio(aibm.getValues(),hw,lw)
  rieee = IbmIeeeAmplitudes.computeLinearRatio(aieee.getValues(),None,None)
  plotFft(aibm,aieee,ribm,rieee,hw,lw)
  #ribm = IbmIeeeAmplitudes.computeSCRatio(aibm.getValues())
  #rieee = IbmIeeeAmplitudes.computeSCRatio(aieee.getValues())
  #ribm = IbmIeeeAmplitudes.computeRatio(aibm.getValues(),p,n)
  #rieee = IbmIeeeAmplitudes.computeRatio(aieee.getValues(),p,n)
  #plotFft(aibm,aieee,ribm,rieee)

def getTraces(si):
  ntrace = si.countTraces()
  itrace = ntrace/2
  fmt = si.getFormat()
  si.setFormat(1) # IBM floats
  fibm = si.getTrace(itrace)
  si.setFormat(5) # IEEE floats
  fieee = si.getTrace(itrace)
  si.setFormat(fmt)
  return fibm, fieee

def testWeights(p,n):
  #w1,w2 = [],[]
  w1,w2 = IbmIeeeAmplitudes.getWeights(p,n)
  #len = 501
  #for i in range(len):
    #d = (i/(len-1.0))*90.0
    #w1.append(cos(d*ArrayMath.DBL_PI/180.0))
    #w2.append(sin(d*ArrayMath.DBL_PI/180.0))
    #x = i/(len-1.0)
    #if (5*x < p):
      #w = (3*pow(5*x/p, 2) - 2*pow(5*x/p, 3)) * pow(x, n)
    #else:
      #w = pow(x, n)
    #print w,
    #w1.append(w)
  pp = PlotPanel()
  pp.addPoints(w1)
  pp.addPoints(w2)
  pf = PlotFrame(pp)
  pf.setVisible(True)

def plotFft(aibm, aieee, ribm, rieee, hw, lw):
  pp2 = PlotPanel(2,1)
  pp2.setTitle("IBM (top r=%f) versus IEEE (bottom r=%f)" % (ribm,rieee))
  pp2.setHLabel(0,"Frequency (cycles/sample)")
  pp2.setVLabel(0,"IBM amplitude (dB)")
  pp2.setVLabel(1,"IEEE amplitude (dB)")
  pp2.addPoints(0,0,aibm.getSampling(),aibm.getValues())
  pp2.addPoints(0,0,hw)
  pp2.addPoints(0,0,lw)
  pp2.addPoints(1,0,aieee.getSampling(),aieee.getValues())
  pp2.addPoints(1,0,hw)
  pp2.addPoints(1,0,lw)
  pf2 = PlotFrame(pp2)
  pf2.setSize(1000,800)
  pf2.setVisible(True)

def plotIbmIeeeFloats(fibm, fieee):
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

################################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
