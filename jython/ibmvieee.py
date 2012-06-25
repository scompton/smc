from imports import *
import DoFft

################################################################################

def main(args):
  basedir = "/Users/scompton/opendtect/F3_Demo/Seismics/"
  sgyfile = basedir+"f3d.sgy"
  si = SegyImage(sgyfile)
  fibm,fieee = getTraces(si)
  plotIbmIeeeFloats(fibm,fieee)
  db = True;
  gibm,aibm = doFft(fibm,db)
  gieee,aieee = doFft(fieee,db)
  ribm = computeAvg(aibm.getValues(), "IBM")
  rieee = computeAvg(aieee.getValues(), "IEEE")
  plotFft(gibm,gieee,aibm,aieee,ribm,rieee)

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

def doFft(f, db):
  fft = DoFft(f, db)
  g = fft.getG()
  a = fft.getAmps()
  return g,a

def computeAvg(amplitudes, label):
  h = DoFft.weightedAvg(amplitudes, True)
  l = DoFft.weightedAvg(amplitudes, False)
  r = h/l
  print "%s: h=%f, l=%f, r=%f" % (label,h,l,r)
  return r

def plotFft(gibm, gieee, gkibm, gkieee, ribm, rieee):
  pp1 = PlotPanel(2,1)
  pp1.setTitle("Frequency - IBM (top) versus IEEE (bottom)")
  pp1.addPoints(0,0,gibm)
  pp1.addPoints(1,0,gieee)
  pf1 = PlotFrame(pp1)
  pf1.setSize(1000,800)
  pf1.setVisible(True)

  pp2 = PlotPanel(2,1)
  pp2.setTitle("IBM (top r=%f) versus IEEE (bottom r=%f)" % (ribm,rieee))
  pp2.setHLabel(0,"Frequency (cycles/sample)")
  pp2.setVLabel(0,"IBM amplitude (dB)")
  pp2.setVLabel(1,"IEEE amplitude (dB)")
  pp2.addPoints(0,0,gkibm.getSampling(),gkibm.getValues())
  pp2.addPoints(1,0,gkieee.getSampling(),gkieee.getValues())
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
