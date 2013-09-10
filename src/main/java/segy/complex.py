from imports import *

################################################################################

nf = 501
df = 0.5/(nf-1)
ff = 0.0
sf = Sampling(nf,df,ff)

af = []

def main(args):
  c = Cdouble(1.0,0.0)
  zinv,czinv,cmcz,pte,poq,apoq = doTest(c)
  plotZinv(zinv,czinv,cmcz,pte,poq,apoq)

def doTest(c):
  zinvList,cZinv,coneMcz,pte,poq,apoq = [],[],[],[],[],[]
  for i in range(0,nf):
    f = ff+i*df
    zinv = Cdouble.polar(1.0,-2.0*DBL_PI*f)
    zinvList.append(zinv)
    cz = c.times(zinv)
    cZinv.append(cz)
    cone = Cdouble(1.0,0.0)
    cmcz = cone.minus(cz)
    coneMcz.append(cmcz)
    p = Cdouble(1.0,0.0)
    p.timesEquals(cmcz)
    pte.append(p)
    q = Cdouble(1.0,0.0)
    h = p.over(q)
    poq.append(h)
    apoq.append(h.abs())
  return zinvList,cZinv,coneMcz,pte,poq,apoq

def plotZinv(zinvList,cZinv,coneMcz,pte,poq,apoq):
  rz,iz,rcz,icz,rcmcz,icmcz,rpte,ipte,rpoq,ipoq = [],[],[],[],[],[],[],[],[],[]
  for j in range(0,len(zinvList)):
    rz.append(zinvList[j].r)
    iz.append(zinvList[j].i)
    rcz.append(cZinv[j].r)
    icz.append(cZinv[j].i)
    rcmcz.append(coneMcz[j].r)
    icmcz.append(coneMcz[j].i)
    rpte.append(pte[j].r)
    ipte.append(pte[j].i)
    rpoq.append(poq[j].r)
    ipoq.append(poq[j].i)
  pp = PlotPanel(5,2)
  pp.addPoints(0,0,rz,iz)
  pp.addPoints(1,0,rcz,icz)
  pp.addPoints(2,0,rcmcz,icmcz)
  pp.addPoints(3,0,rpte,ipte)
  pp.addPoints(4,0,rpoq,ipoq)
  a = Real1(sf, apoq)
  pp.addPoints(0,1,a.getSampling(),a.getValues())
  pf = PlotFrame(pp)
  pf.setSize(500,500)
  pf.setVisible(True)

################################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
