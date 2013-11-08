#############################################################################
# Frequency Analysis/Filtering

from gbcUtils import *

#############################################################################
s1,s2,s3 = getSamplings()
baseDir = getBaseDir()
iClips=[-1.0,1.0]
mClips=[0.0,1.0]
# klower,kupper,kwidth,aerror = 0.0,0.16,0.1,0.01
klower,kupper,kwidth,aerror = 0.0,0.06,0.05,0.01

#############################################################################

def main(args):
  goFilterDesign("ps1",i3=72)
  #goFilterDesign("ps2",i3=72)
  #goFilter("ps2",0.2)
  #goFilterDesign("ps1")
  #viewFkk()

def viewFkk():
  ps1 = getGbcImage(baseDir,"ps1")
  ps1fkk = getGbcImage(baseDir,"ps1_fkk")
  ps1fkk2 = getGbcImage(baseDir,"ps1_fkk2")
  he0 = 320
  plotPP3(ps1,s1=s1,s2=s2,s3=s3,title="PS1",label1="PS time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100,he0=he0)
  plotPP3(ps1fkk,s1=s1,s2=s2,s3=s3,title="PS1 fkk agressive",
          label1="PS time (s)",label2="Crossline (km)",label3="Inline (km)",
          cbar="Amplitude",clips1=iClips,width=900,height=1000,cbw=100,he0=he0)
  plotPP3(ps1fkk2,s1=s1,s2=s2,s3=s3,title="PS1 fkk",label1="PS time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100,he0=he0)
  go3dAmpSpec(ps1,"PS1")
  go3dAmpSpec(ps1fkk,"PS1 fkk aggressive")
  go3dAmpSpec(ps1fkk2,"PS1 fkk")

def goFilterDesign(name,i3=None):
  f = getGbcImage(baseDir,name)
  if i3:
    f = f[i3]
    zm = ZeroMask(f)
    fk = FK(s1,s2,f)
    fk.setZeroMask(zm)
    fk.designFilter()
  else:
    zm = ZeroMask(f)
    fkk = FKK(s1,s2,s3,f)
    fkk.setZeroMask(zm)

def goFilter(name,m):
  f = getGbcImage(baseDir,name)
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  g = zerofloat(n1,n2,n3)
  for i3 in range(len(f)):
    fi3 = f[i3]
    zm = ZeroMask(fi3)
    fk = FK(s1,s2,fi3)
    fk.setZeroMask(zm)
    g[i3] = fk.inverse(fk.filter(m))
  plotPP3(f,s1=s1,s2=s2,s3=s3,title="Input",label1="Time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100,he0=320)
  plotPP3(g,s1=s1,s2=s2,s3=s3,title="Filtered",label1="Time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100,he0=320)
  writeImage(baseDir,name+"_fk",g)

def go3DPP():
  pp,zmpp,ps1,zmps1 = getGbcImages(baseDir)
  plotPP3(pp,s1=s1,s2=s2,s3=s3,title="PP",label1="PP time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100)
  # plotPP3(zmpp.getAsFloats3(),s1=s1,s2=s2,s3=s3,title="PP Mask",
  #         label1="PP time (s)",label2="Crossline (km)",label3="Inline (km)",
  #         cbar="Amplitude",clips1=mClips,width=900,height=1000,cbw=100,
  #         lineColor=Color.YELLOW)
  goFft3(pp,desc=" PP")
  # ppf = BandPass.goBandPass(klower,kupper,kwidth,aerror,pp)
  # zmpp.apply(0.0,ppf)
  # plotPP3(ppf,s1=s1,s2=s2,s3=s3,title="PP BPF",label1="PP time (s)",
  #         label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
  #         clips1=iClips,width=900,height=1000,cbw=100)
  # goFft3(ppf,desc=" PP BPF")
  # writeImage(baseDir,"pp_bpf",ppf)
  
def go3DPS1():
  pp,zmpp,ps1,zmps1 = getGbcImages(baseDir)
  plotPP3(ps1,s1=s1,s2=s2,s3=s3,title="PS1",label1="PS time (s)",
          label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
          clips1=iClips,width=900,height=1000,cbw=100)
  # plotPP3(zmps1.getAsFloats3(),s1=s1,s2=s2,s3=s3,title="PS1 Mask",
  #         label1="PS time (s)",label2="Crossline (km)",label3="Inline (km)",
  #         cbar="Amplitude",clips1=mClips,width=900,height=1000,cbw=100,
  #         lineColor=Color.YELLOW)
  goFft3(ps1,desc=" PS1")
  ps1f = BandPass.goBandPass(klower,kupper,kwidth,aerror,pp)
  # zmps1.apply(0.0,ps1f)
  # plotPP3(ps1f,s1=s1,s2=s2,s3=s3,title="PS1 BPF",label1="PS time (s)",
  #         label2="Crossline (km)",label3="Inline (km)",cbar="Amplitude",
  #         clips1=iClips,width=900,height=1000,cbw=100)
  # goFft3(ps1f,desc=" PS1 BPF")
  # writeImage(baseDir,"ps1_bpf",ps1f)

def filter(f,desc):
  g = goFft2(f,desc=desc)
  h = sub(f,g)
  plot2(f,s1=s1,s2=s2,title=desc,vLabel="PS time (s)",
        hLabel="Crossline (km)",cbar="Amplitude",clips1=iClips,cbw=100)
  plot2(g,s1=s1,s2=s2,title=desc+" removed",vLabel="PS time (s)",
        hLabel="Crossline (km)",cbar="Amplitude",clips1=iClips,cbw=100)
  plot2(h,s1=s1,s2=s2,title=desc+" filtered",vLabel="PS time (s)",
        hLabel="Crossline (km)",cbar="Amplitude",clips1=iClips,cbw=100)
  
def goFft3(f,desc=""):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  fft = Fft(n1,n2,n3)
  fft.setPadding1(n1/2)
  s1 = fft.getFrequencySampling1()
  s2 = fft.getFrequencySampling2()
  s3 = fft.getFrequencySampling3()
  nk1 = s1.getCount()
  nk2 = s2.getCount()
  nk3 = s3.getCount()
  # print "n1=%d, nk1=%d"%(n1,nk1)
  # print "n2=%d, nk2=%d"%(n2,nk2)
  # print "n3=%d, nk3=%d"%(n3,nk3)
  ft = fft.applyForward(f)
  a = zerofloat(nk1)
  for i3 in range(nk3):
    for i2 in range(nk2):
      a = add(a,mul(20.0,log10(cabs(ft[i3][i2]))))
  plot1(a,s=s1,title="Amplitude"+desc,hLabel="Frequency (cycles/sample)",
        vLabel="Amplitude")
        # vLimits=[-4000,8000])
  a3 = mul(20.0,log10(cabs(ft)))
  plotPP3(a3,s1=s1,s2=s2,s3=s3,title="FKK"+desc,
          label1="Frequency (cycles/sample)",label2="K2",label3="K3",
          cbar="Amplitude",width=900,height=1000,cbw=100,cmap1=jet)

def go3dAmpSpec(f,desc=""):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  fft = Fft(n1,n2,n3)
  s1 = fft.getFrequencySampling1()
  s2 = fft.getFrequencySampling2()
  s2 = fft.getFrequencySampling3()
  nk1 = s1.getCount()
  nk2 = s2.getCount()
  nk3 = s3.getCount()
  print "n1=%d, nk1=%d"%(n1,nk1)
  print "n2=%d, nk2=%d"%(n2,nk2)
  print "n3=%d, nk3=%d"%(n3,nk3)
  ft = fft.applyForward(f)
  a = zerofloat(nk1)
  for i3 in range(nk3):
    for i2 in range(nk2):
      a = add(a,mul(20.0,log10(cabs(ft[i3][i2]))))
  plot1(a,s=s1,title="Amplitude"+desc,hLabel="Frequency (cycles/sample)",
        vLabel="Amplitude")

def go2dAmpSpec(f,desc=""):
  n2 = len(f)
  n1 = len(f[0])
  fft = Fft(n1,n2)
  s1 = fft.getFrequencySampling1()
  s2 = fft.getFrequencySampling2()
  nk1 = s1.getCount()
  nk2 = s2.getCount()
  print "n1=%d, nk1=%d"%(n1,nk1)
  print "n2=%d, nk2=%d"%(n2,nk2)
  ft = fft.applyForward(f)
  a = zerofloat(nk1)
  for i2 in range(nk2):
    a = add(a,mul(20.0,log10(cabs(ft[i2]))))
  plot1(a,s=s1,title="Amplitude"+desc,hLabel="Frequency (cycles/sample)",
        vLabel="Amplitude")

def goFft2(f,desc=""):
  n2 = len(f)
  n1 = len(f[0])
  fft = Fft(n1,n2)
  # fft.setPadding1(n1)
  fft.setPadding2(n2)
  s1 = fft.getFrequencySampling1()
  s2 = fft.getFrequencySampling2()
  nk1 = s1.getCount()
  nk2 = s2.getCount()
  # print "n1=%d, nk1=%d"%(n1,nk1)
  # print "n2=%d, nk2=%d"%(n2,nk2)
  ft = fft.applyForward(f)
  ftt = copy(ft)
  # re-arrange
  i2h = s2.getValue(nk2/2)
  h2 = nk2/2
  copy(nk1,h2,0,h2,ftt,0,0,ft)
  copy(nk1,h2,0,0,ftt,0,h2,ft)
  s2 = Sampling(nk2,s2.getDelta(),-0.5)
  ftz = copy(ft)
  for i in range(40,nk2-1-40):
  # for i in range(40):
    for j in range(nk1):
      ftz[i][j] = ftz[i][j]*0.0
  # for i in range(40):
  #   for j in range(nk1):
  #     ftz[nk2-1-i][j] = ftz[nk2-1-i][j]*0.0
  fi = fft.applyInverse(ftz)
  a2 = cabs(ft)
  a2i = cabs(ftz)
  plot2(a2,s1=s1,s2=s2,title="FK"+desc,vLabel="Frequency (cycles/sample)",
        hLabel="K2",cbar="Amplitude",cbw=100,cmap1=jet)
  plot2(a2i,s1=s1,s2=s2,title="FK cut"+desc,vLabel="Frequency (cycles/sample)",
        hLabel="K2",cbar="Amplitude",cbw=100,cmap1=jet)
  return fi

#############################################################################
run(main)
