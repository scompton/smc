###############################################################################
# Structure Tensors and Structure Oriented Smoothing for GBC

from gbcUtils import *
###############################################################################

baseDir = getBaseDir()
s1,s2,s3 = getSamplings()
iClips = [-1.0,1.0]

sigma1,sigma2,sigma3 = 8.0,8.0,8.0
c = 4.0
###############################################################################

def main(args):
  goPP()
  goPS1()

def goPP():
  pp = getGbcImage(baseDir,"pp")
  zm = ZeroMask(pp)
  # goTensors(pp,"pp_tensors")
  # display(pp,name="pp_tensors")
  pps = goSmooth(pp,"pp_tensors")
  DynamicWarpingC.normalize(pps,-1.5,1.5)
  zm.apply(0.0,pps)
  plotPP3(pp,s1=s1,s2=s2,s3=s3,title="PP",label1="PP time (s)",clips1=iClips,
          cbw=100,limits1=[0.0,3.0])
  plotPP3(pps,s1=s1,s2=s2,s3=s3,title="PP Smooth",label1="PP time (s)",
          clips1=iClips,cbw=100,limits1=[0.0,3.0])
  writeImage(baseDir,"pp_smooth",pps)
  # showTwo(pp,pps)
  
def goPS1():
  ps1 = getGbcImage(baseDir,"ps1_fkk")
  zm = ZeroMask(ps1)
  # goTensors(ps1,"ps1_fkk_tensors")
  # display(ps1,name="ps1_fkk_tensors")
  ps1s = goSmooth(ps1,"ps1_fkk_tensors")
  DynamicWarpingC.normalize(ps1s,-1.5,1.5)
  zm.apply(0.0,ps1s)
  plotPP3(ps1,s1=s1,s2=s2,s3=s3,title="PS1",label1="PS time (s)",clips1=iClips,
          cbw=100,limits1=[0.0,3.0])
  plotPP3(ps1s,s1=s1,s2=s2,s3=s3,title="PS1 Smooth",label1="PP time (s)",
          clips1=iClips,cbw=100,limits1=[0.0,3.0])
  writeImage(baseDir,"ps1_fkk_smooth",ps1s)
  # showTwo(ps1,ps1s)
  
def goTensors(f,name):
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  et3 = lof.applyForTensors(f)
  et3.invertStructure(0.0,2.0,4.0)
  writeTensors(name,et3)

def goSmooth(f,name):
  et3 = readTensors(name)
  lsf = LocalSmoothingFilter()
  g = copy(f)
  lsf.apply(c,f,g)
  return g
  
def display(f,g=None,name=None):
  world = World()
  ipg = ImagePanelGroup(s1,s2,s3,f)
  ipg.setClips(-1.0,1.0)
  world.addChild(ipg)
  if g:
    ipg2 = ImagePanelGroup(g)
    ipg2.setClips(-1.0,1.0)
    world.addChild(ipg2)
  # ipg = addImageToWorld(world,f)
  if name:
    et3 = readTensors(name)
    addTensorsInImage(ipg.getImagePanel(Axis.X),et3,30)
    addTensorsInImage(ipg.getImagePanel(Axis.Y),et3,30)
    addTensorsInImage(ipg.getImagePanel(Axis.Z),et3,30)
  frame = makeFrame(world)

def showTwo(g1,g2):
  sf = SimpleFrame()
  for g in [g1,g2]:
    ipg = sf.addImagePanels(s1,s2,s3,g)
    ipg.setClips(-1.0,1.0)
  sf.orbitView.setScale(1.0)
  sf.orbitView.setAxesScale(0.75,0.75,1.5)
  sf.setSize(1250,900)
  
def writeTensors(name,tensors):
  fos = FileOutputStream(baseDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  oos.close()

def readTensors(name):
  fis = FileInputStream(baseDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  ois.close()
  return tensors

###############################################################################
run(main)
