"""
Reads and reformats horizons from GBC.
"""
from imports import *
from gbcUtils import *

# Horizon names with colors, ordered by increasing time.
# horizonColors = {
#   "a_time":[Color.RED,Color.BLUE]
# }
ppHorizonColors = {
  "a_pp_time":[Color.RED,Color.WHITE],
  "m_pp_time":[Color.GREEN,Color.WHITE],
  "b_pp_time":[Color.ORANGE,Color.WHITE],
  "c_pp_time":[Color.MAGENTA,Color.WHITE],
  "d_pp_time":[Color.CYAN,Color.WHITE],
  "n_pp_time":[Color.YELLOW,Color.WHITE],
  "e_pp_time":[Color.PINK,Color.WHITE],
  "f_pp_time":[Color.BLUE,Color.WHITE],
  "l_pp_time":[Color.RED,Color.WHITE],
  "g_pp_time":[Color.GREEN,Color.WHITE],
  "h_pp_time":[Color.ORANGE,Color.WHITE],
  "i_pp_time":[Color.MAGENTA,Color.WHITE],
  "j_pp_time":[Color.CYAN,Color.WHITE],
  "k_pp_time":[Color.YELLOW,Color.WHITE]
}
ppHorizonNames = ppHorizonColors.keys()

ps1HorizonColors = {
  "a_ps1_time":[Color.RED,Color.WHITE],
  "b_ps1_time":[Color.ORANGE,Color.WHITE],
  "c_ps1_time":[Color.MAGENTA,Color.WHITE],
  "d_ps1_time":[Color.CYAN,Color.WHITE],
  "e_ps1_time":[Color.PINK,Color.WHITE],
  "f_ps1_time":[Color.BLUE,Color.WHITE],
  "g_ps1_time":[Color.GREEN,Color.WHITE],
  "h_ps1_time":[Color.ORANGE,Color.WHITE],
  "i_ps1_time":[Color.MAGENTA,Color.WHITE]
}
ps1HorizonNames = ps1HorizonColors.keys()

##############################################################################
def main(args):
  # setGlobals(True); makeHorizons() # make binary files for pp time horizons
  # setGlobals(True); viewHorizons() # view pp time horizons
  # setGlobals(False); makeHorizons() # make binary files for ps1 time horizons
  setGlobals(False); viewHorizons() # view ps1 time horizons

def setGlobals(doPP):
  global horizonDir,horizonNames,horizonColors,datName,s1,s2,s3
  gbcDir = "/data/seis/gbc/"
  horizonDir = gbcDir+"horizons/"
  s1,s2,s3 = getSamplings() # Samplings defined in gbcUtils
  if doPP:
    horizonNames = ppHorizonNames
    horizonColors = ppHorizonColors
    datName = "pp_smooth"
  else:
    horizonNames = ps1HorizonNames
    horizonColors = ps1HorizonColors
    datName = "ps1_fkk_smooth"

def makeHorizons():
  for name in horizonNames:
    h = Horizon.readText(horizonDir+name+".xyz")
    print name," ns =",h.ns," nt=",h.nt
    h.writeBinary(horizonFile(name))

def viewHorizons():
  x = getGbcImage(getBaseDir(),datName)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  world = World()
  world.addChild(ipg)
  addHorizonGroups(world)
  frame = makeFrame(world)

def horizonFile(name):
  return horizonDir+name+".dat"

def horizonRead(name):
  return Horizon.readBinary(horizonFile(name))

def addHorizonGroups(world):
  for hname in horizonNames:
    h = horizonRead(hname)
    colors = horizonColors[hname]
    ct = colors[0] # triangle colors
    cl = colors[1] # line colors
    tg = makeTriangleGroup(h,ct)
    world.addChild(tg)
    lgx3 = makeLineGroup(h,cl,72,True)
    world.addChild(lgx3)
    lgx2 = makeLineGroup(h,cl,75,False)
    world.addChild(lgx2)

def makeFrame(world):
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,2.0)
  view.setScale(1.5)
  view.setAzimuth(-50.0)
  frame.setSize(1200,900)
  frame.setVisible(True)
  return frame

def readImage(filename):
  n1,n2,n3 = s1.count,s2.count,s3.count
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(filename)
  ais.readFloats(x)
  ais.close()
  return x

def makeTriangleGroup(horizon,color):
  ijk = horizon.getIABC()
  xyz = horizon.getX321()
  tg = TriangleGroup(ijk,xyz)
  tg.setColor(color)
  return tg;

def makeLineGroup(horizon,color,index,x3=True):
  xyz = None
  if x3:
    xyz = horizon.getI3X321(s2,s3,index)
  else:
    xyz = horizon.getI2X321(s2,s3,index)
  rgb = makeLineColor(len(xyz),color)
  lg = LineGroup(xyz,rgb)
  return lg;

def makeLineColor(n,color):
  rgb = zerofloat(n)
  r,g,b = color.getRed(),color.getGreen(),color.getBlue()
  j = 0;
  for i in range(n/3):
    rgb[j  ] = r
    rgb[j+1] = g
    rgb[j+2] = b
    j = j+3
  return rgb
    
#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
