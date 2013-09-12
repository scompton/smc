package warp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;

public class KnownShiftUtil {

  /**
   * 
   * @param su
   * @param s1
   */
  public KnownShiftUtil(Sampling su, Sampling s1) {
    this(su,s1,null,null,filldouble(0.0,s1.getCount()));
  }

  /**
   * 
   * @param su
   * @param s1
   */
  public KnownShiftUtil(Sampling su, Sampling s1, double[] r1Min) {
    this(su,s1,null,null,r1Min);
  }

  /**
   * 
   * @param su
   * @param s1
   * @param s2
   */
  public KnownShiftUtil(Sampling su, Sampling s1, Sampling s2) {
    this(su,s1,s2,null,filldouble(0.0,s1.getCount()));
  }

  /**
   * 
   * @param su
   * @param s1
   * @param s2
   */
  public KnownShiftUtil(Sampling su, Sampling s1, Sampling s2, double[] r1Min) {
    this(su,s1,s2,null,r1Min);
  }

  /**
   * 
   * @param su
   * @param s1
   * @param s2
   * @param s3
   */
  public KnownShiftUtil(
      Sampling su, Sampling s1, Sampling s2, Sampling s3) 
  {
    this(su,s1,s2,s3,filldouble(0.0,s1.getCount()));
  }

  /**
   * 
   * @param su
   * @param s1
   * @param s2
   * @param s3
   */
  public KnownShiftUtil(
      Sampling su, Sampling s1, Sampling s2, Sampling s3, double[] r1Min) 
  {
    _errorMap = new HashMap<>();
    _su = su;
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _r1Min = r1Min;
  }

  /**
   * 
   * @param xu
   * @param x1
   * @throws UnphysicalShiftException
   */
  public void add(double xu, double x1) throws UnphysicalShiftException {
    int iu = _su.indexOfNearest(xu);
    int i1 = _s1.indexOfNearest(x1);
    double[] mps = new double[1];
    if (!checkShift(iu,i1,mps))
      throw new UnphysicalShiftException(xu,x1,iu,i1,mps[0]);
    List<Integer> list = _errorMap.get(-1);
    if (list==null)
      list = new LinkedList<>();
    list.add(iu);
    list.add(i1);
    _errorMap.put(-1,list);
  }

  /**
   * 
   * @param xu
   * @param x1
   * @param x2
   * @throws UnphysicalShiftException 
   */
  public void add(double xu, double x1, double x2) throws 
      UnphysicalShiftException
  {
    int iu = _su.indexOfNearest(xu);
    int i1 = _s1.indexOfNearest(x1);
    double[] mps = new double[1];
    if (!checkShift(iu,i1,mps))
      throw new UnphysicalShiftException(xu,x1,iu,i1,mps[0]);
    int i2 = _s2.indexOfNearest(x2);
//    print("ErrorMapUtil.add: iu="+iu+", i1="+i1+", i2="+i2);
    List<Integer> list = _errorMap.get(i2);
    if (list==null)
      list = new LinkedList<>();
    list.add(iu);
    list.add(i1);
    _errorMap.put(i2,list);
  }

  /**
   * 
   * @param xu
   * @param x1
   * @param x2
   * @param x3
   * @throws UnphysicalShiftException 
   */
  public void add(double xu, double x1, double x2, double x3) throws 
  UnphysicalShiftException 
  {
    int iu = _su.indexOfNearest(xu);
    int i1 = _s1.indexOfNearest(x1);
    double[] mps = new double[1];
    if (!checkShift(iu,i1,mps))
      throw new UnphysicalShiftException(xu,x1,iu,i1,mps[0]);
    int i2 = _s2.indexOfNearest(x2);
    int i3 = _s3.indexOfNearest(x3);
    int n2 = _s2.getCount();
    int i23 = i2+(i3*n2);
    print("i2="+i2+", i3="+i3+", i23="+i23);
    List<Integer> list = _errorMap.get(i23);
    if (list==null)
      list = new LinkedList<>();
    list.add(iu);
    list.add(i1);
    _errorMap.put(i23,list);
  }

  /**
   * 
   * @return
   */
  public int[] getCoords() {
    List<Integer> list = _errorMap.get(-1);
    if (list==null) return null;
    int n = list.size();
    int[] l1 = new int[n];
    for (int i=0; i<n; i++)
      l1[i] = list.get(i);
    return l1;
  }

  /**
   * 
   * @return
   */
  public Map<Integer,int[]> getMap() {
    Map<Integer,int[]> map = new HashMap<>();
    for (Integer integer : _errorMap.keySet()) {
      List<Integer> list = _errorMap.get(integer);
      int n = list.size();
      int[] l1 = new int[n];
      for (int i=0; i<n; i++)
        l1[i] = list.get(i);
      map.put(integer,l1);
    }
    return map;
  }

  /**
   * 
   * @param g2
   * @param r2min
   * @param r2max
   * @return
   */
  public int[] getG2(int[] g2, double[] r2min, double[] r2max) {
    if (_errorMap.isEmpty()) return g2;

    List<Integer> g2NewList = new ArrayList<>();
    for (int i2 : _errorMap.keySet()) {
      g2NewList.add(i2);
    }

    List<Integer> g2KeepList = new ArrayList<>();
    int ng2 = g2.length;
    g2KeepList.add(g2[0]);
    g2KeepList.add(g2[ng2-1]);
    for (int ig2=1; ig2<ng2-1; ig2++) {
      int i2 = g2[ig2];
      double rmin = r2min[i2];
      double rmax = r2max[i2];
      boolean add = true;
      for (int i2New : g2NewList) {
        double dg = abs(i2-i2New);
        if (dg*(abs(rmax-rmin))<=1.0) {
          add = false;
          break;
        }
      }
      if (add)
        g2KeepList.add(i2);
    }

    g2NewList.addAll(g2KeepList);
    Collections.sort(g2NewList);
    int ng2New = g2NewList.size();
    int[] g2New = new int[ng2New];
    int ig2=0;
    for (int i2 : g2NewList) {
      g2New[ig2] = i2;
      ig2++;
    }
    return g2New;
  }

  /**
   * 
   * @param g2
   * @param g3
   * @param r2min
   * @param r2max
   * @param r3min
   * @param r3max
   * @return
   */
  public int[][] getG23(
      int[] g2, int[] g3, double[] r2min, double[] r2max,
      double[] r3min, double[] r3max) 
  {
    if (_errorMap.isEmpty()) return new int[][]{g2,g3};

    int n2 = _s2.getCount();
    List<Integer> g2NewList = new ArrayList<>();
    List<Integer> g3NewList = new ArrayList<>();
    for (int i23 : _errorMap.keySet()) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      if (!g2NewList.contains(i2)) g2NewList.add(i2);
      if (!g3NewList.contains(i3)) g3NewList.add(i3);
    }

    List<Integer> g2KeepList = new ArrayList<>();
    int ng2 = g2.length;
    g2KeepList.add(g2[0]);
    g2KeepList.add(g2[ng2-1]);
    for (int ig2=1; ig2<ng2-1; ig2++) {
      int i2 = g2[ig2];
      double rmin = r2min[i2];
      double rmax = r2max[i2];
      boolean add = true;
      for (int i2New : g2NewList) {
        double dg = abs(i2-i2New);
        if (dg*(abs(rmax-rmin))<=1.0) {
          add = false;
          break;
        }
      }
      if (add)
        g2KeepList.add(i2);
    }
    g2NewList.addAll(g2KeepList);
    Collections.sort(g2NewList);
    int ng2New = g2NewList.size();
    int[] g2New = new int[ng2New];
    int ig2=0;
    for (int i2 : g2NewList) {
      g2New[ig2] = i2;
      ig2++;
    }

    List<Integer> g3KeepList = new ArrayList<>();
    int ng3 = g3.length;
    g3KeepList.add(g3[0]);
    g3KeepList.add(g3[ng3-1]);
    for (int ig3=1; ig3<ng3-1; ig3++) {
      int i3 = g3[ig3];
      double rmin = r3min[i3];
      double rmax = r3max[i3];
      boolean add = true;
      for (int i3New : g3NewList) {
        double dg = abs(i3-i3New);
        if (dg*(abs(rmax-rmin))<=1.0) {
          add = false;
          break;
        }
      }
      if (add)
        g3KeepList.add(i3);
    }
    g3NewList.addAll(g3KeepList);
    Collections.sort(g3NewList);
    int ng3New = g3NewList.size();
    int[] g3New = new int[ng3New];
    int ig3=0;
    for (int i3 : g3NewList) {
      g3New[ig3] = i3;
      ig3++;
    }

    return new int[][]{g2New,g3New};
  }

  // Too complicated, not robust enough, and hopefully not necessary.
  // This method uses picks as seed points to create a pseudo-horizon
  // across the 2D image.
//  public float[][] getG1(float[][] g1,  double[] r1min, double[] r1max) {
//    if (_errorMap.isEmpty()) return g1;
//
//    int n2 = g1.length;
//    int n1 = g1[0].length;
//    Check.argument(n2==_s2.getCount(),"n2==_s2.getCount()");
//
//    // Loop through keys, which are trace indices
//    Map<Float,List<CoordX1X2>> newRowSeeds = new HashMap<>(); 
//    for (Integer i2 : _errorMap.keySet()) {
//      List<Float> list = _errorMap.get(i2);
//      int ncoords = list.size();
//
//      // get the i1 coordinates and find which layer it is in
//      for (int i=1; i<ncoords; i+=2) {
//        float i1New = list.get(i);
//        for (int i1=0; i1<n1-1; i1++) {
//          // FIXME assuming only one new layer can be created in an existing
//          //       layer, but there could be multiple points in a layer.
//          if (i1New>g1[i2][i1] && i1New<g1[i2][i1+1]) {
//            float row = i1+0.5f;
//            List<CoordX1X2> layerSeeds = newRowSeeds.get(row);
//            if (layerSeeds==null)
//              layerSeeds = new ArrayList<>();
//            layerSeeds.add(new CoordX1X2(i1New,i2));
//            newRowSeeds.put(row,layerSeeds);
//            break;
//          }
//        }
//      }
//    }
//
//    Map<Float,float[]> g1RowMap = new TreeMap<>();
//    for (int i1=0; i1<n1; i1++) {
//      float[] row = new float[n2];
//      for (int i2=0; i2<n2; i2++) {
//        row[i2] = g1[i2][i1];
//      }
//      g1RowMap.put((float)i1,row);
//    }
//
//    // Interpolate new layers from seed points
//    float[] x2f = rampfloat(0.0f,1.0f,n2);
//    for (Float row : newRowSeeds.keySet()) {
//      List<CoordX1X2> layerSeeds = newRowSeeds.get(row);
//      Collections.sort(layerSeeds);
//      int nseeds = layerSeeds.size();
//      float[] x2 = new float[nseeds];
//      float[] x1 = new float[nseeds];
//      for (int iseed=0; iseed<nseeds; iseed++) {
//        CoordX1X2 c = layerSeeds.get(iseed);
//        x2[iseed] = c._x2;
//        x1[iseed] = c._x1;
//      }
//      CubicInterpolator ci = new CubicInterpolator(x2,x1);
//      float[] newRow = ci.interpolate(x2f);
//      int irowl = (int)floor(row);
//      int irowu = (int) ceil(row);
//      g1RowMap.put(row,newRow);
//      if (irowl!=0) {
//        if (deleteRow(g1RowMap.get((float)irowl),newRow,r1min,r1max))
//          g1RowMap.remove((float)irowl);
//      }
//      if (irowu!=n1-1) {
//        if (deleteRow(g1RowMap.get((float)irowu),newRow,r1min,r1max))
//          g1RowMap.remove((float)irowu);
//      }
//    }
//
//    float[][] g1New = new float[n2][g1RowMap.size()];
//    int i1 = 0;
//    for (Float irow : g1RowMap.keySet()) {
//      float[] row = g1RowMap.get(irow);
//      for (int i2=0; i2<n2; i2++)
//        g1New[i2][i1] = row[i2];
//      i1++;
//    }
//    for (int i2=0; i2<n2; i2++) {
//      Check.argument(isMonotonic(g1New[i2]),"i2 "+i2+" is monotonic");
//    }
//
//    return g1New;
//  }

  public static int[] toSamples(Sampling su, Sampling s1, float[] u1) {
    if (u1==null) return null;
    int n = u1.length;
    int[] l1 = new int[n];
    for (int i=0; i<n; i+=2) {
      l1[i  ] = su.indexOfNearest(u1[i  ]);
      l1[i+1] = s1.indexOfNearest(u1[i+1]);
    }
    return l1;
  }

  /**
   * 
   * @param g1
   * @param l1
   * @param r1min
   * @param r1max
   * @return
   */
  public static int[] getG1(
      int[] g1, int[] l1, double[] r1min, double[] r1max) 
  {
    if (l1==null) return g1;

    List<Integer> g1NewList = new ArrayList<>();
    for (int ix1=1; ix1<l1.length; ix1+=2) {
      g1NewList.add(l1[ix1]);
    }

    List<Integer> g1KeepList = new ArrayList<>();
    int ng1 = g1.length;
    g1KeepList.add(g1[0]);
    g1KeepList.add(g1[ng1-1]);
    for (int ig1=1; ig1<ng1-1; ig1++) {
      int i1 = g1[ig1];
      double rmin = r1min[i1];
      double rmax = r1max[i1];
      boolean add = true;
      for (int i1New : g1NewList) {
        double dg = abs(i1-i1New);
        if (dg*(abs(rmax-rmin))<=1.0) {
          add = false;
          break;
        }
      }
      if (add)
        g1KeepList.add(i1);
    }

    g1NewList.addAll(g1KeepList);
    Collections.sort(g1NewList);
    int ng1New = g1NewList.size();
    int[] g1New = new int[ng1New];
    int ig1=0;
    for (int i1 : g1NewList) {
      g1New[ig1] = i1;
      ig1++;
    }
    return g1New;
  }

  public static float[] getG1(
      Sampling s1, float[] g1, float[] u1, double[] r1min, double[] r1max)
  {
    if (u1==null) return g1;
    List<Float> g1KeepList = new ArrayList<>();
    List<Float> g1KnownList = new ArrayList<>();
    int ng1 = g1.length;
    float first = g1[0    ];
    float last  = g1[ng1-1];
    g1KeepList.add(first);
    g1KeepList.add(last);

    // Add known shift x1 coordinates unless they are out of bounds.
    int nx1 = u1.length/2;
    int c = 0;
    for (int ix1=1; ix1<u1.length; ix1+=2) {
      float x1 = u1[ix1];
      if (x1>first && x1<last) {
        g1KnownList.add(x1);
        c++;
      }
    }
    if (c!=nx1)
      print("WARNING: Only "+c+" known shifts out of "+nx1+" are in bounds.");

    // Determine which g1 coordinates to keep.
    for (int ig1=1; ig1<ng1-1; ig1++) {
      int i1 = s1.indexOfNearest(g1[ig1]);
      double rmin = r1min[i1];
      double rmax = r1max[i1];
      boolean add = true;
      for (float x1New : g1KnownList) {
        int i1New = s1.indexOfNearest(x1New);
        // int i1New = (int)(x1New+0.5f);
        double dg = abs(i1-i1New);
        if (dg*(abs(rmax-rmin))<1.0) {
          add = false;
          break;
        }
      }
      if (add)
        g1KeepList.add(g1[ig1]);
    }

    // Combine the known shift x1 coordinates and the remaining g1 coordinates.
    g1KnownList.addAll(g1KeepList);
    Collections.sort(g1KnownList);
    int ng1New = g1KnownList.size();
    float[] g1New = new float[ng1New];
    int ig1=0;
    for (float x1 : g1KnownList) {
      g1New[ig1] = x1;
      ig1++;
    }
    return g1New;
  }

  public static void main(String[] args) {
    int nl = 169;
    int n1 = 1501;
    int n2 = 1211;
    int n3 = 826;
    double dt = 2.0;
    double dx = 22.5;
    double dy = 12.5;
    Sampling su = new Sampling(nl,dt,-200.0);
    Sampling s1 = new Sampling(n1,dt,0.0);
    Sampling s2 = new Sampling(n2,dx,1001);
    Sampling s3 = new Sampling(n3,dy,566);

    // 2D
    KnownShiftUtil emu = new KnownShiftUtil(su,s1,s2);
    try {
      emu.add(-100.0,5.0,1200.0);
      emu.add(-75.0,12.0,1200.0);
      emu.add(-150.0,300.0,1200.0);
      emu.add(0.0,0.0,1800.0);
    } catch (UnphysicalShiftException e) {
      e.printStackTrace();
    }
    Map<Integer,int[]> map = emu.getMap();
    for (Integer integer : map.keySet()) {
      int[] l1 = map.get(integer);
      print("il coordinate array for index "+integer+":");
      dump(l1);
    }

    //3
    emu = new KnownShiftUtil(su,s1,s2,s3);
    try {
      emu.add(-100.0,5.0,1200.0);
      emu.add(-75.0,12.0,1200.0);
      emu.add(-150.0,300.0,1200.0);
      emu.add(0.0,0.0,1800.0);
    } catch (UnphysicalShiftException e) {
      e.printStackTrace();
    }
    map = emu.getMap();
    for (Integer integer : map.keySet()) {
      int[] l1 = map.get(integer);
      print("il coordinate array for index "+integer+":");
      dump(l1);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private
  private final Map<Integer,List<Integer>> _errorMap;
  private final Sampling _su, _s1, _s2, _s3;
  private final double[] _r1Min;

  private boolean checkShift(int iu, int i1, double[] mps) {
    mps[0] = 0.0;
    for (int i=0; i<i1; i++)
      mps[0] += _r1Min[i];
    if (iu<=mps[0]) return false;
    return true;
  }

  private static void print(String s) {
    System.out.println(s);
  }

//  private boolean deleteRow(
//      float[] row, float[] newRow, double[] r1min, double[] r1max) 
//  {
//    boolean delete = false;
//    int n2 = row.length;
//    for (int i2=0; i2<n2; i2++) {
//      float f1 = row[i2];
//      double rmin = r1min[(int)(f1+0.5f)];
//      double rmax = r1max[(int)(f1+0.5f)];
//      if (!checkDistance(f1,newRow[i2],rmin,rmax)) {
//        delete = true;
//        break;
//      }
//    }
//    return delete;
//  }

//  private static class CoordX1X2 implements Comparable<CoordX1X2> {
//
//    public CoordX1X2(float x1, float x2) {
//      _x1 = x1;
//      _x2 = x2;
//    }
//
//    @Override
//    public int compareTo(CoordX1X2 c) {
//      if (this._x2<c._x2) return -1;
//      if (this._x2>c._x2) return  1;
//      return 0;
//    }
//    
//    private final float _x1, _x2;
//  }

  public class UnphysicalShiftException extends Throwable {

    public UnphysicalShiftException(
        double xu, double x1, int iu, int i1, double mps) 
    {
      super("Sample: "+x1+", Sample index: "+i1+
            ", Shift: "+xu+", Shift index: "+iu+
            ", Max physical shift index: "+mps);
      _xu = xu;
      _x1 = x1;
      _iu = iu;
      _i1 = i1;
      _mps = mps;
    }

    public double getShift() { return _xu; }
    public int getShiftIndex() { return _iu; }
    public double getSample() { return _x1; }
    public int getSampleIndex() { return _i1; }
    public double getMaxPhysicalShift() { return _mps; }

    private static final long serialVersionUID = 1L;
    private final double _xu, _x1, _mps;
    private final int _iu, _i1;
  }

}
