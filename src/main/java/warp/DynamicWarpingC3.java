package warp;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.dsp.SincInterp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping for 3 converted wave data images. Time shifts for PP, PS1,
 * and PS2 images are found by simultaneously warping all three images. Time
 * shifts can be found for 3 image pairs independently with the relationships:
 *     tps1(tpp) = tpp + u1(tpp),
 *     tps2(tpp) = tpp + u2(tpp),
 *     tps2(tps1) = tps1 + uS(tps1),
 * but these shifts should be related such that u2 = u1 + uS. Using this
 * relationship and all 3 images simultaneously should better constrain the time
 * shifts.
 */
public class DynamicWarpingC3 {

  /**
   * Constructs a 1D dynamic warping.
   * @param u1Min minimum shift between PP and PS1.
   * @param u1Max maximum shift between PP and PS1.
   * @param uSMin minimum shift between PS1 and PS2.
   * @param uSMax maximum shift between PS1 and PS2.
   * @param s1 first dimension sampling for warping.
   */
  public DynamicWarpingC3(
      double u1Min, double u1Max,
      double uSMin, double uSMax,
      Sampling s1)
  {
    this(u1Min,u1Max,uSMin,uSMax,s1,null,null);
  }

  /**
   * Constructs a 2D dynamic warping.
   * @param u1Min minimum shift between PP and PS1.
   * @param u1Max maximum shift between PP and PS1.
   * @param uSMin minimum shift between PS1 and PS2.
   * @param uSMax maximum shift between PS1 and PS2.
   * @param s1 first dimension sampling for warping.
   * @param s2 second dimension sampling.
   */
  public DynamicWarpingC3(
      double u1Min, double u1Max,
      double uSMin, double uSMax,
      Sampling s1, Sampling s2)
  {
    this(u1Min,u1Max,uSMin,uSMax,s1,s2,null);
  }

  /**
   * Constructs a 2D dynamic warping.
   * @param u1Min minimum shift between PP and PS1.
   * @param u1Max maximum shift between PP and PS1.
   * @param uSMin minimum shift between PS1 and PS2.
   * @param uSMax maximum shift between PS1  and PS2.
   * @param s1 first dimension sampling for warping.
   * @param s2 second dimension sampling.
   * @param s3 third dimension sampling.
   */
  public DynamicWarpingC3(
      double u1Min, double u1Max,
      double uSMin, double uSMax,
      Sampling s1, Sampling s2, Sampling s3)
  {
    double ds = s1.getDelta(); // shift sampling interval
    int iu1Min = (int) ceil(u1Min/ds);
    int iu1Max = (int)floor(u1Max/ds);
    int iuSMin = (int) ceil(uSMin/ds);
    int iuSMax = (int)floor(uSMax/ds);
    _su1 = new Sampling(1+iu1Max-iu1Min,ds,iu1Min*ds);
    _suS = new Sampling(1+iuSMax-iuSMin,ds,iuSMin*ds);
    print("ns1="+_su1.getCount()+", nsS="+_suS.getCount());
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _r1Min1 =  0.0;
    _r1Max1 =  1.0;
    _r1MinS =  0.0;
    _r1MaxS =  1.0;
    _r2Min  = -1.0;
    _r2Max  =  1.0;
    _r3Min  = -1.0;
    _r3Max  =  1.0;
    _si = new SincInterp();
    // _si.setExtrapolation(Extrapolation.CONSTANT);
    _si.setExtrapolation(Extrapolation.ZERO);
    _ui = new ShiftInterp();
    _ui.setMethod(ShiftInterp.Method.LINEAR);
  }

  public Sampling getShiftSampling1() {
    return _su1;
  }

  public Sampling getShiftSamplingS() {
    return _suS;
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min1 lower bound on strain 1 in 1st dimension.
   * @param r1Max1 upper bound on strain 1 in 1st dimension.
   * @param r1MinS lower bound on strain S in 1st dimension.
   * @param r1MaxS upper bound on strain S in 1st dimension.
   */
  public void setStrainLimits(
      double r1Min1, double r1Max1,
      double r1MinS, double r1MaxS)
  {
    setStrainLimits(r1Min1,r1Max1,r1MinS,r1MaxS,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min1 lower bound on strain in 1st dimension.
   * @param r1Max1 upper bound on strain in 1st dimension.
   * @param r1MinS lower bound on strain S in 1st dimension.
   * @param r1MaxS upper bound on strain S in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
    double r1Min1, double r1Max1,
    double r1MinS, double r1MaxS,
    double r2Min, double r2Max)
  {
    setStrainLimits(r1Min1,r1Max1,r1MinS,r1MaxS,r2Min,r2Max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min1 lower bound on strain in 1st dimension.
   * @param r1Max1 upper bound on strain in 1st dimension.
   * @param r1MinS lower bound on strain S in 1st dimension.
   * @param r1MaxS upper bound on strain S in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   * @param r3Min lower bound on strain in 3rd dimension.
   * @param r3Max upper bound on strain in 3rd dimension.
   */
  public void setStrainLimits(
    double r1Min1, double r1Max1,
    double r1MinS, double r1MaxS,
    double r2Min, double r2Max,
    double r3Min, double r3Max)
  {
    _r1Min1 = r1Min1; _r1Max1 = r1Max1;
    _r1MinS = r1MinS; _r1MaxS = r1MaxS;
    _r2Min = r2Min; _r2Max = r2Max;
    _r3Min = r3Min; _r3Max = r3Max;
  }

  public void setInterpolationMethod(ShiftInterp.Method method) {
    _ui.setMethod(method);
  }

  public float[][] findShifts(
      Sampling sf, float[] f,
      Sampling sg, float[] g,
      Sampling sh, float[] h,
      float[] g1)
  {
    float dfi = 1.0f/(float)sf.getDelta();
    int ng1 = g1.length;
    int[] g1i = new int[ng1];
    for (int ig1=0; ig1<ng1; ig1++)
      g1i[ig1] = (int)(g1[ig1]*dfi+0.5f);
    printInfo(ng1,1,1);

    double[] r1Min1 = getR1Min1();
    double[] r1Max1 = getR1Max1();
    double[] r1MinS = getR1MinS();
    double[] r1MaxS = getR1MaxS();

    Stopwatch sw = new Stopwatch();
    sw.start();
    float[][][] e = computeErrors(sf,f,sg,g,sh,h);
    print("Computed errors in "+sw.time()+" seconds");
    sw.restart();
    float[][][] d = makeD(ng1);
    int[][][][] m = makeM(ng1);
    double de  =  _s1.getDelta();
    double du1 = _su1.getDelta();
    double duS = _suS.getDelta();
    accumulateSparse(1,r1Min1,r1Max1,r1MinS,r1MaxS,g1i,de,du1,duS,e,d,m);
    print("Accumulated errors in "+sw.time()+" seconds");
    sw.restart();
    float[] u1 = new float[ng1];
    float[] uS = new float[ng1];
    backtrack(-1,d[ng1-1],m,_su1,u1,_suS,uS);
    print("Backtracked accumulated errors in "+sw.time()+" seconds");

    // // Check slopes
    // float lastLag = (float)_sl1.getLast();
    // for (int ig1=1; ig1<ng1; ig1++) {
    //   int i1 = g1i[ig1];
    //   int i1m1 = g1i[ig1-1];
    //   float n = u[ig1] - u[ig1-1];
    //   float d = i1-i1m1;
    //   float r = n/d;
    //   assert (r>=r1Min1[i1] && r<=r1Max1[i1]) || u[ig1]==lastLag :
    //     "n="+n+", d="+d+", r="+r;
    // }

    // print("u1:"); dump(u1);
    // print("uS:"); dump(uS);
    float[] u1i = _ui.interpolate(sf,g1,u1);
    float[] uSi = _ui.interpolate(sg,g1,uS);
    return new float[][]{u1i,uSi};
  }

  public float[][][] computeErrors(
      Sampling sf, float[] f,
      Sampling sg, float[] g,
      Sampling sh, float[] h)
  {
    Sampling s1 = _su1;
    Sampling sS = _suS;
    Sampling se = _s1;
    int n1 = s1.getCount();
    int nS = sS.getCount();
    int ne = se.getCount();
    int nf = sf.getCount();
    int ng = sg.getCount();
    int nh = sh.getCount();
    double de = se.getDelta();
    double dh = sh.getDelta();
    double dg = sg.getDelta();
    double fe = se.getFirst();
    double fh = sh.getFirst();
    double fg = sg.getFirst();
    float[][][] e = new float[ne][n1][nS];
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    float[] hi = new float[ne];
    _si.interpolate(sf,f,se,fi);
    for (int iS=0; iS<nS; iS++) {
      double vS = sS.getValue(iS);
      for (int i1=0; i1<n1; i1++) {
        double v1 = s1.getValue(i1);
        _si.interpolate(nh,dh,fh,h,ne,de,fe+vS+v1,hi);
        _si.interpolate(ng,dg,fg,g,ne,de,fe+v1,   gi);
        for (int ie=0; ie<ne; ie++) {
          float e1 = error(fi[ie],gi[ie]);
          float e2 = error(fi[ie],hi[ie]);
          float e3 = error(gi[ie],hi[ie]);
          e[ie][i1][iS] = e1+e2+e3;
        }
      }
    }
    return e;
  }

  /*package*/ static void linearInterp(
      int ie, int je, int is, float r,
      int il1s, int il1e, int ilSs, int ilSe,
      float[][] dm, float[][][] e, boolean interp1)
  {
    for (int x=je+is; x!=ie; x+=is) {
      float ky = r*(ie-x);
      int k1 = (int)ky;
      if (ky<0.0f) --k1;
      int k2 = k1+1;
      float w1 = k2-ky;
      float w2 = 1.0f-w1;
      if (interp1) {
        for (int il1=il1s; il1<il1e; il1++) {
          int l1 = il1+k1;
          int l2 = il1+k2;
          for (int ilS=ilSs; ilS<ilSe; ilS++)
            dm[il1][ilS] += w1*e[x][l1][ilS]+w2*e[x][l2][ilS];
        }
      } else {
        for (int il1=il1s; il1<il1e; il1++)
          for (int ilS=ilSs; ilS<ilSe; ilS++)
            dm[il1][ilS] += w1*e[x][il1][ilS+k1]+w2*e[x][il1][ilS+k2];
      }
    }
  }

  /*package*/ static void bilinearInterp(
      int ie, int je, int is, float r1, float rS,
      int il1s, int il1e, int ilSs, int ilSe,
      float[][] dm, float[][][] e)
  {
    for (int x=je+is; x!=ie; x+=is) {
      float k1y = r1*(ie-x);
      float kSy = rS*(ie-x);
      int k1i1 = (int)k1y;
      int kSi1 = (int)kSy;
      if (k1y<0.0f) --k1i1;
      if (kSy<0.0f) --kSi1;
      int k1i2 = k1i1+1;
      int kSi2 = kSi1+1;
      float d1 = k1i2-k1y;
      float dS = kSi2-kSy;
      float d1S = d1*dS;
      for (int il1=il1s; il1<il1e; il1++) {
        int l1 = il1+k1i1;
        int l2 = il1+k1i2;
        for (int ilS=ilSs; ilS<ilSe; ilS++) {
          float e11 = e[x][l1][kSi1+ilS];
          float e21 = e[x][l1][kSi2+ilS];
          float e12 = e[x][l2][kSi1+ilS];
          float e22 = e[x][l2][kSi2+ilS];
          dm[il1][ilS] += 
            e11 + (e21-e11)*dS + (e12-e11)*d1 + (e11-e21-e12+e22)*d1S;
        }
      }
    }
  }

  /*package*/ static void accumulateSparse(
      int dir,
      double[] r1Min, double[] r1Max,
      double[] rSMin, double[] rSMax,
      int[] g, double de, double du1, double duS,
      float[][][] e, float[][][] d, int[][][][] m)
  {
    int nl1  = e[0].length;
    int nlS  = e[0][0].length;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ibg = (dir>0)?0:ngm1; // beginning index
    int ieg = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int isp = ibg; // sparse grid index
    int ie = g[isp]; // error index

    // Initialize accumulation values
    for (int il1=0; il1<nl1; ++il1)
      for (int ilS=0; ilS<nlS; ++ilS)
        d[isp][il1][ilS] = e[ie][il1][ilS];
    isp += is;
    // Loop over all sparse grid points.
    for (; isp!=ieg; isp+=is) {
      int ispm1 = isp-is; // previous sparse grid index
      ie = g[isp]; // new error index
      int je = g[ispm1]; // min error index for interpolation.
      int dg = ie-je; // sparse grid delta, possibly negative
      int k1Min, k1Max, kSMin, kSMax;
      if (dg>0) { // forward
        k1Min = (int) ceil(-r1Max[ie]*dg*de/du1);
        k1Max = (int)floor(-r1Min[ie]*dg*de/du1);
        kSMin = (int) ceil(-rSMax[ie]*dg*de/duS);
        kSMax = (int)floor(-rSMin[ie]*dg*de/duS);
      } else { // reverse
        k1Min = (int) ceil(-r1Min[ie]*dg*de/du1);
        k1Max = (int)floor(-r1Max[ie]*dg*de/du1);
        kSMin = (int) ceil(-rSMin[ie]*dg*de/duS);
        kSMax = (int)floor(-rSMax[ie]*dg*de/duS);
      }
      k1Min = k1Min>k1Max ? k1Max : k1Min;
      kSMin = kSMin>kSMax ? kSMax : kSMin;
      float[][] dm = new float[nl1][nlS]; // temp accumulate array.
      fill(Float.MAX_VALUE,d[isp]);
      // loop over all slope indices
      for (int k1=k1Min; k1<=k1Max; k1++) {
        int il1s = max(0,-k1);      // lag 1 start
        int il1e = min(nl1,nl1-k1); // lag 1 end
        float r1 = (float)k1/(float)dg; // slope
        for (int kS=kSMin; kS<=kSMax; kS++) {
          int ilSs = max(0,-kS);      // lag S start
          int ilSe = min(nlS,nlS-kS); // lag S end
          float rS = (float)kS/(float)dg; // slope
          // print("dg="+dg+", k1="+k1+", kS="+kS);

          // Start the sum with the previous accumulated error (d) and alignment
          // error (e). d is accumulated on a sparse grid and here we add e at
          // the next sparse grid point. In between we will need to interpolate
          // and keep adding to dm. Our optimum path will follow the minimum dm.
          for (int il1=il1s; il1<il1e; il1++) {
            int il1k1 = il1+k1;
            for (int ilS=ilSs; ilS<ilSe; ilS++)
              dm[il1][ilS] = d[ispm1][il1k1][ilS+kS] + e[ie][il1][ilS];
          }
          if (rS==0 && r1==0) { // zero slopes, no interpolation necessary
            for (int x=je+is; x!=ie; x+=is)
              for (int il1=il1s; il1<il1e; il1++)
                for (int ilS=ilSs; ilS<ilSe; ilS++)
                  dm[il1][ilS] += e[x][il1][ilS];
          } else if (rS==0 && r1!=0) {
            linearInterp(ie,je,is,r1,il1s,il1e,ilSs,ilSe,dm,e,true);
          } else if (rS!=0 && r1==0) {
            linearInterp(ie,je,is,rS,il1s,il1e,ilSs,ilSe,dm,e,false);
          } else {
            bilinearInterp(ie,je,is,r1,rS,il1s,il1e,ilSs,ilSe,dm,e);
          }

          // update previous errors and record moves.
          for (int il1=il1s; il1<il1e; il1++) {
            for (int ilS=ilSs; ilS<ilSe; ilS++) {
              if (dm[il1][ilS]<d[isp][il1][ilS]) {
                d[isp][il1][ilS] = dm[il1][ilS];
                if (m!=null) {
                  m[ispm1][il1][ilS][0] = k1;
                  m[ispm1][il1][ilS][1] = kS;
                }
              }
            }
          }

        } // end kS loop
      } // end k1 loop
    } // end sparse grid loop
  }

  /*package*/ static void backtrack(
      int dir, float[][] d, int[][][][] m,
      Sampling su1, float[] u1, Sampling suS, float[] uS)
  {
    int n1  = m.length;
    int nl1 = m[0].length;
    int nlS = m[0][0].length;
    int n1m  = n1-1;
    int nl1m = nl1-1;
    int nlSm = nlS-1;
    int ib = (dir>0)?0:n1m; // begining index
    int ie = (dir>0)?n1m:0; // end index
    int is = (dir>0)?1:-1;  // stride
    int i1 = ib;
    // Set initial lag indices and minimum value.
    int il1Min = (dir>0)?0:nl1m;
    int ilSMin = (dir>0)?0:nlSm;
    float dm = d[il1Min][ilSMin];
    // Search for minimum lag indices.
    for (int il1=0; il1<nl1; il1++) {
      for (int ilS=0; ilS<nlS; ilS++) {
        if (d[il1][ilS]<dm) {
          dm = d[il1][ilS];
          il1Min = il1;
          ilSMin = ilS;
        }
      }
    }
    u1[i1] = (float)su1.getValue(il1Min);
    uS[i1] = (float)suS.getValue(ilSMin);
    while (i1!=ie) {
      i1 += is;
      int il1 = m[i1][il1Min][ilSMin][0];
      int ilS = m[i1][il1Min][ilSMin][1];
      il1Min += il1;
      ilSMin += ilS;
      u1[i1] = (float)su1.getValue(il1Min);
      uS[i1] = (float)suS.getValue(ilSMin);
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  private final Sampling _su1, _suS;
  private Sampling _s1, _s2, _s3;
  private double _r1Min1, _r1Max1, _r1MinS, _r1MaxS;
  private double _r2Min, _r2Max;
  private double _r3Min, _r3Max;
  private double[] _r1Min1A, _r1Max1A, _r1MinSA, _r1MaxSA;
  private SincInterp _si;
  private ShiftInterp _ui;
  private final static float EPOW = 2.0f;
  private static final float BYTES_TO_MB = 1.0f/1000000.0f;

  private void printInfo(int ng1, int ng2, int ng3) {
    int n1 = _s1.getCount();
    int n2 = _s2==null ? 1 : _s2.getCount();
    int n3 = _s3==null ? 1 : _s3.getCount();
    int nl1 = _su1.getCount();
    int nlS = _suS.getCount();
    float e1Mem = (float)ng1*n2*n3*nl1*nlS*4.0f*BYTES_TO_MB;
    float e2Mem = (float)ng1*ng2*n3*nl1*nlS*4.0f*BYTES_TO_MB;
    float e3Mem = (float)ng1*ng2*ng3*nl1*nlS*4.0f*BYTES_TO_MB;
    print("DynamicWarpingC info:");
    print("  Input data samples (n1,n2,n3): ("+n1+","+n2+","+n3+")");
    print("  Coarse grid samples: (ng1,ng2,ng3): ("+ng1+","+ng2+","+ng3+")");
    print("  Number of lags: nl1="+nl1+", nlS="+nlS);
    print("  Alignment error smooth 1 memory: "+e1Mem+" MB");
    print("  Alignment error smooth 2 memory: "+((n2>1)?(e2Mem+" MB"):"NA"));
    print("  Alignment error smooth 3 memory: "+((n3>1)?(e3Mem+" MB"):"NA"));
  }

  private double[] getR1Min1() {
    double[] r1Min1;
    if (_r1Min1A==null)
      r1Min1 = filldouble(_r1Min1,_s1.getCount());
    else
      r1Min1 = _r1Min1A;
    return r1Min1;
  }

  private double[] getR1MinS() {
    double[] r1MinS;
    if (_r1MinSA==null)
      r1MinS = filldouble(_r1MinS,_s1.getCount());
    else
      r1MinS = _r1MinSA;
    return r1MinS;
  }

  private double[] getR1Max1() {
    double[] r1Max1;
    if (_r1Max1A==null)
      r1Max1 = filldouble(_r1Max1,_s1.getCount());
    else
      r1Max1 = _r1Max1A;
    return r1Max1;
  }

  private double[] getR1MaxS() {
    double[] r1MaxS;
    if (_r1MaxSA==null)
      r1MaxS = filldouble(_r1MaxS,_s1.getCount());
    else
      r1MaxS = _r1MaxSA;
    return r1MaxS;
  }

  private float[][][] makeD(int ng) {
    return new float[ng][_su1.getCount()][_suS.getCount()];
  }

  private int[][][][] makeM(int ng) {
    return new int[ng][_su1.getCount()][_suS.getCount()][2];
  }

  private float error(float f, float g) {
    return pow(abs(f-g),EPOW);
  }

  /**
   * Returns smooth alignment errors on the sparse grid defined
   * by the indices of g.
   * @param e 3D array of alignment errors.
   * @param rmin minimum slope.
   * @param rmax maximum slope.
   * @param g sparse grid indices.
   * @return smoothed alignment errors with size
   *  [g.length][e[0].length].
   */
  private static float[][][] smoothErrors(
      double[] rMin1, double[] rMax1,
      double[] rMinS, double[] rMaxS,
      int[] g, double du1, double duS, double de, float[][][] e)
  {
    int ng = g.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;
    float[][][] ef = new float[ng][nl1][nlS];
    float[][][] er = new float[ng][nl1][nlS];
    float[][][] es = new float[ng][nl1][nlS];
    accumulateSparse( 1,rMin1,rMax1,rMinS,rMaxS,g,de,du1,duS,e,ef,null);
    accumulateSparse(-1,rMin1,rMax1,rMinS,rMaxS,g,de,du1,duS,e,er,null);
    float s = 1.0f/e.length; // scale
    for (int i1=0; i1<ng; i1++) {
      for (int il1=0; il1<nl1; il1++) {
        for (int ilS=0; ilS<nlS; ilS++) {
          float v = s*(ef[i1][il1][ilS]+er[i1][il1][ilS]-e[g[i1]][il1][ilS]);
          es[i1][il1][ilS] = Float.isInfinite(v) ? Float.MAX_VALUE : v;
        }
      }
    }
    return es;
  }

  private static void print(String s) {
    System.out.println(s);
  }

}