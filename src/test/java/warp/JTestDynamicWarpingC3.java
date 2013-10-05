package warp;

import junit.framework.TestCase;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import utils.Synthetic;
import static warp.DynamicWarpingC3.*;

public class JTestDynamicWarpingC3 extends TestCase {

  public void testLinearInterp() {
    print("Testing linear interpolation...");

    // Test lag 1 interpolation
    for (int n1=2; n1<100; n1++) {
      int nl1 = n1; // make lag1 and i1 a square
      int nlS = 1;
      int is;
      float[][] dm = new float[nl1][nlS]; // accumulation array
      float[][][] e = new float[n1][nl1][nlS]; // error array
      for (int i1=0; i1<n1; i1++) {
        for (int il1=0; il1<nl1; il1++) {
          if (i1==il1)
            e[i1][il1][0] = i1; // populate diagonal
          else
            e[i1][il1][0] = 0.0f;
        }
      }
      // Sum along diagonal excluding first and last sample. This is the 
      // answer we should get after interpolation.
      int sum = 0;
      for (int i1=1; i1<n1-1; i1++) sum+=i1;

      // Alternate accumulation direction
      for (int i=0; i<2; i++) {
        if (i==0)
          is = 1;
        else
          is = -1;
        int ie = is>0 ? n1-1 :    0;
        int je = is>0 ?    0 : n1-1;
        float r = -1.0f;
        int il1s = is>0 ? nl1-1 : 0;
        int il1e = is>0 ? nl1   : 1;
        int ilSs = nlS-1;
        int ilSe = nlS;
        linearInterp(ie,je,is,r,il1s,il1e,ilSs,ilSe,dm,e,true);
        assert sum==dm[il1s][0];
      }
    }

    // Test lag S interpolation
    for (int n1=2; n1<100; n1++) {
      int nlS = n1; // make lagS and i1 a square
      int nl1 = 1;
      int is;
      float[][] dm = new float[nl1][nlS]; // accumulation array
      float[][][] e = new float[n1][nl1][nlS]; // error array
      for (int i1=0; i1<n1; i1++) {
        for (int ilS=0; ilS<nlS; ilS++) {
          if (i1==ilS)
            e[i1][0][ilS] = i1; // populate diagonal
          else
            e[i1][0][ilS] = 0.0f;
        }
      }
      // Sum along diagonal excluding first and last sample. This is the 
      // answer we should get after interpolation.
      int sum = 0;
      for (int i1=1; i1<n1-1; i1++) sum+=i1;

      // Alternate accumulation direction
      for (int i=0; i<2; i++) {
        if (i==0)
          is = 1;
        else
          is = -1;
        int ie = is>0 ? n1-1 :    0;
        int je = is>0 ?    0 : n1-1;
        float r = -1.0f;
        int ilSs = is>0 ? nlS-1 : 0;
        int ilSe = is>0 ? nlS   : 1;
        int il1s = nl1-1;
        int il1e = nl1;
        linearInterp(ie,je,is,r,il1s,il1e,ilSs,ilSe,dm,e,false);
        assert sum==dm[0][ilSs];
      }
    }
  }

  public void testComputeErrors() {
    // if (true) return;
    print("Testing computeErrors...");
    // Create samplings and data
    Sampling sf = new Sampling(101,0.002,0.0);
    Sampling sg = new Sampling(101,0.002,0.0);
    Sampling sh = new Sampling(111,0.002,0.0);
    // TODO understand/implement test when samplings are different.
    // Sampling sf = new Sampling(201,0.002,0.0);
    // Sampling sg = new Sampling(101,0.004,0.0);
    // Sampling sh = new Sampling( 82,0.005,0.0);
    float[] f = fillfloat(100.0f,sf.getCount());
    float[] g = fillfloat(200.0f,sg.getCount());
    float[] h = fillfloat(300.0f,sh.getCount());
    // Sampling for warping
    int n1w = 61;
    Sampling s1 = new Sampling(n1w,sf.getDelta(),sf.getFirst());

    // Make min/max shifts
    double s1Min = 0.0;
    double s1Max = sf.getLast()-s1.getLast();
    double sSMin = 0.0;  // 0 ms
    double sSMax = sh.getLast()-sg.getLast();
    print("s1Max="+s1Max+", sSMax="+sSMax);
    DynamicWarpingC3 dw = new DynamicWarpingC3(s1Min,s1Max,sSMin,sSMax,s1);
    Sampling su1 = dw.getSamplingU1();
    Sampling suS = dw.getSamplingUS();
    Almost a = new Almost();

    double ff = sf.getFirst();
    double fg = sg.getFirst();
    double fh = sh.getFirst();
    double df = sf.getDelta();
    double dg = sg.getDelta();
    double dh = sh.getDelta();
    double d1 = s1.getDelta();
    double du1 = su1.getDelta();
    double duS = suS.getDelta();
    for (double v=s1.getFirst(); v<=s1.getLast(); v+=2*d1) {
      for (double u1=s1Min; u1<=s1Max; u1+=2*du1) {
        for (double uS=sSMin; uS<=sSMax; uS+=2*duS) {
          // print("v="+v);
          double vf = v;
          double vg = v+u1;
          double vh = v+u1+uS;
          int i1f = sf.indexOfNearest(vf);
          int i1g = sg.indexOfNearest(vg);
          int i1h = sh.indexOfNearest(vh);
          // print("vf="+vf+", vg="+vg+", vh="+vh+
          //   ", i1f="+i1f+", i1g="+i1g+", i1h="+i1h);
          // double vf = sf.getValue(i1f);
          // double vg = sg.getValue(i1g);
          // double vh = sh.getValue(i1h);
          float[] fc = copy(f);
          float[] gc = copy(g);
          float[] hc = copy(h);
          fc[i1f] = gc[i1g] = hc[i1h] = 1000.0f;
          float[][][] e = dw.computeErrors(sf,fc,sg,gc,sh,hc);

          // Get the indices of the minimum error
          int i1Min, iu1Min, iuSMin;
          i1Min = iu1Min = iuSMin = 0;
          float min = Float.MAX_VALUE;
          for (int i1=0; i1<s1.getCount(); i1++) {
            for (int iu1=0; iu1<su1.getCount(); iu1++) {
              for (int iuS=0; iuS<suS.getCount(); iuS++) {
                if (e[i1][iu1][iuS]<min) {
                  i1Min = i1;
                  iu1Min = iu1;
                  iuSMin = iuS;
                  min = e[i1][iu1][iuS];
                }
              }
            }
          }
          // print("i1Min="+i1Min+",iu1Min="+iu1Min+", iuSMin="+iuSMin);
          double  v1Min =  s1.getValue( i1Min);
          double vu1Min = su1.getValue(iu1Min);
          double vuSMin = suS.getValue(iuSMin);
          // double  v1Min = sf.getValue( i1Min);
          // double vu1Min = sg.getValue(iu1Min);
          // double vuSMin = sh.getValue(iuSMin);
          double vfc = v1Min;
          double vgc = v1Min+vu1Min;
          double vhc = v1Min+vu1Min+vuSMin;
          String info = "v="+v+", u1="+u1+", uS="+uS;
          assert a.equal(vfc,vf):"vfc==vf, vfc="+vfc+", vf="+vf+", "+info;
          assert a.equal(vgc,vg):"vgc==vg, vgc="+vgc+", vg="+vg+", "+info;
          assert a.equal(vhc,vh):"vhc==vh, vhc="+vhc+", vh="+vh+", "+info;
        }
      }
    }
  }

  public void testAccumulateSparse() {
    // if (true) return;
    print("Testing AccumulateSparse...");
    Sampling se  = new Sampling(101,0.002,0.0);
    Sampling su1 = new Sampling(53,0.002,0.0);
    Sampling suS = new Sampling(22,0.002,0.0);
    int ne  =  se.getCount();
    int nu1 = su1.getCount();
    int nuS = suS.getCount();
    double de  =  se.getDelta();
    double du1 = su1.getDelta();
    double duS = suS.getDelta();

    // Test if errors are all 1.0 and strain is 0
    double[] r1Min = filldouble(0.0,ne);
    double[] r1Max = filldouble(0.0,ne);
    double[] rSMin = filldouble(0.0,ne);
    double[] rSMax = filldouble(0.0,ne);
    float[][][] e = fillfloat(1.0f,nuS,nu1,ne);
    for (int dg=1; dg<ne; dg++) {
      int[] g = Subsample.subsample(ne,dg);
      int ng = g.length;
      int ngm = ng-1;
      float[][][] d = new float[ng][nu1][nuS];
      accumulateSparse(1,r1Min,r1Max,rSMin,rSMax,g,de,du1,duS,e,d,null);
      accumulateAssert(nu1,nuS,ne,ngm,d);
      accumulateSparseP(1,r1Min,r1Max,rSMin,rSMax,g,de,du1,duS,e,d,null);
      accumulateAssert(nu1,nuS,ne,ngm,d);
    }
  }

  // public void testAccumulateSparseTime() {
  //   print("Testing AccumulateSparseTime...");
  //   Sampling se  = new Sampling(101,0.002,0.0);
  //   Sampling su1 = new Sampling(53,0.002,0.0);
  //   Sampling suS = new Sampling(22,0.002,0.0);
  //   int ne  =  se.getCount();
  //   int nu1 = su1.getCount();
  //   int nuS = suS.getCount();
  //   double de  =  se.getDelta();
  //   double du1 = su1.getDelta();
  //   double duS = suS.getDelta();
  //   double[] r1Min = filldouble(-1.0,ne);
  //   double[] r1Max = filldouble( 1.0,ne);
  //   double[] rSMin = filldouble(-1.0,ne);
  //   double[] rSMax = filldouble( 1.0,ne);
  //   float[][][] e = fillfloat(1.0f,nuS,nu1,ne);
  //   Stopwatch sw  = new Stopwatch();
  //   Stopwatch swp = new Stopwatch();
  //   for (int dg=1; dg<ne; dg++) {
  //     int[] g = Subsample.subsample(ne,dg);
  //     int ng = g.length;
  //     int ngm = ng-1;
  //     float[][][] d = new float[ng][nu1][nuS];
  //     sw.start();
  //     accumulateSparse(1,r1Min,r1Max,rSMin,rSMax,g,de,du1,duS,e,d,null);
  //     sw.stop();
  //     // accumulateAssert(nu1,nuS,ne,ngm,d);
  //     swp.start();
  //     accumulateSparseP(1,r1Min,r1Max,rSMin,rSMax,g,de,du1,duS,e,d,null);
  //     swp.stop();
  //     // accumulateAssert(nu1,nuS,ne,ngm,d);
  //   }
  //   print("accumulate serial: "+sw.time()+" seconds");
  //   print("accumulate parallel: "+swp.time()+" seconds");
  // }

  public void testFindShifts1D() {
    print("Testing findShifts1D...");
    float u1c = 0.024f;
    float uSc = 0.008f;
    float fpeak = 0.06f; // cycles/sample

    // Create samplings and data
    Sampling sf = new Sampling(101,0.002,0.0);
    Sampling sg = new Sampling(101,0.002,0.0);
    Sampling sh = new Sampling(111,0.002,0.0);
    float[] fe = Synthetic.makeRandomEvents(sf.getCount(),1);
    SincInterp si = new SincInterp();
    float[] ge = new float[sg.getCount()];
    float[] he = new float[sh.getCount()];
    si.interpolate(
      sf.getCount(),sf.getDelta(),sf.getFirst(),fe,
      sg.getCount(),sg.getDelta(),-u1c,ge);
    si.interpolate(
      sg.getCount(),sg.getDelta(),sg.getFirst(),ge,
      sh.getCount(),sh.getDelta(),-uSc,he);
    float[] f = Synthetic.addRickerWavelet(fpeak,fe);
    float[] g = Synthetic.addRickerWavelet(fpeak,ge);
    float[] h = Synthetic.addRickerWavelet(fpeak,he);
    // Sampling for warping
    int n1w = sg.getCount()-(int)(ceil(u1c/sf.getDelta()))+2;
    Sampling s1 = new Sampling(n1w,sf.getDelta(),sf.getFirst());

    // Make min/max shifts
    double s1Min = 0.0;
    double s1Max = u1c+2.0*s1.getDelta();
    double sSMin = 0.0;
    // double sSMax = 0.0;
    double sSMax = uSc+2.0*s1.getDelta();
    print("s1Max="+s1Max+", sSMax="+sSMax);
    DynamicWarpingC3 dw = new DynamicWarpingC3(s1Min,s1Max,sSMin,sSMax,s1);
    Sampling su1 = dw.getSamplingU1();
    Sampling suS = dw.getSamplingUS();
    for (int dg=1; dg<n1w; dg++) {
      int[] g1i = Subsample.subsample(n1w,dg);
      int ng = g1i.length;
      float[] g1 = new float[ng];
      float d1 = (float)s1.getDelta();
      for (int ig=0; ig<ng; ig++)
        g1[ig] = g1i[ig]*d1;
      float[][] u1S = dw.findShifts(sf,f,sg,g,sh,h,g1);
      float[] u1 = u1S[0];
      float[] uS = u1S[1];
      for (int i1=0; i1<u1.length; i1++)
        assert u1[i1]==u1c : "expected "+u1c+", got "+u1[i1];
      for (int iS=0; iS<uS.length; iS++)
        assert uS[iS]==uSc : "expected "+uSc+", got "+uS[iS];
      // print("u1:"); dump(u1S[0]);
      // print("uS:"); dump(u1S[1]);
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private void accumulateAssert(
      int nu1, int nuS, int ne, int ngm, float[][][] d)
  {
    for (int iu1=0; iu1<nu1; iu1++)
      for (int iuS=0; iuS<nuS; iuS++)
        assert d[ngm][iu1][iuS]==ne : "iu1="+iu1+", iuS="+iuS+": "+
          "expected "+ne+" but got "+d[ngm][iu1][iuS];
  }

  private static void print(String s) {
    System.out.println(s);
  }

}