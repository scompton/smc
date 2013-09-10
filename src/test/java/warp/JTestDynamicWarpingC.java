package warp;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;

import junit.framework.TestCase;
import junit.framework.TestSuite;

public class JTestDynamicWarpingC extends TestCase {

  /*
   * Tests sparse accumulation for the simplest case when rmin and rmax
   * are equal to zero, and all errors are 1.0. This test verifies 
   * returned values for accumulating forward, accumulating reverse,
   * and smoothing (forward + reverse - duplicated error).
   */
  public void testAccumulateSparseSimple() {
    int nl = 50;
    int n1 = 100;
    float[][] e = new float[n1][nl];
    fill(1.0f,e);
    double[] rmin = new double[n1];
    double[] rmax = new double[n1];
    // Test Forward - After normalization, errors should equal 0,
    // at the zero index, and 1 at the last index. Errors should always
    // be increasing from left to right.
    for (int ig=1; ig<=n1; ig++) {
      int[] g = Subsample.subsample(n1,ig);
      int ng = g.length;
      float[][][] dmf = DynamicWarpingC.accumulateForwardSparse(e,rmin,rmax,g);
      float[][] df = dmf[0];
      float[][] mf = dmf[1];
//      Viewer vef = new Viewer(DynamicWarpingC.transposeLag(df),
//          Orientation.X1RIGHT_X2UP);
//      vef.setTitle("Accumulated Errors Forward");
//      vef.setSize(900,600);
//      vef.setColorModel1(ColorMap.JET);
//      vef.addColorBar("Error");
//      vef.show();
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          float dv = df[i1][il];
          assert mf[i1][il]==0.0f;
          assert dv==g[i1]+1 : "i1,il="+i1+","+il+", dv="+dv+"!="+(g[i1]+1);
        }
      }
      // Test Reverse - After normalization, errors should equal 0,
      // at the last index, and 1 at the zero index. Errors should always
      // be increasing from right to left.
      float[][][] dmr = DynamicWarpingC.accumulateReverseSparse(e,rmin,rmax,g);
      float[][] dr = dmr[0];
      float[][] mr = dmr[1];
//      Viewer ver = new Viewer(DynamicWarpingC.transposeLag(dr),
//          Orientation.X1RIGHT_X2UP);
//      ver.setTitle("Accumulated Errors Reverse");
//      ver.setSize(900,600);
//      ver.setColorModel1(ColorMap.JET);
//      ver.addColorBar("Error");
//      ver.show();  
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          float dv = dr[i1][il];
          assert mr[i1][il]==0.0f;
          assert dv==n1-g[i1] : "i1,il="+i1+","+il+", dv="+dv+"!="+(n1-g[i1]);
        }
      }
      // Test Smoothing - For these simple alignment errors, the smoothed
      // errors should be constant everywhere.
      float[][] es = 
          DynamicWarpingC.smoothErrors(e,rmin,rmax,g);
//      Viewer ves = new Viewer(DynamicWarpingC.transposeLag(es),
//          Orientation.X1RIGHT_X2UP);
//      ves.setTitle("Smoothed Errors");
//      ves.setSize(900,600);
//      ves.setColorModel1(ColorMap.JET);
//      ves.addColorBar("Error");
//      ves.show();
      float sc = 1.0f;
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          assert es[i1][il]==sc :
            "es["+i1+"]["+il+"]="+es[i1][il]+", expected "+sc;
        }
      }
    }
  }
  
  /** 
   * This automatically generates a suite of all "test" methods
   * @return new Test
   */
  public static junit.framework.Test suite() {
    return new TestSuite(JTestDynamicWarpingC.class);
  }

  /** 
   * Run all tests with text gui if this class main is invoked
   * @param args command line 
   */
  public static void main (String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
