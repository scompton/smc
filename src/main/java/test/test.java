package test;

import edu.mines.jtk.util.Almost;
import edu.mines.jtk.util.ArrayMath;
import edu.mines.jtk.util.MathPlus;
import edu.mines.jtk.util.Stopwatch;
import static test.Operations.*;
import static test.Util.plot;

public class test {

  public static void main(String[] args) {
    float[][] d = new float[200][200];
    d[3][7] = -0.25f;
    d[4][7] = -0.25f;
    d[5][7] = -0.25f;
    d[5][6] = -0.25f;
    d[4][6] = 1f;
    
    d[5][4] = -0.25f;
    d[6][4] = -0.25f;
    d[7][4] = -0.25f;
    d[7][3] = -0.25f;
    d[6][3] = 1f;
    
    d[15][7] = -1f;
    d[18][15] = 1f;
    
    int n1 = 40000;
    float[] dd = new float[n1];
    float[] yy = new float[n1];
    convert(d,200,200,dd);
    int[] nhf = {3,2};
    int[] n = {200,200};
    HelixFilter h = new HelixFilter(nhf,n,n1,n1);
    h.helixLag(nhf,n);
    float[] a = {1,-0.25f,-0.25f,-0.25f,-0.25f};
    h.coef = a;
    float[] ddo = ArrayMath.copy(dd);
    plot(dd,200,200,"0");
    boolean doParallel = false;
    Stopwatch sw = new Stopwatch(); 
    sw.restart();
    for( int i=0; i<10000; i++) {
      helixDecon(h,dd,dd);
      //plot(dd,20,20,"1");
      helixDeconAdjoint(h,dd,dd);
      //plot(dd,20,20,"2");
      helixConvAdjoint(h,yy,dd,doParallel);
      //plot(yy,20,20,"3");
      helixConv(h,yy,dd,doParallel);
    }
    sw.stop();
    System.out.println("Filter time: "+sw.time()+" seconds");
    plot(dd,200,200,"4");
    float[] diff = ArrayMath.sub(dd,ddo);
    plot(diff,200,200,"4");
    Almost almost = new Almost(MathPlus.FLT_EPSILON*10.0f,0.0001);
    for (int i1=0; i1<n1; i1++) {
      assert almost.zero(diff[i1]) : 
        "Arrays are not equal! diff["+i1+"]="+diff[i1];
    }
  }
  
}
