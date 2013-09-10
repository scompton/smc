package test;

import static edu.mines.jtk.util.ArrayMath.*;

public class XCorTest {

  /**
   * @param args
   */
  public static void main(String[] args) {
    int n = 100;
    float[] x = new float[n];
    float[] y = new float[n];
//    fill(1.0f,x);
    rand(x);
    rand(y);
    
    float xy = 0, xx = 0 , yy = 0;
    for (int i=0; i<n; i++) {
      xy += x[i]*y[i];
      xx += x[i]*x[i];
      yy += y[i]*y[i];
    }
    float e = 1E-10f;
    float cc1 = xy/(sqrt(xx*yy));
    float cc2 = xy/(sqrt(xx*yy)+e);
    System.out.println("cc1="+cc1+", cc2="+cc2);
  }

}
