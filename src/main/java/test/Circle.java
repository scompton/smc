package test;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.SwingUtilities;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;

import static edu.mines.jtk.util.ArrayMath.*;

public class Circle {

  /**
   * @param args
   */
  public static void main(String[] args) {
    int n2 = 1001;
    int n1 = 1001;
    final Sampling s1 = new Sampling(n1,0.1,-50.0);
    final Sampling s2 = new Sampling(n2,0.1,-50.0);
    final float[][] f = new float[n2][n1];

//    float r = 10.0f;
//    float rr = r*r;
//    double[] v2 = s2.getValues();
//    List<Float> xList = new LinkedList<>();
//    List<Float> yList = new LinkedList<>();
//    for (int i2=0; i2<n2; i2++) {
//      double yy = v2[i2]*v2[i2];
//      if (rr-yy<0.0) continue;
//      float y = (float)v2[i2];
//      float x = (float)sqrt(rr-yy);
//      yList.add(y);
//      xList.add(x);
//    }
//    for (int i=yList.size()-2; i>=0; i--) {
//      yList.add( yList.get(i));
//      xList.add(-xList.get(i));
//    }
//
//    int nx = xList.size();
//    final float[] xa = new float[nx];
//    for (int ix=0; ix<nx; ix++) {
//      xa[ix] = xList.get(ix);
//    }
//    
//    int ny = yList.size();
//    final float[] ya = new float[ny];
//    for (int iy=0; iy<ny; iy++) {
//      ya[iy] = yList.get(iy);
//    }

    float r = 10.0000567f;
    float rr = r*r;
    double[] v2 = s2.getValues();
    int hc = (int)ceil(2.0*DBL_PI*r/(s2.getDelta()*0.5))+1;
    double dc = r*2.0/(hc-1);
    System.out.println("hc="+hc+", dc="+dc);
//    int hc = (int)(s2.valueOfNearest(r)*2.0/s2.getDelta())+1;
    int nc = hc*2-1;
    final float[] xa = new float[nc];
    final float[] ya = new float[nc];
    int i = 0;
//    for (int i2=s2.indexOfNearest(-r); i<hc; i++, i2++) {
    for (; i<hc; i++) {
      double y = -r+i*dc;
//      xa[i] = (float)sqrt(rr-v2[i2]*v2[i2]);
//      ya[i] = (float)v2[i2];
      xa[i] = (float)sqrt(rr-y*y);
      ya[i] = (float)y;
    }
    for (int ih=hc-2; ih>=0; ih--, i++) {
      xa[i] = -xa[ih];
      ya[i] =  ya[ih];
    }
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        SimplePlot sp = new SimplePlot();
        sp.addPixels(s1,s2,f);
        sp.addPoints(xa,ya);
        sp.setSize(800,800);
        sp.setVisible(true);
      }
    });
  }

}
