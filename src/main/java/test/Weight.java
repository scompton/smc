package test;

import javax.swing.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import warp.*;

public class Weight {

  public static void main(String[] args) {
    final Sampling s = new Sampling(2001,0.002,0.0);
    final float[] y = VpVsWarp.getU(0.01f,0.07f,s);
    final float[] x = getWeights(10.0f,0.5f,s);
    final float[] xy = mul(x,y);
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
      SimplePlot sp = new SimplePlot();
      // PointsView pv1 = sp.addPoints(s,x);
      PointsView pv2 = sp.addPoints(s,y);
      PointsView pv3 = sp.addPoints(s,xy);
      // pv1.setStyle("k-");
      pv2.setStyle("b-");
      pv3.setStyle("r-");
      sp.setVisible(true);
    }});
  }

  public static float[] getWeights(float a, float p, Sampling s) {
    int n1 = s.getCount();
    float[] x = new float[n1];
    x[0] = 0.0f;
    for (int i1=1; i1<n1; i1++) {
      double v = s.getValue(i1);
      if ((v*5.0f)<a) {
        x[i1] = (float)(pow(v,p)*(3.0f*(pow(5.0f*v/a,2))-2.0f*pow(5.0f*v/a,3)));
      }
      else {
        x[i1] = (float)pow(v,p);
      }
    }
    return x;
  }

}