package test;

import static edu.mines.jtk.util.ArrayMath.*;

import javax.swing.SwingUtilities;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;

public class CosineTaper {

  /**
   * @param args
   */
  public static void main(String[] args) {
    if (args.length!=1) {
      System.out.println("Usage: java CosineTaper max");
      System.exit(0);
    }
    int n = 2501;
    double d = 0.0004;
    double f = -0.5;
    final Sampling s = new Sampling(n,d,f);
    double max = Double.valueOf(args[0]);
    double min = -max;
    double fcs = 0.5/(max-min+d);
    final float[] t = new float[n];
    double twoPi = 2.0*PI;
    for (int i=0; i<n; i++) {
      double v = f+i*d;
      if (v<min || v>max) {
        t[i] = 0.0f;
        continue;
      }
      t[i] = (float)cos(twoPi*fcs*v);
    }

    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        SimplePlot sp = new SimplePlot();
        sp.addPoints(s,t);
        sp.setHLimits(s.getFirst(),s.getLast());
        sp.setVLimits(0.0,1.0);
        sp.setVisible(true);        
      }
    });
  }
}
