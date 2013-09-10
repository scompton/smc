package utils;

import java.util.Random;
import javax.swing.SwingUtilities;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.SimplePlot;

import static edu.mines.jtk.util.ArrayMath.*;

public class Synthetic {

  /**
   * Get a synthetic trace with the given Sampling and peak frequency.
   * @param s Sampling of synthetic trace.
   * @param fpeak peak frequency in Hz.
   * @param seed for random events.
   * @return a synthetic trace.
   */
  public static float[] getTrace(Sampling s, float fpeak, long seed) {
    int n = s.getCount();
    double d = s.getDelta();
    fpeak = fpeak*(float)d; // cycles/s * s/sample
    float[] f = makeRandomEvents(n,seed);
    return addRickerWavelet(fpeak,f);
  }

  /**
   * Get random events for a synthetic trace. These are just impulses without
   * a wavelet.
   * @param n lenght of the output array.
   * @param seed for random events.
   * @return an array of random impulse responses.
   */
  public static float[] makeRandomEvents(int n, long seed) {
    Random r = new Random(seed);
    return pow(mul(2.0f,sub(randfloat(r,n),0.5f)),15.0f);
  }

  /**
   * Adds a ricker wavelet with specified peak frequency to array <code>f</code>.
   * @param fpeak peak frequency in Cycles/Sample.
   * @param f array of impulse responses.
   * @return an array of random impulse responses convolved with a ricker
   *   wavelet.
   */
  public static float[] addRickerWavelet(float fpeak, float[] f) {
    int n = f.length;
    int ih = (int)(3.0f/fpeak);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; jh++)
      h[jh] = ricker(fpeak,jh-ih);
    float[] g = new float[n];
    Conv.conv(nh,-ih,h,n,0,f,n,0,g);
    return g;
  }

  public static void main(String[] args) {
    final Sampling s = new Sampling(1001,0.002,0.0);
    float fpeak = 30; // Hz
    long seed = 200;
    final float[] f = Synthetic.getTrace(s,fpeak,seed);
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        SimplePlot.asPoints(s,f);
      }
    });
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private static float ricker(float fpeak, int time) {
    float x = FLT_PI*fpeak*time;
    return (1.0f-2.0f*x*x)*exp(-x*x);
  }

}