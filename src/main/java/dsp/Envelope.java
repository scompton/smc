package dsp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Envelope {

  public static float[][][] getEnvelope(final float[][][] f) {
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] e = new float[n3][n2][];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++)
        e[i3][i2] = getEnvelope(f[i3][i2]);
    }});
    return e;
  }

  public static float[][] getEnvelope(final float[][] f) {
    int n2 = f.length;
    final float[][] e = new float[n2][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      e[i2] = getEnvelope(f[i2]);
    }});
    return e;
  }

  public static float[] getEnvelope(float[] f) {
    int n1 = f.length;
    float[] g = copy(f);
    float[] e = new float[n1];
    _htf.apply(n1,f,g);
    for (int i1=0; i1<n1; i1++)
      e[i1] = sqrt(f[i1]*f[i1]+g[i1]*g[i1]);
    return e;
  }

  private static HilbertTransformFilter _htf = new HilbertTransformFilter();

}