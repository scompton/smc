package test;

import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.ArrayMath.*;

public class MaxSlopeTaperTest {

  public static void main(String[] args) {
    float vpvsMax = 4.0f;
    float vpvsMin = 1.5f;
    int nt = 2501;
    
    float rMax = getSlope(vpvsMax);
    float rMin = getSlope(vpvsMin);
    
    float twoPi = 2.0f*FLT_PI;
    float f = 1.0f/(4.0f*nt);
    float[] r = new float[nt]; 
    for (int i=0; i<nt; i++)
      r[i] = cos(twoPi*f*i);
    normalize(r,rMin,rMax);
    
    for (int i=0; i<nt; i++) {
      assert r[i]<=rMax && r[i]>=rMin;
    }
    
    float[] vpvs = new float[nt];
    for (int i=0; i<nt; i++) {
      vpvs[i] = getVpVs(r[i]);
      assert vpvs[i]<=vpvsMax && vpvs[i]>=vpvsMin;
    }
    
    SimplePlot sp = new SimplePlot();
    sp.addPoints(vpvs);
    sp.setVisible(true);
  }
  
  private static float getSlope(float vpvs) {
    return (vpvs-1.0f)/2.0f;
  }

  private static float getVpVs(float slope) {
    return 1.0f+(slope*2.0f);
  }
  
  private static void normalize(
      float[] f, final float nmin, final float nmax) 
  {
    final int n1 = f.length;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    for (int i1=0; i1<n1; ++i1) {
      float vi = f[i1];
      f[i1] = nrange*(vi-vmin)/range + nmin;
    }
  }
  
  
}
