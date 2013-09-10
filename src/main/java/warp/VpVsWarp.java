package warp;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.Check;

import static edu.mines.jtk.interp.CubicInterpolator.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class VpVsWarp {

  /**
   * Computes Vp/Vs values at every time sample. The values
   * correspond to a scaled cosine function with a constant
   * frequency b.
   * @param m centered Vp/Vs values of the cosine function.
   * @param a scale factor between -1.0 and 1.0.
   * @param b frequency in Hz.
   * @param s Sampling.
   * @return synthetic Vp/Vs.
   */
  public static float[] getVpVs(float m, float a, float b, Sampling s) {
    Check.argument(a>=-1.0 & a<=1.0,"a>=-1.0 & a<=1.0");
    int n1 = s.getCount();
    double ft = s.getFirst();
    double dt = s.getDelta();
    float[] vpvs = new float[n1];
    for (int i1=0; i1<n1; i1++) {
      double t = ft+i1*dt;
      vpvs[i1] = m + a*(float)cos(2.0*PI*b*t);
    }
    return vpvs;
  }

  public static float[] getU(float a, float b, Sampling s) {
    Check.argument(a>=-1.0 & a<=1.0,"a>=-1.0 & a<=1.0");
    int n1 = s.getCount();
    double ft = s.getFirst();
    double dt = s.getDelta();
    float[] u = new float[n1];
    for (int i1=0; i1<n1; i1++) {
      double t = ft+i1*dt;
      u[i1] = a*(float)sin(2.0*PI*b*t);
    }
    return u;
  }

  /**
   * Computes the shifts u, from the vpvs values.
   * @param vpvs array of computed vpvs values.
   * @param du shift delta.
   * @return shift from synthetic Vp/Vs.
   */
  public static float[] getU(float[] vpvs, float du) {
    int n1 = vpvs.length;
    float[] u = new float[n1];
    u[0] = (vpvs[0]-1.0f)*0.5f*du;
    for (int i1=1; i1<n1; i1++)
      u[i1] = u[i1-1] + (vpvs[i1]-1.0f)*0.5f*du;
    return u;
  }

  public static float[] getU2(
      float[] u1, Sampling s1, float[] uS, Sampling sS)
  {
    int n1 = u1.length;
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();

    int nS = uS.length;
    double fS = sS.getFirst();
    double dS = sS.getDelta();
    float[] xS = new float[nS];
    for (int iS=0; iS<nS; iS++)
      xS[iS] = (float)(fS+iS*dS);
    CubicInterpolator ci = new CubicInterpolator(Method.LINEAR,xS,uS);

    float[] u2 = new float[n1];
    for (int i1=0; i1<n1; i1++)
      u2[i1] = u1[i1]+ci.interpolate((float)(f1+i1*d1)+u1[i1]);
    return u2;
  }

  public static float getAverage(float[] vpvs) {
    int n1 = vpvs.length;
    return sum(vpvs)/n1;
  }

  /**
   * Returns g warped to f, with shifts computed from vpvs array.
   * @param g
   * @param vpvs
   * @param sg
   * @param sf
   * @return g warped to f.
   */
  public static float[] warp(
      float[] g, float[] u, Sampling sg, Sampling sf)
  {
    int ng = g.length;
    int nf = u.length;
    double fg = sg.getFirst();
    double dg = sg.getDelta();
    double ff = sf.getFirst();
    double df = sf.getDelta();
    SincInterp si = new SincInterp();
    float[] f = new float[nf];
    for (int i=0; i<nf; ++i) {
      double y = ff+i*df;
      f[i] = si.interpolate(ng,dg,fg,g,y+u[i]);
    }
    return f;
  }

  public static float[][] shift2(int n2, float[] x, double scale, double f) {
    int n1 = x.length;
    float[][] s = new float[n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    CubicInterpolator ci =
        new CubicInterpolator(CubicInterpolator.Method.LINEAR,r,x);
    for (int i2=0; i2<n2; i2++) {
      float s2 = (float)cos(2.0*PI*f*i2);
      float s2s = (float)(s2*scale);
      assert(abs(s2s)>=0 && abs(s2s)<=scale):"|s2s|="+abs(s2s)+", scale="+scale;
      for (int i1=0; i1<n1; i1++) {
        s[i2][i1] = ci.interpolate(s2s+i1);
      }
    }
    return s;
  }

  public static float[][][] shift3(
      int n3, float[][] x, double scale, double f)
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][][] s = new float[n3][n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci =
            new CubicInterpolator(CubicInterpolator.Method.LINEAR,r,x[i2]);
        float s2 = (float)sin(2.0*PI*f*i2);
        for (int i1=0; i1<n1; i1++) {
          s[i3][i2][i1] = ci.interpolate((float)(s2*scale)+i1);
        }
      }
    }
    return s;
  }

}
