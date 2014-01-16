package utils;

import javax.swing.SwingUtilities;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.interp.CubicInterpolator.Method;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Class for numerical integration.
 */
public class Quadrature {

  /**
   * Construct and instance for numerical integration of the function f.
   * @param f function to integrate.
   * @param s sampling of the function f.
   */
  public Quadrature(float[] f, Sampling s) {
    _f = f;
    _n = _f.length;
    float[] sf = new float[_n];
    double[] sv = s.getValues();
    for (int i=0; i<_n; i++)
      sf[i] = (float)sv[i];
    _ci = new CubicInterpolator(Method.LINEAR,sf,_f);
    //_ci = new CubicInterpolator(Method.MONOTONIC,sf,_f);
    //_ci = new CubicInterpolator(Method.SPLINE,sf,_f);
  }
  
  /**
   * Numerical integration using the midpoint method.
   * @param a lower integration limit.
   * @param b upper integration limit.
   * @param j number of subintervals.
   * @return approximate integrated value.
   */
  public double computeMidpoint(double a, double b, int j) {
    Check.argument(j>0,"j>0");
    double h = (b-a)/j;
    //print("h="+h);
    double v = 0.0f;
    double yjm1 = a;
    for (int i=1; i<=j; i++) {
      double yj = yjm1+h;
      double m = (yjm1+yj)/2.0f;
      v += _ci.interpolate((float)m)*h;
      yjm1 = yj;
    }
    return v;
  }

  /**
   * Numerical integration using the trapezoid method.
   * @param a lower integration limit.
   * @param b upper integration limit.
   * @param j number of subintervals.
   * @return approximate integrated value.
   */
  public double computeTrapezoidal(double a, double b, int j) {
    Check.argument(j>0,"j>0");
    double h = (b-a)/j;
    double v = 0.0f;
    double yjm1 = a;
    for (int i=1; i<=j; i++) {
      double yj = yjm1+h;
      double fyjm1 = _ci.interpolate((float)yjm1);
      double fyj = _ci.interpolate((float)yj);
      v += 0.5f*h*(fyjm1+fyj);
      yjm1 = yj;
    }
    return v;
  }

  /**
   * Numerical integration using the Simpson method.
   * @param a lower integration limit.
   * @param b upper integration limit.
   * @param j number of subintervals.
   * @return approximate integrated value.
   */
  public double computeSimpson(double a, double b, int j) {
    Check.argument(j>0,"j>0");
    double h = (b-a)/j;
    double v = 0.0f;
    double yjm1 = a;
    double ms = 2.0f/3.0f;
    double ts = 1.0f/3.0f;
    for (int i=1; i<=j; i++) {
      double yj = yjm1+h;
      double m = (yjm1+yj)/2.0f;
      double vm = _ci.interpolate((float)m);
      double fyjm1 = _ci.interpolate((float)yjm1);
      double fyj = _ci.interpolate((float)yj);
      double vt = 0.5f*(fyjm1+fyj);
      v += (ms*vm+ts*vt)*h;
      yjm1 = yj;
    }
    return v;
  }
  
  public static float[] integrate(float[] f, float d) {
    int n = f.length;
    float[] g = new float[n];
    g[0] = f[0]*d;
    for (int i=1; i<n; i++) {
      g[i] = g[i-1]+f[i]*d;
    }
    return g;
  }
  
  public static void main(String[] args) {
    int n = 501;
    final Sampling s = new Sampling(n,2.0*DBL_PI/(n-1),0.0);
    final float[] f = new float[n];
    for (int i=0; i<n; i++)
      f[i] = (float)cos(s.getValue(i));
    final float[] gi = integrate(f,(float)s.getDelta());
    
    Quadrature quad = new Quadrature(f,s);
    double a = DBL_PI/4.0;
    //double a = 0.0;
    double b = 3.0*DBL_PI/2.0;
    int j = 10;
    double mid = quad.computeMidpoint(a,b,j);
    double trp = quad.computeTrapezoidal(a,b,j);
    double smp = quad.computeSimpson(a,b,j);
    System.out.format("integral=%f, mid=%f, trp=%f, smp=%f%n",
        gi[s.indexOf(b)],mid,trp,smp);

    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        SimplePlot spf = new SimplePlot();
        spf.addPoints(s,f);
        spf.setVLimits(-1.0,1.0);
        spf.setTitle("Cosine Function");
        //SimplePlot spi = new SimplePlot();
        //spi.addPoints(s,gi);
        //spi.setTitle("Sum");
        //spi.setVLimits(-1.0,1.0);
      }
    });
  }
  
  ////////////////////////////////////////////////////////////////////////////
  // Private

  private float[] _f;
  private int _n;
  private CubicInterpolator _ci;
  
  private static void print(String s) {
    System.out.println(s);
  }
  
}
