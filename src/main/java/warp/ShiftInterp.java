package warp;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.CubicInterpolator;

/*package*/ class ShiftInterp {

  /**
   * Constants for interpolation.
   */
  // /*package*/ enum Method {
  public enum Method {
    LINEAR,
    MONOTONIC,
    SPLINE
  }

  /**
   * Default constructor.
   */
  /*package*/ ShiftInterp() {}

  /**
   * Set the interpolation method for this instance.
   */
  /*package*/ void setMethod(Method method) {
    _method = method;
  }

  /**
   * Interpolates subsampled shifts <code>u[ng]</code> to uniformly sampled
   * shifts <code>ui[n]</code>, where <code>n</code> is the length of the given
   * sampling <code>s</code>.
   * @param s Sampling for interpolated shifts.
   * @param g sparse grid coordinates.
   * @param u sparse shifts.
   * @return the interpolated shifts.
   */
  /*package*/ float[] interpolate(Sampling s, float[] g, float[] u) {
    CubicInterpolator.Method m1;
    switch (_method) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
        _method.toString()+" is not a recognized interpolation method.");
    }
    int ng = g.length;
    int n = s.getCount();
    float[] ui = new float[n];
    CubicInterpolator ci = new CubicInterpolator(m1,ng,g,u);
    double d = s.getDelta();
    double f = s.getFirst();
    double v = f;
    for (int i=0; i<n; i++, v=f+i*d)
      ui[i] = ci.interpolate((float)v);
    return ui;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private Method _method;

}