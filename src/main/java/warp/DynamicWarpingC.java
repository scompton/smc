package warp;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.CancellationException;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.dsp.SincInterp.Extrapolation;
import edu.mines.jtk.interp.BicubicInterpolator2;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static warp.WarpUtils.*;

/**
 * Smooth dynamic warping for converted wave images.
 */
public class DynamicWarpingC {

  /**
   * Constants for interpolation.
   */
  public enum Interp {
    LINEAR,
    MONOTONIC,
    SPLINE
  }

  /**
   * Constructor for smooth dynamic warping of converted wave images. Required
   * parameters are minimum and maximum lag. The number of lags is nl =
   * {@code lMax}-{@code lMin}+1. For warping PS to PP images the typical
   * minimum lag is 0, but negative values may need to be allowed due to
   * statics corrections applied to the data. The maximum lag can be computed
   * from the {@link #computeMaxLag(int, float)} method.
   * @param lMin the minimum lag.
   * @param lMax the maximum lag.
   */
  public DynamicWarpingC(int lMin, int lMax) {
    Check.argument(lMin<=lMax,"lMin<=lMax");
    _nl = lMax-lMin+1;
    _sl1 = new Sampling(_nl,1.0,lMin);
    _r1Min =  0.0;
    _r1Max =  1.0;
    _r2Min = -1.0;
    _r2Max =  1.0;
    _r3Min = -1.0;
    _r3Max =  1.0;
    _si = new SincInterp();
    _si.setExtrapolation(Extrapolation.ZERO);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   */
  public void setStrainLimits(double r1Min, double r1Max) {
    setStrainLimits(r1Min,r1Max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
    double r1Min, double r1Max,
    double r2Min, double r2Max)
  {
    setStrainLimits(r1Min,r1Max,r2Min,r2Max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   * @param r3Min lower bound on strain in 3rd dimension.
   * @param r3Max upper bound on strain in 3rd dimension.
   */
  public void setStrainLimits(
    double r1Min, double r1Max,
    double r2Min, double r2Max,
    double r3Min, double r3Max)
  {
    _r1Min = r1Min; _r1Max = r1Max;
    _r2Min = r2Min; _r2Max = r2Max;
    _r3Min = r3Min; _r3Max = r3Max;
  }

  public void setStrainLimits(double[] r1Min, double[] r1Max) {
    setStrainLimits(r1Min,r1Max,-1.0,1.0,-1.0,1.0);
  }

  public void setStrainLimits(
      double[] r1Min, double[] r1Max,
      double r2Min, double r2Max)
  {
    setStrainLimits(r1Min,r1Max,r2Min,r2Max,-1.0,1.0);
  }

  public void setStrainLimits(
      double[] r1Min, double[] r1Max,
      double r2Min, double r2Max,
      double r3Min, double r3Max)
  {
    _r1MinA = r1Min; _r1MaxA = r1Max;
    _r2Min  = r2Min; _r2Max  = r2Max;
    _r3Min  = r3Min; _r3Max  = r3Max;
  }

  /**
   * Set the {@link WarperWorkTracker} implementation that can report and
   * monitor progress for completed task.
   * @param wwt
   */
  public void setWorkTracker(WarperWorkTracker wwt) {
    _wwt = wwt;
  }

  /**
   * Returns the number of lags.
   * @return the number of lags.
   */
  public int getNumberOfLags() {
    return _nl;
  }

  /**
   * Find shifts for 1D traces. For the input trace {@code f[n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng1]}. This length matches the length of array {@code g1}
   * and the indices of the subsampled shifts are specified by contents of
   * this array. The contents of the {@code g1} array are rounded to the
   * nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is
   * {@code ui[n1]}.
   * @param f the PP trace.
   * @param g the PS trace.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(float[] f, float[] g, float[] g1, Interp interp1) {
    return findShifts(f,g,g1,interp1,null);
  }

  /**
   * Find shifts for 1D traces. For the input trace {@code f[n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng1]}. This length matches the length of array {@code g1}
   * and the indices of the subsampled shifts are specified by contents of
   * this array. The contents of the {@code g1} array are rounded to the
   * nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is
   * {@code ui[n1]}.
   * @param f the PP trace.
   * @param g the PS trace.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param xl coordinate array of lag values. These values along with the
   *  {@code x1} coordinate array set a hard constraint in the the error
   *  function. Can be {@code null}.
   * @param x1 coordinate array of i1 values. These values along with the
   *  {@code xl} coordinate array set a hard constraint in the the error
   *  function. Can be {@code null}.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(
      float[] f, float[] g, float[] g1, Interp interp1, int[] l1)
  {
    int n1 = f.length;
    int ng1 = g1.length;
    int[] g1i = new int[ng1];
    for (int ig1=0; ig1<ng1; ig1++)
      g1i[ig1] = (int)(g1[ig1]+0.5f);
    printInfo(n1,1,1,ng1,1,1);

    double[] r1Min = getR1Min(n1);
    double[] r1Max = getR1Max(n1);

    float[][] e = computeErrors(f,g);
    fixShifts(e,l1);
    g1i = KnownShiftUtil.getG1(g1i,l1,r1Min,r1Max);
    // float[][][] dm = accumulateForwardSparse(e,r1Min,r1Max,g1i);
    float[][] es = smoothErrors(e,r1Min,r1Max,g1i);
    float[][][] dm = accumulateForward(es,r1Min,r1Max,g1i);
    float[] u = backtrackReverse(dm[0],dm[1]);

    // Check slopes
    float lastLag = (float)_sl1.getLast();
    for (int ig1=1; ig1<ng1; ig1++) {
      int i1 = g1i[ig1];
      int i1m1 = g1i[ig1-1];
      float n = u[ig1] - u[ig1-1];
      float d = i1-i1m1;
      float r = n/d;
      assert (r>=r1Min[i1] && r<=r1Max[i1]) || u[ig1]==lastLag :
        "n="+n+", d="+d+", r="+r;
    }

    return interpolate(n1,g1,u,interp1);
  }

  /**
   * Find shifts for 2D images. For the input image {@code f[n2][n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays
   * {@code g1,g2} and the indices of the subsampled shifts are specified by
   * contents of these arrays. The contents of the {@code g1} array are rounded
   * to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n2][n1]}.
   * @param f the PP traces.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return shifts for 2D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][] findShifts(
      float[][] f, float[][] g, final float[][] g1, final int[] g2,
      Interp interp1, Interp interp2) throws CancellationException
  {
    return findShifts(f,g,g1,g2,interp1,interp2,new HashMap<Integer,int[]>());
  }

  /**
   * Find shifts for 2D images. For the input image {@code f[n2][n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays
   * {@code g1,g2} and the indices of the subsampled shifts are specified by
   * contents of these arrays. The contents of the {@code g1} array are rounded
   * to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n2][n1]}.
   * @param f the PP traces.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return shifts for 2D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][] findShifts(
      float[][] f, float[][] g, final float[][] g1, final int[] g2,
      Interp interp1, Interp interp2, Map<Integer,int[]> l1Map)
      throws CancellationException
  {
    int n2 = f.length;
    int n1 = f[0].length;
    int ng2 = g2.length;
    int ng1 = g1[0].length;
    Check.argument(n2==g1.length,"n2==g1.length");
    printInfo(n1,n2,1,ng1,ng2,1);

    // Round to integer indices, use float values for interpolation only.
    final int[][] g1i = new int[n2][ng1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        g1i[i2][i1] = (int)(g1[i2][i1]+0.5f);
      }
    }
    final double[] r1Min = getR1Min(n1);
    final double[] r1Max = getR1Max(n1);
    final double[] r2Min = getR2Min(n2);
    final double[] r2Max = getR2Max(n2);

    // Initialize ProgressTracker: smooth1=n2 units, smooth2=ng1 units,
    //   find shifts=ng2 units, interpolation=1 unit.
    int totalWorkUnits = n2+ng1+ng2+1;
    final ProgressTracker pt = new ProgressTracker(totalWorkUnits);
    startProgressThread(pt);

    // Smooth1
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][] es1 = smoothErrors1(f,g,r1Min,r1Max,g1i,l1Map,pt);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    // Smooth2
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,r2Min,r2Max,g2,pt);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);

    // Find shifts on coarse grid
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],r1Min,r1Max,g1i[g2[i2]]);
      u[i2] = backtrackReverse(dm[0],dm[1]);
      pt.worked();
    }});

    // Check slopes
    float lastLag = (float)_sl1.getLast();
    for (int ig2=0; ig2<ng2; ig2++) {
      for (int ig1=1; ig1<ng1; ig1++) {
        int i1 = g1i[g2[ig2]][ig1];
        int i1m1 = g1i[g2[ig2]][ig1-1];
        float n = u[ig2][ig1] - u[ig2][ig1-1];
        float d = i1-i1m1;
        float r = n/d;
        assert (r>=r1Min[i1] && r<=r1Max[i1]) || u[ig2][ig1]==lastLag :
          "n="+n+", d="+d+", r="+r;
      }
    }

    // Interpolate shifts to original grid.
    float[][] ui = interpolate(n1,n2,g1,g2,u,interp1,interp2);
    pt.worked();
    assert pt.getPercentComplete()==100.0f;
    return ui;
  }

  /**
   * Find shifts for 3D images. For the input image {@code f[n3][n2][n1]},
   * shifts {@code u} are computed on a subsampled grid such that the computed
   * shifts are {@code u[ng3][ng2][ng1]}. These lengths match the length of
   * arrays {@code g1,g2,g3} and the indices of the subsampled shifts are
   * specified by contents of these arrays. The contents of the {@code g1}
   * array are rounded to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n3][n2][ne1]}.
   * @param f the PP traces.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param g3 array of size [ng3] specifying third dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return shifts for 3D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][][] findShifts(
      float[][][] f, float[][][] g,
      final float[][][] g1, final int[] g2, final int[] g3,
      Interp interp1, Interp interp23) throws CancellationException
  {
    return findShifts(f,g,g1,g2,g3,interp1,interp23,
        new HashMap<Integer,int[]>());
  }

  public float[][][] findShifts(
      float[][][] f, float[][][] g,
      final float[][][] g1, final int[] g2, final int[] g3,
      Interp interp1, Interp interp23, Map<Integer,int[]> l1Map)
      throws CancellationException
  {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    final int ng3 = g3.length;
    final int ng2 = g2.length;
    final int ng1 = g1[0][0].length;
    Check.argument(n3==g1.length,"n3==g1.length");
    Check.argument(n2==g1[0].length,"n2==g1[0].length");
    printInfo(n1,n2,n3,ng1,ng2,ng3);

    // Round to integer indices, use float values for interpolation only.
    final int[][][] g1i = new int[n3][n2][ng1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          g1i[i3][i2][i1] = (int)(g1[i3][i2][i1]+0.5f);
        }
      }
    }
    final double[] r1Min = getR1Min(n1);
    final double[] r1Max = getR1Max(n1);
    final double[] r2Min = getR2Min(n2);
    final double[] r2Max = getR2Max(n2);
    final double[] r3Min = getR3Min(n3);
    final double[] r3Max = getR3Max(n3);

    // Initialize ProgressTracker: smooth1=n2*n3 units, smooth2=n3*ng1 units,
    //   smooth3=ng1*ng2 units, find shifts=ng2*ng3 units, interpolation=1 unit.
    int totalWorkUnits = (n2*n3)+(n3*ng1)+(ng1*ng2)+(ng2*ng3)+1;
    final ProgressTracker pt = new ProgressTracker(totalWorkUnits);
    startProgressThread(pt);

    // Smooth1
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][][] es1 = smoothErrors1(f,g,r1Min,r1Max,g1i,l1Map,pt);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    // Smooth2
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,r2Min,r2Max,g2,pt);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);

    // Smooth3
    s.restart();
    print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,r3Min,r3Max,g3,pt);
    normalizeErrors(es);
    print("Finished 3rd dimension smoothing in "+s.time()+" seconds");

    // Find shifts on coarse grid
    final float[][][] u = new float[ng3][ng2][];
    int ng23 = ng3*ng2;
    Parallel.loop(ng23,new Parallel.LoopInt() {
    public void compute(int i23) {
      int i2 = i23%ng2;
      int i3 = i23/ng2;
      float[][][] dm = accumulateForward(
          es[i3][i2],r1Min,r1Max,g1i[g3[i3]][g2[i2]]);
      u[i3][i2] = backtrackReverse(dm[0],dm[1]);
      pt.worked();
    }});

    // Check slopes
    float lastLag = (float)_sl1.getLast();
    for (int ig3=0; ig3<ng3; ig3++) {
      for (int ig2=0; ig2<ng2; ig2++) {
        for (int ig1=1; ig1<ng1; ig1++) {
          int i1 = g1i[g3[ig3]][g2[ig2]][ig1];
          int i1m1 = g1i[g3[ig3]][g2[ig2]][ig1-1];
          float n = u[ig3][ig2][ig1] - u[ig3][ig2][ig1-1];
          float d = i1-i1m1;
          float r = n/d;
          assert (r>=r1Min[i1] && r<=r1Max[i1]) || u[ig3][ig2][ig1]==lastLag :
            "n="+n+", d="+d+", r="+r;
        }
      }
    }

    // Interpolate
    float[][][] ui = interpolate(n1,n2,n3,g1,g2,g3,u,interp1,interp23);
    pt.worked();
    assert pt.getPercentComplete()==100.0f;
    return ui;
  }

  /**
   * Applies the shifts {@code u} to the PS trace {@code g}.
   * @param n1f the length of PP trace. This is the length of the returned
   *  warped trace.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS trace to the PP trace.
   * @return the warped PS trace.
   */
  public float[] applyShifts(int n1f, float[] g, float[] u) {
    int n1g = g.length;
    int nu = u.length;
    int num = nu-1;
    float[] h = new float[n1f];
    for (int iu=0; iu<nu; ++iu) {
      h[iu] = _si.interpolate(n1g,1.0,0.0,g,iu+u[iu]);
    }
    for (int i1=nu; i1<n1f; ++i1) {
      h[i1] = _si.interpolate(n1g,1.0,0.0,g,i1+u[num]);
    }
    return h;
  }
  
  /**
   * Applies the shifts {@code u} to the PS image {@code g}.
   * @param n1f the length of PP traces. This is the length of the returned
   *  warped traces.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][] applyShifts(
      final int n1f, final float[][] g, final float[][] u)
  {
    final int n2 = g.length;
    final int n1g = g[0].length;
    final int n1u = u[0].length;
    final int n1um = n1u-1;
    final float[][] hf = new float[n2][n1f];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i2],i1+u[i2][i1]);
      }
      for (int i1=n1u; i1<n1f; ++i1) {
        hf[i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i2],i1+u[i2][n1um]);
      }
    }});
    return hf;
  }

  /**
   * Applies the shifts {@code u} to the PS image {@code g}.
   * @param n1f the length of PP traces. This is the length of the returned
   *  warped traces.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][][] applyShifts(
      final int n1f, final float[][][] g, final float[][][] u)
  {
    final int n3 = g.length;
    final int n2 = g[0].length;
    final int n1g = g[0][0].length;
    final int n1u = u[0][0].length;
    final int n1um = n1u-1;
    final float[][][] hf = new float[n3][n2][n1f];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1u; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i3][i2],
              i1+u[i3][i2][i1]);
        }
        for (int i1=n1u; i1<n1f; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i3][i2],
              i1+u[i3][i2][n1um]);
        }
      }
    }});
    return hf;
  }

  /**
   * Compute the maximum length of PP traces from a guess of the average Vp/Vs
   * ratio. This method can be used to find the maximum PP length useful for
   * warping. Truncating the PP data to this length will improve efficiency.
   * This calculation uses the entire length of the PS trace {@code n1PS} and
   * determines the maximum sample in the PP trace that could correspond to the
   * final PS sample. The length is [2.0/({@code vpvsAvg}+1.0)]*{@code n1PS} or
   * {@code n1PP} if the computed value is larger than {@code n1PP}.
   * @param n1PP the length of PP traces.
   * @param n1PS the length of PS traces.
   * @param vpvsAvg a guess of the average Vp/Vs ratio at the PP sample that
   *  corresponds with the maximum PS sample. Usually a value of 2.0 is a good
   *  starting point.
   * @return the maximum length of PP traces useful for warping, based on the
   *  {@code vpvsAvg} and {@code n1PS}.
   */
  public static int computeMaxLength(int n1PP, int n1PS, float vpvsAvg) {
    float scale = getScale(vpvsAvg);
    int n1PPMax = (int)ceil(n1PS*scale);
    return (n1PPMax>n1PP)?n1PP:n1PPMax;
  }

  /**
   * Computes the maximum lag for warping PS to PP traces from a guess of the
   * average Vp/Vs ratio. This value can be used to construct a
   * {@link DynamicWarpingC} instance. This method uses the approach described
   * in the {@link #computeMaxLength(int, int, float)} method. The maximum lag
   * is simply the difference between the computed PP trace length useful for
   * warping and the PS trace length.
   * @param n1PS the length of PS traces.
   * @param vpvsAvg a guess of the average Vp/Vs ratio at the PP sample that
   *  corresponds with the maximum PS sample. Usually a value of 2.0 is a good
   *  starting point.
   * @return the maximum lag (in samples) for warping, based on the
   *  {@code vpvsAvg} and {@code n1PS}.
   */
  public static int computeMaxLag(int n1PS, float vpvsAvg) {
    float scale = getScale(vpvsAvg);
    int n1PPMax = (int)ceil(n1PS*scale);
    int lMax = n1PS-n1PPMax;
    Check.argument(lMax>0,"lMax>0");
    return lMax;
  }

  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications

  /**
   * Find shifts for 2D images from averaged alignment errors. These
   * shifts are constant for all traces.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param g1 array of subsampled indices.
   * @return shifts for 2D images from averaged alignment errors.
   */
  public float[] findShifts2(float[][] f, float[][] g, int[] g1) {
    final float[][] e = computeErrorsSum2(f,g);
    int n1 = f[0].length;
    double[] r1Min = getR1Min(n1);
    double[] r1Max = getR1Max(n1);
    float[][][] dm = accumulateForwardSparse(e,r1Min,r1Max,g1);
    return backtrackReverse(dm[0],dm[1]);
  }

  public float[] findShifts3(float[][][] f, float[][][] g, int[] g1) {
    float[][] e3Avg = computeErrorsSum3(f,g);
    int n1 = f[0][0].length;
    double[] r1Min = getR1Min(n1);
    double[] r1Max = getR1Max(n1);
    float[][][] dm = accumulateForwardSparse(e3Avg,r1Min,r1Max,g1);
    return backtrackReverse(dm[0],dm[1]);
  }

  public float[][] applyShifts(
      float[] u1, float[] uS, float[] ps1, float[] ps2)
  {
    int n1 = ps1.length;
    float[] ps1w = new float[n1];
    float[] ps2w = new float[n1];
    int nu = u1.length;
    int num = nu-1;
    for (int iu=0; iu<nu; ++iu) {
      ps1w[iu] = _si.interpolate(n1,1.0,0.0,ps1,iu+u1[iu]);
      ps2w[iu] = _si.interpolate(n1,1.0,0.0,ps2,iu+u1[iu]+uS[iu]);
    }
    for (int i1=nu; i1<n1; ++i1) {
      ps1w[i1] = _si.interpolate(n1,1.0,0.0,ps1,i1+u1[num]);
      ps2w[i1] = _si.interpolate(n1,1.0,0.0,ps2,i1+u1[num]+uS[num]);
    }
    return new float[][]{ps1w, ps2w};
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u
   * using a backward difference approximation.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values.
   */
  public static float[] vpvsBd(float[] u) {
    int n = u.length;
    float[] vpvs = new float[n];
    vpvs[0] = 1.0f + 2.0f*(u[1]-u[0]); // at i1=0, forward difference
    for (int i1=1; i1<n; ++i1)
      vpvs[i1] = 1.0f + 2.0f*(u[i1]-u[i1-1]);
    return vpvs;
  }

  /**
   * Compute alignment errors for 1D traces.
   * @return alignment errors for 1D traces.
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    float[][] e = new float[n1][_nl];
    computeErrors(f,g,e);
    return e;
  }

  /**
   * Compute alignment errors for 2D traces.
   * @return alignment errors for 2D traces.
   */
  public float[][][] computeErrors2(final float[][] f, final float[][] g) {
    int n2 = f.length;
    int n1 = f[0].length;
    final float[][][] e = new float[n2][n1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(f[i2],g[i2],e[i2]);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 2D images.
   * @return summed alignment errors for 2D images.
   */
  public float[][] computeErrorsSum2(final float[][] f, final float[][] g) {
    final int n2 = f.length;
    final int n1 = f[0].length;
    float[][] e = Parallel.reduce(n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[n1][_nl];
      computeErrors(f[i2],g[i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 3D images.
   * @return summed alignment errors for 3D images.
   */
  public float[][] computeErrorsSum3(
      final float[][][] f, final float[][][] g)
  {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e = new float[n1][_nl];
      computeErrors(f[i3][i2],g[i3][i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  public void fixShifts(float[][][] e, Map<Integer,int[]> l1Map) {
    for (Integer i : l1Map.keySet()) {
      fixShifts(e[i],l1Map.get(i));
    }
  }

  public void fixShifts(float[][] e, int[] l1) {
    if (l1==null) return;
    int n1 = e.length;
    int nl = e[0].length;
//    int ob = 0; // out of bounds
    for (int i=0; i<l1.length; i+=2) {
      int lag = l1[i  ];
      int  i1 = l1[i+1];
      if (i1<0 || i1>=n1 || lag<0 || lag>=nl) {
//        ob++;
        continue;
      }
      for (int il=0; il<nl; il++)
        e[i1][il] = Float.MAX_VALUE;
      e[i1][lag] = 0.0f;
    }
//    print("Fixed shifts complete: nx="+nc+", out of bounds count="+ob);
  }

  public float[][][] accumulateForward(
      float[][] e, double[] rMin, double[] rMax, int[] g)
  {
    float[][] d = new float[e.length][e[0].length];
    float[][] m = new float[e.length][e[0].length];
    accumulateFromSparse(1,e,d,m,g,rMin,rMax);
    return new float[][][]{d,m};
  }

  public static float[][][] accumulateForwardSparse(
      float[][] e, double[] rmin, double[] rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }

  public static float[][][] accumulateReverseSparse(
      float[][] e, double[] rmin, double[] rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(-1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }

  public float[] backtrackReverse(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(-1,_sl1,d,m,u);
    return u;
  }

  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(1,_sl1,d,m,u);
    return u;
  }

  public static float[][][] shiftVolume(
      float[] u1, float[] uS, int fr, float[][][] e) {
    int n1 = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;
    Check.argument(n1==u1.length, "n1==u1.length");
    Check.argument(n1==uS.length, "n1==uS.length");
    float[][][] sv = zerofloat(nlS,nl1,n1);
    for (int i1=0; i1<n1; ++i1) {
      int il1 = (int)u1[i1];
      int ilS = (int)uS[i1]*fr;
      sv[i1][il1][ilS] = 1.0f;
    }
    return sv;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _r1Min, _r1Max;
  private double _r2Min, _r2Max;
  private double _r3Min, _r3Max;
  private double[] _r1MinA, _r1MaxA;
  private int _nl; // number of lags
  private Sampling _sl1; // sampling of pp to ps1 lags
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
  private WarperWorkTracker _wwt;

  private static final float BYTES_TO_MB = 1.0f/1000000.0f;

  /**
   * Computes scale to apply to PS traces. n1PP = scale*n1PS
   * @param vpvsAvg
   * @return a scaler to apply to PS traces.
   */
  private static float getScale(float vpvsAvg) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    return 2.0f/(vpvsAvg+1.0f);
  }

  private void printInfo(
      int n1, int n2, int n3, int ng1, int ng2, int ng3)
  {
    float e1Mem = (float)ng1*n2*n3*_nl*4.0f*BYTES_TO_MB;
    float e2Mem = (float)ng1*ng2*n3*_nl*4.0f*BYTES_TO_MB;
    float e3Mem = (float)ng1*ng2*ng3*_nl*4.0f*BYTES_TO_MB;
    print("DynamicWarpingC info:");
    print("  Input data samples (n1,n2,n3): ("+n1+","+n2+","+n3+")");
    print("  Coarse grid samples: (ng1,ng2,ng3): ("+ng1+","+ng2+","+ng3+")");
    print("  Number of lags: "+_nl);
    print("  Alignment error smooth 1 memory: "+e1Mem+" MB");
    print("  Alignment error smooth 2 memory: "+((n2>1)?(e2Mem+" MB"):"NA"));
    print("  Alignment error smooth 3 memory: "+((n3>1)?(e3Mem+" MB"):"NA"));
  }

  private double[] getR1Min(int n1) {
    double[] r1Min;
    if (_r1MinA==null)
      r1Min = filldouble(_r1Min,n1);
    else
      r1Min = _r1MinA;
    return r1Min;
  }

  private double[] getR1Max(int n1) {
    double[] r1Max;
    if (_r1MaxA==null)
      r1Max = filldouble(_r1Max,n1);
    else
      r1Max = _r1MaxA;
    return r1Max;
  }

  private double[] getR2Min(int n2) {
    return filldouble(_r2Min,n2);
  }

  private double[] getR2Max(int n2) {
    return filldouble(_r2Max,n2);
  }

  private double[] getR3Min(int n3) {
    return filldouble(_r3Min,n3);
  }

  private double[] getR3Max(int n3) {
    return filldouble(_r3Max,n3);
  }

  /**
   * Starts a thread to monitor progress if a {@link #_wwt} instance was set
   * from the {@link #setWorkTracker(WarperWorkTracker)} method.
   * @param pt the {@link ProgressTracker} instance that reports completed
   *   work.
   */
  private void startProgressThread(final ProgressTracker pt) {
    if (_wwt!=null) {
      Thread thread = new Thread(new Runnable() {
        @Override
        public void run() {
          _wwt.setTotalWorkUnits(pt._totalWorkUnits);
          int wus = pt.getCompletedWorkUnits();
          while (pt.getPercentComplete()!=100.0f) {
            try {
              if (_wwt.isCanceled()) {
                pt.setCanceled(true);
                break;
              }
              int wu = pt.getCompletedWorkUnits();
              _wwt.worked(wu-wus);
              wus = wu;
              Thread.sleep(50);
            } catch (InterruptedException e) {
              Thread.currentThread().interrupt();
            }
          }
        }
      });
      thread.start();
    }
  }

  /**
   * Interpolates subsampled shifts u[ng1] to uniformly sampled shifts
   * ui[_ne1].
   * @param g1 sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for linear interpolation,
   *  {@code false} for cubic.
   * @return the interpolated shifts.
   */
  public float[] interpolate(int n1, float[] g1, float[] u, Interp interp1) {
    int ng1 = g1.length;
    float[] ui = new float[n1];
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }
    CubicInterpolator ci = new CubicInterpolator(m1,ng1,g1,u);
    for (int i=0; i<n1; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }

  /**
   * Interpolates subsampled shifts u[ng2][ng1] to uniformly sampled
   * shifts ui[_n2][_ne1]. The locations of the subsampled shifts u
   * are defined by the input arrays g1 and g2. The g2 array is
   * assumed to contain integer indices consistent for all n1 locations.
   * That is, g2 indices are the same at every i1, for i1=0,...,n1-1.
   * </p>
   * The g1 coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, for all n2 indices.
   * </p>
   * The interpolation is done in two passes. First interpolation
   * is done in the second dimension where subsampling must be
   * regular. The second pass does interpolation in the first
   * dimension where subsampling may be irregular.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return the interpolated shifts ui[_n2][_ne1].
   */
  private float[][] interpolate(
      int n1, int n2, float[][] g1, int[] g2, float[][] u,
      Interp interp1, Interp interp2)
  {
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    float[] g2f = new float[ng2];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = g2[ig2];

    float[][] ui = new float[n2][n1];
    CubicInterpolator.Method m2;
    switch (interp2) {
      case LINEAR:    m2 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m2 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m2 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp2.toString()+" is not a recognized interpolation method.");
    }
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }

    // interpolate in the second dimension.
    float[] u2 = new float[ng2];
    float[][] ui2 = new float[n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++)
        u2[i2] = u[i2][i1];
      CubicInterpolator ciu = new CubicInterpolator(m2,g2f,u2);
      for (int i2=0; i2<n2; i2++)
        ui2[i2][i1] = ciu.interpolate(i2);
    }

    // interpolate in the first dimension.
    for (int i2=0; i2<n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(m1,g1[i2],ui2[i2]);
      for (int i1=0; i1<n1; i1++)
        ui[i2][i1] = ci.interpolate(i1);
    }

    return ui;
  }

  /**
   * Interpolates subsampled shifts u[ng3][ng2][ng1] to uniformly sampled
   * shifts ui[_n3][_n2][_ne1]. The locations of the subsampled shifts u are
   * defined by the input arrays g1, g2, and g3. The indices in the g3 array
   * must be the same at every i2, for i2=0,...,n2-1 and i1, for i1=0,...n1-1.
   * The indices in the g2 array must be the same at every i3, for
   * i3=0,...,n3-1, and i1, for i1=0,...,n1-1.
   * </p>
   * The g1 coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, for all n3 and n2 indices.
   * </p>
   * The interpolation is done in two passes. First interpolation
   * is done in the second and third dimension where subsampling
   * must be regular. The second pass does interpolation in the
   * first dimension where subsampling may be irregular.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid indices.
   * @param g3 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return the interpolated shifts ui[_n3][_n2][_ne1].
   */
  private float[][][] interpolate(
      int n1, int n2, int n3, float[][][] g1, int[] g2, int[] g3, float[][][] u, 
      Interp interp1, Interp interp23)
  {
    int ng1 = g1[0][0].length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng3+"=="+u.length);
    float[] g2f = new float[ng2];
    float[] g3f = new float[ng3];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = g2[ig2];
    for (int ig3=0; ig3<ng3; ig3++)
      g3f[ig3] = g3[ig3];

    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }

    BicubicInterpolator2.Method m23 = null;
    boolean doLinear;
    switch (interp23) {
      case LINEAR:    doLinear = true;
                      break;
      case MONOTONIC: doLinear = false;
                      m23 = BicubicInterpolator2.Method.MONOTONIC;
                      break;
      case SPLINE:    doLinear = false;
                      m23 = BicubicInterpolator2.Method.SPLINE;
                      break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }

    // interpolate the second and third dimension.
    float[][] u23 = new float[ng3][ng2];
    float[][][] ui23 = new float[n3][n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++)
        for (int i2=0; i2<ng2; i2++)
          u23[i3][i2] = u[i3][i2][i1];
      if (doLinear) {
        BilinearInterpolator2 bli = new BilinearInterpolator2(g2f,g3f,u23);
        for (int i3=0; i3<n3; i3++)
          for (int i2=0; i2<n2; i2++)
            ui23[i3][i2][i1] = bli.interpolate(i2,i3);
      } else {
        BicubicInterpolator2 bci =
            new BicubicInterpolator2(m23,m23,g2f,g3f,u23);
        for (int i3=0; i3<n3; i3++)
          for (int i2=0; i2<n2; i2++)
            ui23[i3][i2][i1] = bci.interpolate(i2,i3);
      }
    }

    float[][][] ui = new float[n3][n2][n1];
    // interpolate the first dimension
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci =
            new CubicInterpolator(m1,g1[i3][i2],ui23[i3][i2]);
        for (int i1=0; i1<n1; i1++)
          ui[i3][i2][i1] = ci.interpolate(i1);
      }
    }
    return ui;
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  /**
   * Computes alignment errors for {@code f} and {@code g}.
   * @param f
   * @param g
   * @param e
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1Max = e.length;
    int nl = e[0].length;
    int ng = g.length;
    float[] gi = new float[n1Max];
    for (int il=0; il<nl; il++) {
      _si.interpolate(ng,1.0,0.0,g,n1Max,1.0,_sl1.getValue(il),gi);
      for (int i1=0; i1<n1Max; i1++) {
        e[i1][il] = error(f[i1],gi[i1]);
      }
    }
  }

  private void computeErrors(
      float[] f, float[] g, float[] h, float[][][] e)
  {
    int n1max = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;

    // Compute errors for f, g, and h.
    for (int i1=0; i1<n1max; i1++) {
      for (int il1=0,j1=i1; il1<nl1; il1++,j1++) {
        float e1 = error(f[i1],g[j1]);
        for (int ilS=0,jS=j1; ilS<nlS; ilS++,jS++) {
          float e2 = error(f[i1],h[jS]);
          float e3 = error(g[j1],h[jS]);
          e[i1][il1][ilS] = e1+e2+e3;
        }
      }
    }
  }

  /**
   * Returns smooth alignment errors on the sparse grid defined
   * by the indices of g.
   * @param e 2D array of alignment errors.
   * @param rmin minimum slope.
   * @param rmax maximum slope.
   * @param g sparse grid indices.
   * @return smoothed alignment errors with size
   *  [g.length][e[0].length].
   */
  public static float[][] smoothErrors(
      float[][] e, double[] rmin, double[] rmax, int[] g)
  {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    accumulateSparse( 1,rmin,rmax,g,e,ef,null);
    accumulateSparse(-1,rmin,rmax,g,e,er,null);
    float scale = 1.0f/e.length;
    for (int i1=0; i1<ng; i1++) {
      for (int il=0; il<nel; il++) {
        float v = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
        es[i1][il] = Float.isInfinite(v) ? Float.MAX_VALUE : v;
      }
    }
    return es;
  }

  public float[][][] getSmoothErrors(float[][] e, int[] g, int n1) {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    double[] r1Min = getR1Min(n1);
    double[] r1Max = getR1Max(n1);
    accumulateSparse( 1,r1Min,r1Max,g,e,ef,null);
    accumulateSparse(-1,r1Min,r1Max,g,e,er,null);
    float scale = 1.0f/e.length;
    for (int i1=0; i1<ng; i1++) {
      for (int il=0; il<nel; il++) {
        float v = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
        es[i1][il] = Float.isInfinite(v) ? Float.MAX_VALUE : v;
      }
    }
    return new float[][][] {ef,er,es};
  }

  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n2 (size[_n2][])
   * @return smoothed alignment errors with size
   *  [_n2][g1.length][_nel].
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  private float[][][] smoothErrors1(
      final float[][] pp, final float[][] ps,
      final double[] r1min, final double[] r1max, final int[][] g1,
      final Map<Integer,int[]> l1Map, final ProgressTracker pt)
      throws CancellationException
  {
    final int n2 = pp.length;
    final int n1 = pp[0].length;
    final int ng1 = g1[0].length;
    final float[][][] es1 = new float[n2][ng1][_nl];
    final Parallel.Unsafe<float[][]> eu = new Parallel.Unsafe<>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      if (pt.getCanceled()) {
        Thread.currentThread().interrupt();
        throw new CancellationException();
      }
      float[][] e = eu.get();
      if (e==null) eu.set(e=new float[n1][_nl]);
      computeErrors(pp[i2],ps[i2],e);
      int[] l1 = l1Map.get(i2);
      fixShifts(e,l1);
      int[] g1ks = KnownShiftUtil.getG1(g1[i2],l1,r1min,r1max);
      es1[i2] = smoothErrors(e,r1min,r1max,g1ks);
      pt.worked();
    }});
    return es1;
  }

  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second and third dimension. Alignment
   * errors are computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n3 and _n2 (size[_n3][_n2][]).
   * @return smoothed alignment errors with size
   *  [_n3][_n2][g1.length][_nel].
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  private float[][][][] smoothErrors1(
      final float[][][] pp, final float[][][] ps,
      final double[] r1min, final double[] r1max, final int[][][] g1,
      final Map<Integer,int[]> ilMap, final ProgressTracker pt)
      throws CancellationException
  {
    final int n3 = pp.length;
    final int n2 = pp[0].length;
    final int n1 = pp[0][0].length;
    final int ng1 = g1[0][0].length;
    final float[][][][] es1 = new float[n3][n2][ng1][_nl];
    final Parallel.Unsafe<float[][]> eu = new Parallel.Unsafe<>();
    int n23 = n2*n3;
    Parallel.loop(n23,new Parallel.LoopInt() {
    public void compute(int i23) {
      if (pt.getCanceled()) {
        Thread.currentThread().interrupt();
        throw new CancellationException();
      }
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e = eu.get();
      if (e==null) eu.set(e=new float[n1][_nl]);
      computeErrors(pp[i3][i2],ps[i3][i2],e);
      int[] l1 = ilMap.get(i23);
      fixShifts(e,l1);
      int[] g1ks = KnownShiftUtil.getG1(g1[i3][i2],l1,r1min,r1max);
      es1[i3][i2] = smoothErrors(e,r1min,r1max,g1ks);
      pt.worked();
    }});
    return es1;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension.
   * @param e alignment errors.
   * @param r2min minimum slope in the second dimension.
   * @param r2max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size
   *  [g2.length][e[0].length][e[0][0].length].
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  private static float[][][] smoothErrors2(
      final float[][][] e,
      final double[] r2min, final double[] r2max, final int[] g2,
      final ProgressTracker pt) throws CancellationException
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      if (pt.getCanceled()) {
        Thread.currentThread().interrupt();
        throw new CancellationException();
      }
      float[][]  e2 = new float[n2][nl]; // errors at index i1
      for (int i2=0; i2<n2; ++i2)
        e2[i2] = e[i2][i1];
      float[][] es2 = smoothErrors(e2,r2min,r2max,g2);
      for (int i2=0; i2<ng2; i2++)
        es[i2][i1] = es2[i2];
      pt.worked();
    }});
    return es;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension and third dimension.
   * @param e alignment errors.
   * @param r2Min minimum slope in the second dimension.
   * @param r2Max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size
   *  [e.length][g2.length][e[0][0].length][e[0][0][0].length].
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  private static float[][][][] smoothErrors2(
      final float[][][][] e,
      final double[] r2Min, final double[] r2Max, final int[] g2,
      final ProgressTracker pt) throws CancellationException
  {
    final int n3 = e.length;
    final float[][][][] es = new float[n3][][][]; // smoothed errors
    for (int i3=0; i3<n3; i3++)
      es[i3] = smoothErrors2(e[i3],r2Min,r2Max,g2,pt);
    return es;
  }

  /**
   * Returns alignment errors smoothed in the third dimension.
   * Returned errors are sparse in the third dimension, and
   * unchanged in the first and second dimension.
   * @param e alignment errors.
   * @param r3Min minimum slope in the third dimension.
   * @param r3Max maximum slope in the third dimension.
   * @param g3 third dimension sparse grid indices.
   * @return smoothed alignment errors with size
   *  [g3.length][e[0].length][e[0][0].length][e[0][0][0].length].
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  private static float[][][][] smoothErrors3(
      final float[][][][] e,
      final double[] r3Min, final double[] r3Max, final int[] g3,
      final ProgressTracker pt) throws CancellationException
  {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final int ng3 = g3.length;
    final float[][][][] es = new float[ng3][n2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<n2; i2++) {
        if (pt.getCanceled()) {
          Thread.currentThread().interrupt();
          throw new CancellationException();
        }
        float[][]  e3 = new float[n3][nl]; // smooth errors at index i1,i2
        for (int i3=0; i3<n3; i3++)
          e3[i3] = e[i3][i2][i1];
        float[][] es3 = smoothErrors(e3,r3Min,r3Max,g3);
        for (int i3=0; i3<ng3; i3++)
          es[i3][i2][i1] = es3[i3];
        pt.worked();
      }
    }});
    return es;
  }

  private static void accumulateSparse(
      int dir, double[] rMin, double[] rMax, int[] g,
      float[][] e, float[][] d, float[][] m)
  {
    int nl   = e[0].length;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ibg = (dir>0)?0:ngm1; // beginning index
    int ieg = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int isp = ibg; // sparse grid index
    int ie = g[isp]; // error index

    // Initialize accumulation values
    for (int il=0; il<nl; ++il)
      d[isp][il] = e[ie][il];
    isp += is;
    // Loop over all sparse grid points.
    for (; isp!=ieg; isp+=is) {
      int ispm1 = isp-is; // previous sparse grid index
      ie = g[isp]; // new error index
      int je = g[ispm1]; // min error index for interpolation.
      int dg = ie-je; // sparse grid delta, possibly negative
      int kmin, kmax;
      if (dg>0) { // forward
        kmin = (int) ceil(-rMax[ie]*dg);
        kmax = (int)floor(-rMin[ie]*dg);
      } else { // reverse
        kmin = (int) ceil(-rMin[ie]*dg);
        kmax = (int)floor(-rMax[ie]*dg);
      }
      kmin = kmin>kmax ? kmax : kmin;
      float[] dm = new float[nl];
      fill(Float.MAX_VALUE,d[isp]);
      // loop over all slope indices
      for (int k=kmin; k<=kmax; k++) {
        int ils = max(0,-k);
        int ile = min(nl,nl-k);
        for (int il=ils; il<ile; il++)
          dm[il] = d[ispm1][il+k] + e[ie][il];
        float r = (float)k/(float)dg; // slope
        if (r==0) { // zero slope, no interpolation necessary
          for (int x=je+is; x!=ie; x+=is)
            for (int il=ils; il<ile; il++)
              dm[il] += e[x][il];
        } else { // linearly interpolate
          for (int x=je+is; x!=ie; x+=is) {
            float ky = r*(ie-x);
            int k1 = (int)ky;
            if (ky<0.0f) --k1;
            int k2 = k1+1;
            float w1 = k2-ky;
            float w2 = 1.0f-w1;
            for (int il=ils; il<ile; il++)
              dm[il] += w1*e[x][k1+il]+w2*e[x][k2+il];
          }
        }
        // update previous errors and record moves.
        for (int il=ils; il<ile; il++) {
          if (dm[il]<d[isp][il]) {
            d[isp][il] = dm[il];
            if (m!=null)
              m[ispm1][il] = k;
          }
        }
      }
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   * @param m output array[ni][nl] of recorded moves.
   * @param g
   * @param rMin
   * @param rMax
   */
  private static void accumulateFromSparse(
      int dir, float[][] e, float[][] d, float[][] m, int[] g,
      double[] rMin, double[] rMax)
  {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    int ii=ib;
    for (int il=0; il<nl; ++il)
      d[ii][il] = e[ii][il];
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      int iemax = g[ii];
      int iemin = g[ji];
      int dg = abs(iemax-iemin); // sparse grid delta
      int kmin = (int)ceil( rMin[iemax]*dg);
      int kmax = (int)floor(rMax[iemax]*dg);
      for (int il=0; il<nl; ++il) {
        float dmin = Float.MAX_VALUE;
        int mi = 0;
        for (int k=kmin; k<=kmax; k++) {
          int rk = k*ic;
          int ik = il+rk;
          if (ik<0 || ik>nlm1)
            continue;
          float dc = d[ji][ik];
          if (dc<dmin) {
            dmin = dc;
            mi = rk;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dmin+e[ii][il];
      }
    }
  }

  private static void backtrack(
      int dir, Sampling shifts, float[][] d, float[][] m, float[] u)
  {
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // begining index
    int ie = (dir>0)?nim1:0; // end index
    int is = (dir>0)?1:-1;   // stride
    int ii = ib;
    // Set initial lag for the case that all errors at ii are equal.
    int il = (dir>0)?0:nlm1;
    float dl = d[ii][il]; // Current accumulated error value.

    // Find minimum lag value(dl) and index(il) at trace index ii.
    for (int jl=0; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = (float)shifts.getValue(il);
    while (ii!=ie) {
      ii += is;
      il += (int)m[ii][il];
      u[ii] = (float)shifts.getValue(il);
    }
  }

  private static void print(String s) {
    System.out.println(s);
  }

  /**
   * Class for tracking progress using atomic fields so that progress can be
   * updated from parallel threads.
   */
  private static class ProgressTracker {

    /**
     * Constructor that sets the total number of units to be completed.
     * @param totalWorkUnits
     */
    private ProgressTracker(int totalWorkUnits) {
      _totalWorkUnits = totalWorkUnits;
      _completedWorkUnits = new AtomicInteger(0);
    }

    /**
     * Increments the completed number of work units by one.
     */
    private void worked() {
      _completedWorkUnits.incrementAndGet();
    }

    /**
     * Get the number of completed work units.
     * @return the number of completed work units.
     */
    private int getCompletedWorkUnits() {
      return _completedWorkUnits.get();
    }

    /**
     * Get the percentage complete.
     * @return the percentage complete.
     */
    private float getPercentComplete() {
      return (float)_completedWorkUnits.get()/_totalWorkUnits*100.0f;
    }

    /**
     * Set the value of an AtomicBoolean to indicate whether or not this
     * instance should be canceled.
     * @param cancel
     */
    private void setCanceled(boolean cancel) {
      _canceled.set(cancel);
    }

    /**
     * Get the canceled state of this instance.
     * @return the canceled state of this instance.
     */
    private boolean getCanceled() {
      return _canceled.get();
    }

    private AtomicBoolean _canceled = new AtomicBoolean(false);
    private int _totalWorkUnits;
    private AtomicInteger _completedWorkUnits;

  }

}