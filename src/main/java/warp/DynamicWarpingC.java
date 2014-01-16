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
   * Constructs a 1D dynamic warping.
   * @param uMin minimum shift between PP and PS1.
   * @param uMax maximum shift between PP and PS1.
   * @param s1 first dimension sampling for warping.
   */
  public DynamicWarpingC(double uMin, double uMax, Sampling s1) {
    this(uMin,uMax,s1,null,null);
  }

  /**
   * Constructs a 2D dynamic warping.
   * @param uMin minimum shift between PP and PS1.
   * @param uMax maximum shift between PP and PS1.
   * @param s1 first dimension sampling for warping.
   * @param s2 second dimension sampling.
   */
  public DynamicWarpingC(
      double uMin, double uMax,Sampling s1, Sampling s2)
  {
    this(uMin,uMax,s1,s2,null);
  }

  /**
   * Constructs a 3D dynamic warping.
   * @param uMin minimum shift between PP and PS1.
   * @param uMax maximum shift between PP and PS1.
   * @param s1 first dimension sampling for warping.
   * @param s2 second dimension sampling.
   * @param s3 third dimension sampling.
   */
  public DynamicWarpingC(
      double uMin, double uMax, Sampling s1, Sampling s2, Sampling s3)
  {
    double ds = s1.getDelta(); // shift sampling interval
    int iuMin = (int) ceil(uMin/ds);
    int iuMax = (int)floor(uMax/ds);
    _su = new Sampling(1+iuMax-iuMin,ds,iuMin*ds);
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _r1Min =  0.0;
    _r1Max =  1.0;
    _r2Min = -1.0;
    _r2Max =  1.0;
    _r3Min = -1.0;
    _r3Max =  1.0;
    _c = 1.0f; // Set default compression to do nothing
    _si = new SincInterp();
    //_si.setExtrapolation(Extrapolation.CONSTANT);
    _si.setExtrapolation(Extrapolation.ZERO);
    _ui = new ShiftInterp();
    _ui.setMethod(ShiftInterp.Method.LINEAR);
  }

  /**
   * Constructs a 1D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS,
      float vpvsAvg, double uMin)
  {
    return fromVpVs(sPP,sPS,vpvsAvg,uMin,null,null);
  }

  /**
   * Constructs a 2D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @param s2 second dimension sampling.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS,
      float vpvsAvg, double uMin,
      Sampling s2)
  {
    return fromVpVs(sPP,sPS,vpvsAvg,uMin,s2,null);
  }

  /**
   * Constructs a 3D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @param s2 second dimension sampling.
   * @param s3 third dimension sampling.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS,
      float vpvsAvg, double uMin,
      Sampling s2, Sampling s3)
  {
    double psLast = sPS.getLast();
    double scale = getScale(vpvsAvg);
    double ppLast = psLast*scale;
    double ppl = sPP.valueOfNearest(ppLast);
    ppl = ppl<ppLast ? ppl += sPP.getDelta() : ppl;
    double uMax = psLast-ppl;
    Check.argument(uMax>0,"uMax>0");
    int ippl = sPP.indexOfNearest(ppl);
    Sampling s1 = new Sampling(ippl+1,sPP.getDelta(),sPP.getFirst());
    print("ns1="+(ippl+1)+", nps="+sPS.getCount()+", psLast="+psLast+
      ", ppLast="+ppLast+", ppl="+ppl);
    return new DynamicWarpingC(uMin,uMax,s1,s2,s3);
  }

  /**
   * Constructs a 1D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping. In addition this Factory method computes an initial scale
   * factor for compressing the PS data before warping. This compression
   * constant effectively reduces the number of time shifts the algorithm
   * searches over which can reduce the cost. The {@code vpvsMin} and 
   * {@code vpvsMax} are used to design a fairway of potential time shifts
   * around the constant slope determined from {@code vpvsAvg}.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param vpvsMin
   * @param vpvsMax
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS,
      float vpvsMin, float vpvsMax, float vpvsAvg, float uMin)
  {
    return fromVpVs(sPP,sPS,null,null,vpvsMin,vpvsMax,vpvsAvg,uMin);
  }

  /**
   * Constructs a 2D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping. In addition this Factory method computes an initial scale
   * factor for compressing the PS data before warping. This compression
   * constant effectively reduces the number of time shifts the algorithm
   * searches over which can reduce the cost. The {@code vpvsMin} and 
   * {@code vpvsMax} are used to design a fairway of potential time shifts
   * around the constant slope determined from {@code vpvsAvg}.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param s2 second dimension sampling.
   * @param vpvsMin
   * @param vpvsMax
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS, Sampling s2,
      float vpvsMin, float vpvsMax, float vpvsAvg, float uMin)
  {
    return fromVpVs(sPP,sPS,s2,null,vpvsMin,vpvsMax,vpvsAvg,uMin);
  }

  /**
   * Constructs a 3D dynamic warping from PP/PS samplings and a guess of the
   * average Vp/Vs ratio which is used to compute the maximum PP time used
   * for warping. In addition this Factory method computes an initial scale
   * factor for compressing the PS data before warping. This compression
   * constant effectively reduces the number of time shifts the algorithm
   * searches over which can reduce the cost. The {@code vpvsMin} and 
   * {@code vpvsMax} are used to design a fairway of potential time shifts
   * around the constant slope determined from {@code vpvsAvg}.
   * @param sPP the PP sampling.
   * @param sPS the PS sampling.
   * @param s2 second dimension sampling.
   * @param s3 third dimension sampling.
   * @param vpvsMin
   * @param vpvsMax
   * @param vpvsAvg
   * @param uMin minimum shift between PP and PS. Theoretically, this should
   *   be zero, however this cannot be assumed to be true after processing.
   * @return a DynamicWarpingC instance.
   */
  public static DynamicWarpingC fromVpVs(
      Sampling sPP, Sampling sPS, Sampling s2, Sampling s3,
      float vpvsMin, float vpvsMax, float vpvsAvg, float uMin)
  {
    double psLast = sPS.getLast();
    double scale = getScale(vpvsAvg);
    double ppLast = psLast*scale;
    double ppl = sPP.valueOfNearest(ppLast);
    ppl = ppl<ppLast ? ppl += sPP.getDelta() : ppl;
    double uMax = psLast-ppl;
    Check.argument(uMax>0,"uMax>0");
    float r = (float)((uMax-uMin)/(ppl-sPP.getFirst()));
    float rMin = (vpvsMin-1.0f)*0.5f;
    float rMax = (vpvsMax-1.0f)*0.5f;
    float uMinNew = (float)(ppl*rMin+uMin-uMax);
    float uMaxNew = (float)(ppl*rMax+uMin-uMax);
    //print("r="+r+", rMin="+rMin+", rMax="+rMax);
    float c = 1.0f+r; // compression constant
    int ippl = sPP.indexOfNearest(ppl);
    Sampling s1 = new Sampling(ippl+1,sPP.getDelta(),sPP.getFirst());
    print("DynamicWarpingC fromVpVs:\n"+
          "  PP max time:      "+sPP.getLast()+"\n"+
          "  PS max time:      "+sPS.getLast()+"\n"+
          "  Warping max time: "+s1.getLast()+"\n"+
          "  Minimum shift:    "+uMin+"\n"+
          "  Maximum shift:    "+uMax+"\n"+
          "  Scale factor c:   "+c+"\n"+
          "  Minimum shift c:  "+uMinNew+"\n"+
          "  Maximum shift c:  "+uMaxNew+"\n");
    DynamicWarpingC dw = new DynamicWarpingC(uMinNew,uMaxNew,s1,s2,s3);
    dw.setCompression(c);
    return dw;
  }

  /**
   * Set the constant scalar used for an initial compression of the PS data.
   * If not set, the default value is 1.0f, which has no effect.
   * @param c the compreesion constant.
   */
  public void setCompression(float c) {
    _c = c;
  }

  /**
   * Get the constant scalar used for an initial compression of the PS data.
   * @return the compreesion constant.
   */
  public float getCompression() {
    return _c;
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

  /**
   * Sets bounds on strains for this dynamic warping. These bounds can vary 
   * from sample to sample, but are ultimately only enforced where the bounds
   * coincide with coarse samples.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   */
  public void setStrainLimits(double[] r1Min, double[] r1Max) {
    setStrainLimits(r1Min,r1Max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping. These bounds can vary 
   * from sample to sample in the first dimension, 1, but are ultimately only
   * enforced where the bounds coincide with coarse samples.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
      double[] r1Min, double[] r1Max,
      double r2Min, double r2Max)
  {
    setStrainLimits(r1Min,r1Max,r2Min,r2Max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping. These bounds can vary 
   * from sample to sample in the first dimension, 1, but are ultimately only
   * enforced where the bounds coincide with coarse samples.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   * @param r3Min lower bound on strain in 3rd dimension.
   * @param r3Max upper bound on strain in 3rd dimension.
   */
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
   * Set the method used for interpolating shifts. LINEAR interpolation is most
   * consistent with the warping algorithm used. If smoother shifts are desired
   * then MONOTONIC can be used. SPLINE is also a valid option, but is not
   * recommended since returned time shifts cannot be guaranteed to be
   * monotonically increasing.
   */
  public void setInterpolationMethod(ShiftInterp.Method method) {
    _ui.setMethod(method);
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
   * Get the shift sampling.
   * @return the shift sampling.
   */
  public Sampling getSamplingU() {
    return _su;
  }

  /**
   * Get the sampling for warping in the first (fast) dimension. This can be
   * different than the sampling of the input data.
   * @return the sampling for warping in the first dimension.
   */
  public Sampling getSampling1() {
    return _s1;
  }

  /**
   * Get the sampling in the second (slowest dimension in 2D, middle in 3D)
   * dimension.
   * @return the sampling in the second dimension.
   */
  public Sampling getSampling2() {
    return _s2;
  }

  /**
   * Get the sampling in the third (slowest) dimension.
   * @return the sampling in the third dimension.
   */
  public Sampling getSampling3() {
    return _s3;
  }

  /**
   * Find shifts for 1D traces. For the input trace {@code f[n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng1]}. This length matches the length of array {@code g1}
   * and the coordinates of the subsampled shifts are specified by contents of
   * this array.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is
   * {@code ui[n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP trace.
   * @param sg the PS trace sampling.
   * @param g the PS trace.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  coordinates.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(
      Sampling sf, float[] f,
      Sampling sg, float[] g,
      float[] g1)
  {
    return findShifts(sf,f,sg,g,g1,null);
  }

  /**
   * Find shifts for 1D traces. For the input trace {@code f[n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng1]}. This length matches the length of array {@code g1}
   * and the coordinates of the subsampled shifts are specified by contents of
   * this array. Optionally, known shifts can be specified as a shift,sample
   * coordinate pair in the {@code u1} array.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is
   * {@code ui[n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP trace.
   * @param sg the PS trace sampling.
   * @param g the PS trace.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  coordinates.
   * @param u1 array of shift,sample coordinate pairs that set hard constraints
   *   on the error function.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(
      Sampling sf, float[] f,
      Sampling sg, float[] g,
      float g1[], float[] u1)
  {
    int ne = _s1.getCount();
    double[] r1Min = getR1Min(ne);
    double[] r1Max = getR1Max(ne);
    g1 = KnownShiftUtil.getG1(_s1,g1,u1,r1Min,r1Max);
    int[] g1i = gridCoordsToSamples(_s1,g1);
    printInfo(g1i.length,1,1);
    float[][] e = computeErrors(sf,f,sg,g);
    fixShifts(e,u1); // NOP if u1 is null
    float[][] es = smoothErrors(e,r1Min,r1Max,g1i);
    float[][][] dm = accumulateForward(es,r1Min,r1Max,g1i);
    float[] u = backtrackReverse(dm[0],dm[1]);
    checkSlopes(u,g1i,r1Min,r1Max);
    return _ui.interpolate(sf,g1,u);
  }

  /**
   * Find shifts for previously computed alignment errors {@code e}. Shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng1]}. This length matches the length of array {@code g1}
   * and the coordinates of the subsampled shifts are specified by contents of
   * this array.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the length {@code n1=sf.getCount()}, that 
   * is {@code ui[n1]}.
   * @param sf the PP trace sampling.
   * @param e previously computed alignment errors.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  coordinates.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(Sampling sf, float[][] e, float[] g1) {
    int ne = _s1.getCount();
    Check.argument(ne==e.length,"ne==e.length");
    int[] g1i = gridCoordsToSamples(_s1,g1);
    printInfo(g1i.length,1,1);
    double[] r1Min = getR1Min(ne);
    double[] r1Max = getR1Max(ne);
    float[][] es = smoothErrors(e,r1Min,r1Max,g1i);
    float[][][] dm = accumulateForward(es,r1Min,r1Max,g1i);
    float[] u = backtrackReverse(dm[0],dm[1]);
    checkSlopes(u,g1i,r1Min,r1Max);
    return _ui.interpolate(sf,g1,u);
  }

  /**
   * Find shifts for 2D images. For the input image {@code f[n2][n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays
   * {@code g1,g2} and the coordinates of the subsampled shifts are specified 
   * by contents of these arrays.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n2][n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP traces.
   * @param sg the PS trace sampling.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @return shifts for 2D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][] findShifts(
      Sampling sf, float[][] f,
      Sampling sg, float[][] g,
      float[][] g1, float[] g2) throws CancellationException
  {
    return findShifts(sf,f,sg,g,g1,g2,new HashMap<Integer,float[]>());
  }

  /**
   * Find shifts for 2D images. For the input image {@code f[n2][n1]}, shifts
   * {@code u} are computed on a subsampled grid such that the computed shifts
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays
   * {@code g1,g2} and the coordinates of the subsampled shifts are specified 
   * by contents of these arrays.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n2][n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP traces.
   * @param sg the PS trace sampling.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param u1Map
   * @return shifts for 2D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][] findShifts(
      final Sampling sf, final float[][] f,
      final Sampling sg, final float[][] g,
      final float[][] g1, final float[] g2,
      final Map<Integer,float[]> u1Map)
      throws CancellationException
  {
    final Sampling su = _su;
    final Sampling se = _s1;
    final Sampling s2 = sampling2(f);
    final int nu = su.getCount();
    final int ne = se.getCount();
    final int n2 = s2.getCount();
    final int nf = sf.getCount();
    final int ng = sg.getCount();
    final double de = se.getDelta();
    final double df = sf.getDelta();
    final double dg = sg.getDelta();
    final double fe = se.getFirst();
    final double ff = sf.getFirst();
    final double fg = sg.getFirst();
    final double[] uv = su.getValues();
    final int ng2 = g2.length;
    final int ng1 = g1[0].length;
    Check.argument(n2==g1.length,"n2==g1.length");
    printInfo(ng1,ng2,1);

    // Round to integer indices, use float values for interpolation only.
    final int[][] g1i = gridCoordsToSamples(_s1,g1);
    final int[] g2i = gridCoordsToSamples(s2,g2);
    final double[] r1Min = getR1Min(ne);
    final double[] r1Max = getR1Max(ne);
    final double[] r2Min = getR2Min(n2);
    final double[] r2Max = getR2Max(n2);

    /*
     * Initialize ProgressTracker:
     *     smooth1=n2 units,
     *     smooth2=ng1 units,
     *     find shifts=ng2 units,
     *     interpolation=1 unit.
     */
    int totalWorkUnits = n2+ng1+ng2+1;
    final ProgressTracker pt = new ProgressTracker(totalWorkUnits);
    startProgressThread(pt);

    // Smooth1
    //Stopwatch s = new Stopwatch();
    //s.start();
    //print("Smoothing 1st dimension...");
    final float[][][] es1 = new float[n2][ng1][nu];
    final Parallel.Unsafe<float[][]> eu = new Parallel.Unsafe<>();
    final Parallel.Unsafe<float[]> fiu = new Parallel.Unsafe<>();
    final Parallel.Unsafe<float[]> giu = new Parallel.Unsafe<>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      if (pt.getCanceled()) {
        Thread.currentThread().interrupt();
        throw new CancellationException();
      }
      float[][] e = eu.get();
      float[] fi = fiu.get();
      float[] gi = giu.get();
      if (e==null) eu.set(e=new float[ne][nu]);
      if (fi==null) fiu.set(fi=new float[ne]);
      if (gi==null) giu.set(gi=new float[ne]);
      computeErrors(nf,df,ff,f[i2],fi,
                    ng,dg,fg,g[i2],gi,
                    ne,de,fe,e,nu,uv);
      // float[] u1 = u1Map.get(i2);
      // fixShifts(e,u1);
      // int[] g1ks = KnownShiftUtil.getG1(g1[i2],l1,r1min,r1max);
      // es1[i2] = smoothErrors(e,r1Min,r1Max,g1ks);
      es1[i2] = smoothErrors(e,r1Min,r1Max,g1i[i2]);
      pt.worked();
    }});
    //print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    // Smooth2
    //s.restart();
    //print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,r2Min,r2Max,g2i,pt);
    //print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);

    // Find shifts on coarse grid
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],r1Min,r1Max,g1i[g2i[i2]]);
      u[i2] = backtrackReverse(dm[0],dm[1]);
      pt.worked();
    }});
    checkSlopes(u,g1i,g2i,r1Min,r1Max);

    // Interpolate shifts to original grid.
    float[][] ui = _ui.interpolate(sf,s2,g1,g2,u);
    pt.worked();
    assert pt.getPercentComplete()==100.0f;
    return ui;
  }

  /**
   * Find shifts for 3D images. For the input image {@code f[n3][n2][n1]},
   * shifts {@code u} are computed on a subsampled grid such that the computed
   * shifts are {@code u[ng3][ng2][ng1]}. These lengths match the length of
   * arrays {@code g1,g2,g3} and the coordinates of the subsampled shifts are
   * specified by contents of these arrays.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n3][n2][n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP traces.
   * @param sg the PS trace sampling.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param g3 array of size [ng3] specifying third dimension sparse grid
   *  locations.
   * @return shifts for 3D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][][] findShifts(
      Sampling sf, float[][][] f,
      Sampling sg, float[][][] g,
      float[][][] g1, float[] g2, float[] g3) throws CancellationException
  {
    return findShifts(sf,f,sg,g,g1,g2,g3,new HashMap<Integer,float[]>());
  }

  /**
   * Find shifts for 3D images. For the input image {@code f[n3][n2][n1]},
   * shifts {@code u} are computed on a subsampled grid such that the computed
   * shifts are {@code u[ng3][ng2][ng1]}. These lengths match the length of
   * arrays {@code g1,g2,g3} and the coordinates of the subsampled shifts are
   * specified by contents of these arrays.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is
   * {@code ui[n3][n2][n1]}.
   * @param sf the PP trace sampling.
   * @param f the PP traces.
   * @param sg the PS trace sampling.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid
   *  locations.
   * @param g3 array of size [ng3] specifying third dimension sparse grid
   *  locations.
   * @param u1Map
   * @return shifts for 3D images.
   * @throws CancellationException if a {@link WarperWorkTracker} instance is
   *  set with the {@link #setWorkTracker(WarperWorkTracker)} method, and this
   *  the {@link WarperWorkTracker#isCanceled()} method returns {@code true}.
   */
  public float[][][] findShifts(
      final Sampling sf, final float[][][] f,
      final Sampling sg, final float[][][] g,
      final float[][][] g1, final float[] g2, final float[] g3,
      final Map<Integer,float[]> u1Map)
      throws CancellationException
  {
    final Sampling su = _su;
    final Sampling se = _s1;
    final Sampling s2 = sampling2(f);
    final Sampling s3 = sampling3(f);
    final int nu = su.getCount();
    final int ne = se.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    final int nf = sf.getCount();
    final int ng = sg.getCount();
    final double de = se.getDelta();
    final double df = sf.getDelta();
    final double dg = sg.getDelta();
    final double fe = se.getFirst();
    final double ff = sf.getFirst();
    final double fg = sg.getFirst();
    final double[] uv = su.getValues();
    final int ng1 = g1[0][0].length;
    final int ng2 = g2.length;
    final int ng3 = g3.length;
    Check.argument(n2==g1[0].length,"n2==g1[0].length");
    Check.argument(n3==g1.length,"n2==g1.length");
    printInfo(ng1,ng2,ng3);

    // Round to integer indices, use float values for interpolation only.
    final int[][][] g1i = gridCoordsToSamples(_s1,g1);
    final int[] g2i = gridCoordsToSamples(s2,g2);
    final int[] g3i = gridCoordsToSamples(s3,g3);
    final double[] r1Min = getR1Min(ne);
    final double[] r1Max = getR1Max(ne);
    final double[] r2Min = getR2Min(n2);
    final double[] r2Max = getR2Max(n2);
    final double[] r3Min = getR3Min(n3);
    final double[] r3Max = getR3Max(n3);

    /*
     * Initialize ProgressTracker:
     *     smooth1=n2*n3 units,
     *     smooth2=n3*ng1 units,
     *     smooth3=ng1*ng2 units,
     *     find shifts=ng2*ng3 units,
     *     interpolation=1 unit.
     */
    int totalWorkUnits = (n2*n3)+(n3*ng1)+(ng1*ng2)+(ng2*ng3)+1;
    final ProgressTracker pt = new ProgressTracker(totalWorkUnits);
    startProgressThread(pt);

    // Smooth1
    //Stopwatch s = new Stopwatch();
    //s.start();
    //print("Smoothing 1st dimension...");
    final float[][][][] es1 = new float[n3][n2][ng1][nu];
    final Parallel.Unsafe<float[][]> eu = new Parallel.Unsafe<>();
    final Parallel.Unsafe<float[]> fiu = new Parallel.Unsafe<>();
    final Parallel.Unsafe<float[]> giu = new Parallel.Unsafe<>();
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
      float[] fi = fiu.get();
      float[] gi = giu.get();
      if (e==null) eu.set(e=new float[ne][nu]);
      if (fi==null) fiu.set(fi=new float[ne]);
      if (gi==null) giu.set(gi=new float[ne]);
      computeErrors(nf,df,ff,f[i3][i2],fi,
                    ng,dg,fg,g[i3][i2],gi,
                    ne,de,fe,e,nu,uv);
      // int[] l1 = ilMap.get(i23);
      // fixShifts(e,l1);
      // int[] g1ks = KnownShiftUtil.getG1(g1[i3][i2],l1,r1min,r1max);
      // es1[i3][i2] = smoothErrors(e,r1Min,r1Max,g1ks);
      es1[i3][i2] = smoothErrors(e,r1Min,r1Max,g1i[i3][i2]);
      pt.worked();
    }});
    //print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    // Smooth 2
    //s.restart();
    //print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,r2Min,r2Max,g2i,pt);
    //print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);

    // Smooth 3
    //s.restart();
    //print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,r3Min,r3Max,g3i,pt);
    normalizeErrors(es);
    //print("Finished 3rd dimension smoothing in "+s.time()+" seconds");

    // Find shifts on coarse grid
    final float[][][] u = new float[ng3][ng2][];
    int ng23 = ng3*ng2;
    Parallel.loop(ng23,new Parallel.LoopInt() {
    public void compute(int i23) {
      int i2 = i23%ng2;
      int i3 = i23/ng2;
      float[][][] dm = accumulateForward(
          es[i3][i2],r1Min,r1Max,g1i[g3i[i3]][g2i[i2]]);
      u[i3][i2] = backtrackReverse(dm[0],dm[1]);
      pt.worked();
    }});
    checkSlopes(u,g1i,g2i,g3i,r1Min,r1Max);

    // Interpolate shifts to original grid.
    float[][][] ui = _ui.interpolate(sf,s2,s3,g1,g2,g3,u);
    pt.worked();
    assert pt.getPercentComplete()==100.0f;
    return ui;
  }

  /**
   * Computes alignment errors for traces {@code f} and {@code g}.
   * @return 2D array of alignment errors.
   */
  public float[][] computeErrors(
      Sampling sf, float[] f,
      Sampling sg, float[] g)
  {
    Sampling su = _su;
    Sampling se = _s1;
    int nu = su.getCount();
    int ne = se.getCount();
    int nf = sf.getCount();
    int ng = sg.getCount();
    double de = se.getDelta();
    double df = sf.getDelta();
    double dg = sg.getDelta();
    double fe = se.getFirst();
    double ff = sf.getFirst();
    double fg = sg.getFirst();
    float[][] e = new float[ne][nu];
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    computeErrors(
      nf,df,ff,f,fi,
      ng,dg,fg,g,gi,
      ne,de,fe,e,
      nu,su.getValues());
    return e;
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

  ///////////////////////////////////////////////////////////////////////////
  // These methods should typically only be called for research applications

  /**
   * Compute alignment errors for 2D traces.
   * @return alignment errors for 2D traces.
   */
  public float[][][] computeErrors2(
      final Sampling sf, final float[][] f,
      final Sampling sg, final float[][] g)
  {
    int n2 = f.length;
    final float[][][] e = new float[n2][][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      e[i2] = computeErrors(sf,f[i2],sg,g[i2]);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 2D images.
   * @return summed alignment errors for 2D images.
   */
  public float[][] computeErrorsSum(
      final Sampling sf, final float[][] f,
      final Sampling sg, final float[][] g)
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    float[][] e = Parallel.reduce(n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      return computeErrors(sf,f[i2],sg,g[i2]);
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
  public float[][] computeErrorsSum(
      final Sampling sf, final float[][][] f,
      final Sampling sg, final float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      return computeErrors(sf,f[i3][i2],sg,g[i3][i2]);
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  public float[][][] getSmoothErrors(float[][] e, int[] g, int n1) {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    double[] r1Min = getR1Min(_s1.getCount());
    double[] r1Max = getR1Max(_s1.getCount());
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
    // backtrack(-1,_sl1,d,m,u);
    backtrack(-1,_su,d,m,u);
    return u;
  }

  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    // backtrack(1,_sl1,d,m,u);
    backtrack(1,_su,d,m,u);
    return u;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private
  private Sampling _su; // Shift sampling
  private Sampling _s1, _s2, _s3;
  private double _r1Min, _r1Max;
  private double _r2Min, _r2Max;
  private double _r3Min, _r3Max;
  private double[] _r1MinA, _r1MaxA;
  private float _c; // compression factor before warping
  private int _nl; // number of lags
  private Sampling _sl1; // sampling of pp to ps1 lags
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
  private ShiftInterp _ui;
  private WarperWorkTracker _wwt;
  private Almost a = new Almost();

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

  private void printInfo(int ng1, int ng2, int ng3) {
    int n1 = _s1.getCount();
    int n2 = _s2==null ? 1 : _s2.getCount();
    int n3 = _s3==null ? 1 : _s3.getCount();
    int nu = _su.getCount();
    float e1Mem = (float)ng1*n2*n3*nu*4.0f*BYTES_TO_MB;
    float e2Mem = (float)ng1*ng2*n3*nu*4.0f*BYTES_TO_MB;
    float e3Mem = (float)ng1*ng2*ng3*nu*4.0f*BYTES_TO_MB;
    print("findShifts info:");
    print("  Input data samples (n1,n2,n3): ("+n1+","+n2+","+n3+")");
    print("  Coarse grid samples: (ng1,ng2,ng3): ("+ng1+","+ng2+","+ng3+")");
    print("  Number of lags: "+nu);
    print("  Alignment error smooth 1 memory: "+e1Mem+" MB");
    print("  Alignment error smooth 2 memory: "+((n2>1)?(e2Mem+" MB"):"NA"));
    print("  Alignment error smooth 3 memory: "+((n3>1)?(e3Mem+" MB"):"NA"));
  }

  private Sampling sampling2(float[][] f) {
    if (_s2!=null) {
      Check.argument(f.length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
  }

  private Sampling sampling2(float[][][] f) {
    if (_s2!=null) {
      Check.argument(f[0].length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f[0].length,1.0,0.0);
    }
  }

  private Sampling sampling3(float[][][] f) {
    if (_s3!=null) {
      Check.argument(f.length==_s3.getCount(),"valid sampling3");
      return _s3;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
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

  private void checkSlopes(
      float[][][] u, int[][][] g1, int[] g2, int[] g3,
      double[] r1Min, double[] r1Max)
  {
    int n3 = u.length;
    for (int i3=0; i3<n3; i3++) {
      print("checking i3 "+i3);
      checkSlopes(u[i3],g1[g3[i3]],g2,r1Min,r1Max);
    }
  }

  private void checkSlopes(
      float[][] u, int[][] g1, int[] g2, double[] r1Min, double[] r1Max)
  {
    int n2 = u.length;
    for (int i2=0; i2<n2; i2++)
      checkSlopes(u[i2],g1[g2[i2]],r1Min,r1Max);
  }

  private void checkSlopes(
      float[] u, int[] g1, double[] r1Min, double[] r1Max)
  {
    int ng = g1.length;
    float uLast = (float)_su.getLast();
    float dui = 1.0f/(float)_su.getDelta();
    for (int ig=1; ig<ng; ig++) {
      int igm1 = ig-1;
      int ii = g1[ig  ]; // current index in finely sampled array
      int ji = g1[igm1]; // previous index in finely sampled array
      float n = (int)((u[ig]-u[igm1])*dui+0.5f); // numerator
      float d = ii-ji; // denominator
      float r = n/d; // slope
      assert (a.ge(r,r1Min[ii]) && a.le(r,r1Max[ii])) || u[ig]==uLast :
        "Slope constraints violated:\n"+
        "  n="+n+", d="+d+", r="+r+",\n"+
        "  rmin="+r1Min[ii]+", rmax="+r1Max[ii]+",\n"+
        "  u["+ig+"]="+u[ig]+", u["+igm1+"]="+u[igm1]+", uLast="+uLast+"\n"+
        "  g1["+ig+"]="+ii+", g1["+ig+"]="+ji;
    }
  }

  private void computeErrors(
      int nf, double df, double ff, float[] f, float[] fi,
      int ng, double dg, double fg, float[] g, float[] gi,
      int ne, double de, double fe, float[][] e,
      int nu, double[] uv)
  {
    double lg = (ng-1)*dg;
    _si.interpolate(nf,df,ff,f,ne,de,fe,fi);
    for (int iu=0; iu<nu; iu++) {
      double u = uv[iu];
      _si.interpolate(ng,dg,fg,g,ne,de,fe+u,gi);
      for (int ie=0; ie<ne; ie++) {
        double gv = fe+ie*de+u;
        e[ie][iu] = (gv<fg || gv>lg) ? Float.NaN : error(fi[ie],gi[ie]);
      }
    }
    reflectErrors(e);
  }

  /**
   * Replaces NaNs from out-of-bounds errors with reflected error values.
   * @param e alignment errors.
   */
  private void reflectErrors(float[][] e) {
    reflect(e, 1);
    reflect(e,-1);
  }

  /**
   * Replaces NaNs from out-of-bounds errors with reflected error values.
   * This method checks out-of-bounds errors from one direction given by the
   * {@code dir} argument. In the forward direction ({@code dir>0}) NaNs
   * are removed in the direction of increasing shifts. In the reverse
   * direction ({@code dir<0}) NaNs are removed in the direction of decreasing
   * shifts. The method terminates when the first or last value is not a NaN.
   * @param e alignment errors.
   * @param dir {@code >0} forward direction, {@code <0} reverse direction.
   */
  private void reflect(float[][] e, int dir) {
    int nu = e[0].length;
    int ne = e.length;
    int ieb = dir>0 ?  0 : ne-1; // beginning error index
    int iub = dir>0 ?  0 : nu-1; // beginning shift index
    int iee = dir>0 ? ne :   -1; // end error index
    int iue = dir>0 ? nu :   -1; // end shift index
    int ire = dir>0 ? -1 :   nu; // end shift reflected index
    int s   = dir>0 ?  1 :   -1; // stride
    for (int ie=ieb; ie!=iee; ie+=s) {
      if (Float.isNaN(e[ie][iub])) {
        int iu = iub+s;
        while (iu!=iue && Float.isNaN(e[ie][iu])) iu += s;
        int ir = iu-s; // reflected index, is 1 sample back/forward.
        for (iu+=s; iu!=iue && ir!=ire; iu+=s, ir-=s)
          e[ie][ir] = e[ie][iu];
        for (; ir!=ire; ir-=s) // extrapolate where we couldn't reflect
          e[ie][ir] = e[ie][ir+s];
      } else {
        break;
      }
    }
  }

  private void fixShifts(float[][] e, float[] u1) {
    if (u1==null) return;
    int n1 = e.length;
    int nu = e[0].length;
    for (int i=0; i<u1.length; i+=2) {
      float xu = u1[i  ];
      float x1 = u1[i+1];
      int lag = _su.indexOfNearest(xu);
      int i1  = _s1.indexOfNearest(x1);
      if (i1<0 || i1>=n1 || lag<0 || lag>= nu)
        continue;
      for (int il=0; il<nu; il++)
        e[i1][il] = Float.MAX_VALUE;
      e[i1][lag] = 0.0f;
    }
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
          int wu = pt.getCompletedWorkUnits();
          _wwt.worked(wu-wus);
        }
      });
      thread.start();
    }
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
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
      // kmin = kmin>kmax ? kmax : kmin;
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
        } else { // interpolate
          linearInterp(ie,je,is,r,ils,ile,dm,e);
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

  private static void linearInterp(
      int ie, int je, int is, float r, int ils, int ile,
      float[] dm, float[][] e) 
  {
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
