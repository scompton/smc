package warp;

import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

public class WarpUtils {

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP trace.
   * @param h the warped PS trace.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(int n1Max, float[] f, float[] h) {
    int n1f = f.length;
    int n1h = h.length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    float[] fs = copy(n1Max,f);
    float[] hs = copy(n1Max,h);
    float scale = 1.0f/n1Max;
    float[] d = sub(hs,fs);
    float frms = sqrt(sum(mul(fs,fs))*scale);
    float hrms = sqrt(sum(mul(hs,hs))*scale);
    float drms = sqrt(sum(mul( d, d))*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][] f, final float[][] h)
  {
    int n1f = f[0].length;
    int n1h = h[0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    int n2 = f.length;
    float scale = 1.0f/(n1Max*n2);
    float[] rms = Parallel.reduce(n2,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i2) {
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i2][i1];
          float hv = h[i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][][] f, final float[][][] h)
  {
    int n1f = f[0][0].length;
    int n1h = h[0][0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n23 = n2*n3;
    float scale = 1.0f/(n1Max*n23);
    float[] rms = Parallel.reduce(n23,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i23) {
        int i2 = i23%n2;
        int i3 = i23/n2;
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i3][i2][i1];
          float hv = h[i3][i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[] vpvs(float[] u) {
    int n = u.length;
    int nm1 = n-1;
    float[] vpvs = new float[n];
    vpvs[ 0 ] = 1.0f + 2.0f*(u[ 1 ]-u[  0  ]); // at i1=0, forward diff
    vpvs[nm1] = 1.0f + 2.0f*(u[nm1]-u[nm1-1]); // at i1=nm1, backward diff
    for (int i1=1; i1<nm1; i1++)
      vpvs[i1] = 1.0f + (u[i1+1]-u[i1-1]);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[][] vpvs(float[][] u) {
    int n2 = u.length;
    float[][] vpvs = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      vpvs[i2] = vpvs(u[i2]);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[][][] vpvs(float[][][] u) {
    int n3 = u.length;
    float[][][] vpvs = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      vpvs[i3] = vpvs(u[i3]);
    return vpvs;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[] extrapolate(int n1, float[] f) {
    int n1f = f.length;
    float v = f[n1f-1];
    float[] ef = new float[n1];
    copy(n1f,f,ef);
    for (int i1=n1f; i1<n1; i1++)
      ef[i1] = v;
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][] extrapolate(int n1, float[][] f) {
    int n2 = f.length;
    float[][] ef = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      ef[i2] = extrapolate(n1,f[i2]);
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][][] extrapolate(int n1, float[][][] f) {
    int n3 = f.length;
    float[][][] ef = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      ef[i3] = extrapolate(n1,f[i3]);
    return ef;
  }

  /**
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin =  Float.MAX_VALUE;
    float emax = -Float.MAX_VALUE;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(
      float[][] e, float ignoreMin, float ignoreMax)
  {
    int nl = e[0].length;
    int n1 = e.length;
    float emin =  Float.MAX_VALUE;
    float emax = -Float.MAX_VALUE;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin && ei>ignoreMin) emin = ei;
        if (ei>emax && ei<ignoreMax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  public static void normalizeErrors(
      float[][][] e, final float ignoreMin, final float ignoreMax)
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin && ei>ignoreMin) emin = ei;
          if (ei>emax && ei<ignoreMax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i3) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            for (int il=0; il<nl; ++il) {
              float ei = ef[i3][i2][i1][il];
              if (ei<emin) emin = ei;
              if (ei>emax) emax = ei;
            }
          }
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  public static void normalizeErrors(
      float[][][][] e, final float ignoreMin, final float ignoreMax)
  {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i3) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            for (int il=0; il<nl; ++il) {
              float ei = ef[i3][i2][i1][il];
              if (ei<emin && ei>ignoreMin) emin = ei;
              if (ei>emax && ei<ignoreMax) emax = ei;
            }
          }
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  public static void normalize(
      float[][] f, final float nmin, final float nmax)
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
	float vi = ff[i2][i1];
	ff[i2][i1] = nrange*(vi-vmin)/range + nmin;
      }
    }});
  }

  public static void normalize(
      float[][][] f, final float nmin, final float nmax)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float vi = ff[i3][i2][i1];
          ff[i3][i2][i1] = nrange*(vi-vmin)/range + nmin;
        }
      }
    }});
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  public static float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n2][n1][nl] of errors.
   * @return transposed array[nl][n2][n1] of errors.
   */
  public static float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }

  public static float[][][] transposeLag12(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n2][nl][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int il=0; il<nl; ++il) {
        for (int i1=0; i1<n1; ++i1) {
          t[i2][il][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }

  public static float[][][] transposeLag23(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n1][n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        t[i1][i2] = e[i2][i1];
      }
    }
    return t;
  }

  public static float[][][][] transposeLag12(float[][][][] e) {
    int nl = e[0][0][0].length;
    int n1 = e[0][0].length;
    int n2 = e[0].length;
    int n3 = e.length;
    float[][][][] t = new float[n3][n2][nl][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          for (int i1=0; i1<n1; ++i1) {
            t[i3][i2][il][i1] = e[i3][i2][i1][il];
          }
        }
      }
    }
    return t;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            ef[i3][i2][i1][il] = (ef[i3][i2][i1][il]-eshift)*escale;
          }
        }
      }
    }});
  }

  private static void print(String s) {
    System.out.println(s);
  }

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

}