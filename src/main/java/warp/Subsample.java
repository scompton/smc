package warp;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Subsample {

  /**
   * Subsamples indices of an array of length n, ensuring that
   * the subsampling includes the first and last index. To meet
   * this requirement, the interval of the sparse indices may 
   * be greater than or equal to {@code d}, 
   * @param n the length of the array to subsample
   * @param d the minimum subsampled interval
   * @return an integer array of sparse grid indices.
   */ 
  public static int[] subsample(int n, int d) {
    if (d>=n)
      d = max(n-1,1);
    int m = 1+(n-1)/d;
    double dd = (double)(n-1)/(double)(m-1);
    int[] g = new int[m];
    for (int ig=0; ig<m; ++ig)
      g[ig] = (int)(ig*dd+0.5);
    return g;
  }
  
  /**
   * Subsamples indices of an array of floats preferentially
   * selecting the indices of f with the highest absolute value
   * amplitudes. This selection ensures that the first and 
   * last sample are always included in the subsampled indices,
   * and that the interval between selected indices is greater
   * than or equal to d.
   * @param f the input array to preferentially subsample.
   * @param d the constraint on the subsample interval.
   * @return an integer array of sparse grid indices. The 
   *  length of this array is unknown. It will include as 
   *  many indices as possible that satisfy the interval
   *  constraint d.
   */
  public static int[] subsample(float[] f, int d) {
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(f);
    quickIndexSort(fa,i);
//    StringBuffer sb = new StringBuffer();
//    for (int ii=0; ii<10; ii++) {
//      sb.append(i[nm-ii]+",");
//    }
//    System.out.println("Top 10: "+sb.toString());
    
    // Get maximum amplitute indices as long as they
    // are greater than or equal to d.
    List<Integer> gl = new ArrayList<Integer>();
    gl.add(0); gl.add(nm);
    for (int ii=0; ii<n; ii++) {
      int im = i[nm-ii]; // next highest amplitude index
      int s = gl.size();
      boolean add = true;
      for (int ig=0; ig<s; ig++) {
        if (abs(im-gl.get(ig))<d) {
          add = false;
          break;
        }
      }
      if (add)
        gl.add(im);
    }
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }

  /**
   * Subsamples indices of an array of floats preferentially
   * selecting the indices of f with the highest absolute value
   * amplitudes. This selection ensures that the first and 
   * last sample are always included in the subsampled indices,
   * and that the number of subsamples is equal to {@code ng}. The 
   * interval between selected indices decreased from {@code d } as
   * needed to meet this requirement.
   * @param f the input array to preferentially subsample.
   * @param d the constraint on the subsample interval.
   * @param ng the number of samples included in the subsampled array.
   * @return an integer array of sparse grid indices. The 
   *  length of this array is {@code ng} and the interval {@code d} 
   *  may be modified to meet this requirement. 
   */
  public static int[] subsample(float[] f, int d, int ng) {
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(f);
    quickIndexSort(fa,i);
//    StringBuffer sb = new StringBuffer();
//    for (int ii=0; ii<10; ii++) {
//      sb.append(i[nm-ii]+",");
//    }
//    System.out.println("Top 10: "+sb.toString());
    
    // Get maximum amplitude indices, adjusting the interval
    // d so that the subsampled array contains ng samples.
    List<Integer> gl = new ArrayList<Integer>();
    for (; d>=1; d--) {
      gl.clear();
      gl.add(0); gl.add(nm);
      for (int ii=0; ii<n; ii++) {
        int im = i[nm-ii]; // next highest amplitude index
        int s = gl.size();
        if (s==ng) // size matches desired sparse array size, quit
          break;
        boolean add = true;      
        for (int ig=0; ig<s; ig++) {
          if (abs(im-gl.get(ig))<d) {
            add = false;
            break;
          }
        }
        if (add)
          gl.add(im);
      }
      if (gl.size()==ng) { // size matches desired sparse array size, quit
        System.out.println("final d="+d);
        break;
      }
    }
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }

  public static void main(String[] args) {
    if (args.length != 2) {
      System.out.println("Usage: java Decompose numberOfSamples delta");
      System.exit(0);
    }
    int n = Integer.parseInt(args[0]);
    int d = Integer.parseInt(args[1]);
    int[] g = subsample(n,d); 
    dump(g);
  }

}
