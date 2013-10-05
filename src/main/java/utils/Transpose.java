package utils;

import edu.mines.jtk.util.*;

public final class Transpose {

  /**
   * Transposes array {@code f} from 21 to 12 (order is slow to fast).
   * @param f array[n2][n1].
   * @return transposed array[n1][n2].
   */
  public static float[][] transpose12(final float[][] f) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] t = new float[n1][n2];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; i1++)
        t[i1][i2] = f[i2][i1];
    }});
    return t;
  }

  /**
   * Transposes array {@code f} from 321 to 312 (order is slow to fast).
   * @param f array[n3][n2][n1].
   * @return transposed array[n3][n1][n2].
   */
  public static float[][][] transpose312(final float[][][] f) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] t = new float[n3][n1][n2];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++)
        for (int i1=0; i1<n1; i1++)
          t[i3][i1][i2] = f[i3][i2][i1];
    }});
    return t;
  }

  /**
   * Transposes array {@code f} from 321 to 132 (order is slow to fast).
   * @param f array[n3][n2][n1].
   * @return transposed array[n1][n3][n2].
   */
  public static float[][][] transpose132(final float[][][] f) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] t = new float[n1][n3][n2];
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
          t[i1][i3][i2] = f[i3][i2][i1];
    }});
    return t;
  }

  /**
   * Transposes array {@code f} from 321 to 231 (order is slow to fast).
   * @param f array[n3][n2][n1].
   * @return transposed array[n2][n3][n1].
   */
  public static float[][][] transpose231(final float[][][] f) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] t = new float[n2][n3][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i3=0; i3<n3; ++i3)
        t[i2][i3] = f[i3][i2];
    }});
    return t;
  }

  /**
   * Transposes array {@code f} from 4321 to 4312 (order is slow to fast).
   * @param f array[n4][n3][n2][n1].
   * @return transposed array[n4][n3][n1][n2].
   */
  public static float[][][][] transpose4312(final float[][][][] f) {
    final int n1 = f[0][0][0].length;
    final int n2 = f[0][0].length;
    final int n3 = f[0].length;
    final int n4 = f.length;
    final float[][][][] t = new float[n4][n3][n1][n2];
    Parallel.loop(n4,new Parallel.LoopInt() {
    public void compute(int i4) {
      for (int i3=0; i3<n3; i3++)
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            t[i4][i3][i1][i2] = f[i4][i3][i2][i1];
    }});
    return t;
  }

}