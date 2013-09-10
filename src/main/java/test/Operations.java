package test;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.Arrays;
import java.util.concurrent.CancellationException;

import edu.mines.jtk.util.Parallel;

public class Operations {

  public static void CGsolver(float[] m, float[] r, float[] dm,
      float[] dr, float[] s, float[] S){
      int n12 = m.length;
      int nr = r.length;
      float a = -dot(r,dr)/dot(dr,dr);
      for (int i = 0; i < n12; i++){
        s[i] = dm[i]*a;
        S[i] = dr[i]*a;
        m[i] += s[i]; 
        r[i] += S[i];
      }
      for (int i = n12; i < nr; i++){
        S[i] = dr[i]*a;
        r[i] += S[i];
      }
  }
  
  public static void CGsolver(final float[] m, final float[] r, final float[] dm,
      final float[] dr, final float[] s, final float[] S, float GG, float SS, float GS, float GR, float SR){
    int n12 = m.length;
    int nr = r.length;
    for (int i =0; i < nr; i++){
      GG += dr[i]*dr[i];
      SS +=  S[i]*S[i];
      GS += dr[i]*S[i];
      GR += dr[i]*r[i];
      SR +=  S[i]*r[i];
    }
      float det = - GG*SS*max(1.0f-(GS/GG)*(GS/SS),1E-12f);
      final float a = ( SS*GR - GS*SR)/det;
      final float b = (-GS*GR + GG*SR)/det;
      Parallel.loop(n12,new Parallel.LoopInt() {
        public void compute(int i) {
          m[i] += (a*dm[i]+b*s[i]); //m[i] += s[i];
          r[i] += (a*dr[i]+b*S[i]); //r[i] += S[i];
        }
      });
      for (int i = n12; i < nr; i++)
        r[i] += (a*dr[i]+b*S[i]);
  }

  
  public static void convert(float[][] im2D, int _n1, int _n2, float[] im){
    // from 2D to 1D array
//    if(im == null) 
//      im = new float[_n1*_n2];
    for (int i=0; i < _n2; i++)
      System.arraycopy(im2D[i], 0, im, i*_n1, _n1);
  }
  
  
  public static void convert(float[] im, int _n1, int _n2, float[][] model2D){
    // from 1D to 2D array
    for (int i=0; i < _n2; i++)
      System.arraycopy(im, i*_n1, model2D[i], 0, _n1);
  }
  
  public static float dot(float[] a, float[] b) throws CancellationException{
    int na = a.length;
    int nb = b.length;
    if(na!=nb) throw new CancellationException( "Dimensions don't match" ); 
    float ab = 0;
    for(int i = 0; i < na; i++)
      ab += a[i]*b[i];
    return ab;
  }
  
  public static float l2Norm(float[] r){
    int n = r.length;
    float rr = 0;
    for(int i = 0; i < n; i++)
      rr += r[i]*r[i];
    return rr;
  }
  
  enum maskType{
    J,
    K};
  
    public static int[] createMask(float[] im, int n12){
      return createMask(im, n12, maskType.J);
    }
    
  public static int[] createMask(float[] im, int n12,  maskType M){
    int num = 0;
    for(int i = 0; i < n12; i++)
      if ( Float.isNaN( im[i] ) || im[i] == 0 ){
        num++;
      }
    int[] _maskK = new int[num];
    num = n12-num;
    int[] _maskJ = new int[num];
    for(int i = 0, j = 0, k = 0; i < n12; i++)
      if ( !Float.isNaN( im[i] ) && im[i] != 0 ){
        _maskJ[j] = i;
        j++;
      }
      else{
        im[i] = 0;
        _maskK[k] = i;
        k++;
      }
    if(M == maskType.J)
      return _maskJ;
    else if(M == maskType.K)
      return _maskK;
    return null;
  }
  
  
  public static void helixConv(
      HelixFilter hf, float[] x, float[] y, boolean doParallel)
  {
    zero(y);
    if (doParallel)
      helixConvParallel(hf,x,y);
    else
      helixConv(hf,x,y);
  }
  
  public static void helixConv(HelixFilter hf, float[] x, float[] y) {
    // pad zeros to the end of input array x if nx < ny 
    int nx = x.length;
    int ny = y.length;
    int nh = hf.nh;
    if (nx<ny) {
      x = Arrays.copyOf(x,ny);
    }
    boolean checkMissing = hf.mis==null ? false : true;
    if (!checkMissing) {
      for (int i=0, ix=0; i<nh; i++) {
        int lag = hf.lag[i];
        float c = hf.coef[i];
        for (int iy=lag; iy<ny; iy++) {
          ix = iy-lag;
          y[iy] += (x[ix]*c);
        }
      }  
    } else {
      for (int i=0, ix=0; i<nh; i++) {
        int lag = hf.lag[i];
        float c = hf.coef[i];
        for (int iy=lag; iy<ny; iy++) {
          if (hf.mis[iy]) continue;
          ix = iy-lag;
          y[iy] += (x[ix]*c);
        }
      }
    }
    
  }
  
  /*Using Jon Clarebout's idea in his book 
  * Image Estimation By example: Geophysical Sounding Image Construction  Chapter 4*/
  public static void helixConvParallel(
      final HelixFilter hf, final float[] x, final float[] y)
  {
    // pad zeros to the end of input array x if nx < ny 
    final int nx = x.length;
    final int ny = y.length;
    int nh = hf.nh;
    final float[]  xx;
    // if the output array is longer than input
    if( nx < ny)
      xx = Arrays.copyOf(x, ny);
    else
      xx = x;
    if (hf.mis==null){
//      for(int i = 0; i < nh; i++){
//        final int lag = hf.lag[i];
//        final float coef = hf.coef[i];
      Parallel.loop(nh,new Parallel.LoopInt() {      
        public void compute(int i) {
//          Parallel.loop(lag,ny,new Parallel.LoopInt() {
//            public void compute(int iy) {        
          int lag = hf.lag[i];
          float coef = hf.coef[i];
          for (int iy=lag; iy<ny; iy++)
            y[iy] += (xx[iy-lag]*coef);
        }
      });
    } else {
      for (int i=0; i<nh; i++) {
//      Parallel.loop(nh,new Parallel.LoopInt() {      
//        public void compute(int i) {
        final float coef = hf.coef[i];
        final int lag = hf.lag[i];
        Parallel.loop(lag,ny,new Parallel.LoopInt() {
          public void compute(int iy) {
//          for (int iy=lag; iy<ny; iy++) {
            if (!hf.mis[iy])
              //TODO write another method with missing data hf.mis 
              y[iy] += (xx[iy-lag]*coef);
//          }
          }
        });
      }
    }
  }
  
//  
//  (ny - start)%
//  Parallel.loop(start,ny,new Parallel.LoopInt() {
//    public void compute(int iy) {        
//      int ix;
//  //for (int iy = hf.lag[i]; iy < ny; iy++) {
//    //if( hf.mis != null && hf.mis[iy]) continue;
//      //TODO write another method with missing data hf.mis 
//    ix = iy - start;
//    y[iy    ] += (xx[ix    ] * coef);
//    y[iy + 1] += (xx[ix - 1] * coef);
//    y[iy + 2] += (xx[ix - 2] * coef);
//    y[iy + 3] += (xx[ix - 3] * coef);
//    }
//  });

  public static void helixConvAdjoint(
      HelixFilter hf, float[] x, float[] y, boolean doParallel)
  {
    zero(x);
    if (doParallel)
      helixConvAdjointParallel(hf,x,y);
    else
      helixConvAdjoint(hf,x,y);
  }
  
  public static void helixConvAdjointParallel(
      final HelixFilter hf, final float[] x, final float[] y)
  {
    // pad zeros to the end of input array x if nx < ny 
    final int nx = x.length;
    final int ny = y.length;
    int nh = hf.nh;
    final int nxm1 = nx-1;

    if (hf.mis==null ) {
      Parallel.loop(nh,new Parallel.LoopInt() {      
        public void compute(int i) {
          //for(int i = 0; i < nh; i++){
          int lag = hf.lag[i];
          float coef = hf.coef[i];
          for (int iy=lag; iy<ny; iy++) {
            int ix = iy-lag;
            if (ix>nxm1) continue;
            x[ix] += (y[iy]*coef);
          }
//          int end = (ny-lag<nx) ? (ny-lag) : nx;
          //Parallel.loop(end,new Parallel.LoopInt() {
          //public void compute(int ix) {        
//          for (int ix=0; ix<end; ix++) {
//            x[ix] += (y[ix+lag] * coef);
//          }
        }
      });
    }
    else{
      //for(int i = 0; i < nh; i++){
      Parallel.loop(nh,new Parallel.LoopInt() {      
        public void compute(int i) {
          int lag = hf.lag[i];
          float coef = hf.coef[i];
          for (int iy=lag; iy<ny; iy++) {
            if (hf.mis[iy]) continue;
            int ix = iy-lag;
            if (ix>nxm1) continue;
            x[ix] += (y[iy]*coef);  
          }
//          int end = min(ny-lag,nx);
          //          Parallel.loop(end,new Parallel.LoopInt() {
          //            public void compute(int ix) {        
          
//          for (int ix=0; ix<end; ix++) {
//            int iy = ix+lag;
//            if(!hf.mis[iy])
//              x[ix] += (y[iy]*coef);
//          }
        }
      });
    }
  }

  public static void helixConvAdjoint(HelixFilter hf, float[] x, float[] y) {
    // pad zeros to the end of input array x if nx < ny 
    int nx = x.length;
    int ny = y.length;
    int nh = hf.nh;
    int nxm1 = nx-1;
    boolean checkMissing = hf.mis==null ? false : true;
    if (!checkMissing) {
      for (int i=0, ix=0; i<nh; i++) {
        int lag = hf.lag[i];
        float c = hf.coef[i];
        for (int iy=lag; iy<ny; iy++) {
          ix = iy-lag;
          if (ix>nxm1) continue;
          x[ix] += (y[iy]*c);  
        }
      }  
    } else {
      for (int i=0, ix=0; i<nh; i++) {
        int lag = hf.lag[i];
        float c = hf.coef[i];
        for (int iy=lag; iy<ny; iy++) {
          if (hf.mis[iy]) continue;
          ix = iy-lag;
          if (ix>nxm1) continue;
          x[ix] += (y[iy]*c);  
        }
      }
    }
    
  }

    
//    for(int i = 0; i < nh; i++){
//      final int start = hf.lag[i];
//      final float coef = hf.coef[i];
//      Parallel.loop(start,ny,new Parallel.LoopInt() {
//        public void compute(int iy) {        
//          int ix;
            //k=%4; 
//          for(int iy = hf.lag[i]; iy < k; iy++){
//            ix = iy - start;
//            if(ix<nx) 
//    x[ix] += (y[iy] * coef);
//              }
//        //for (int iy = hf.lag[i]; iy < ny; iy++) {
//          ix = iy - start;
//          if(ix < nx) 
//            x[ix]     += (y[iy    ] * coef);
//            x[ix + 1] += (y[iy + 1] * coef);
//            x[ix + 2] += (y[iy + 2] * coef);
//            x[ix + 3] += (y[iy + 3] * coef);    
//        }
//      });
//    }
//  }
    
  public static void helixDecon(
      HelixFilter hf, float[] x, float[] y, String zero)
  { 
    zero(y);
    helixDecon(hf, x, y);  
  }
  
  public static void helixDecon( HelixFilter hf, float[] x, float[] y){ 
//     be aware that this take all the diagnals (The first filter coef) = 1 as granted
//     this method exclude the first 1  
//     can be in place, but can't be parallel
//     "any input affects all subsequent outputs"
    int nx = x.length;
    int ny = y.length;
    float[] tmp = new float[nx];
    for (int iy=0, ix=0; iy<ny; iy++) {
      tmp[iy] = x[iy];
      for (int i=1; i<hf.nh; i++) {
        ix = iy - hf.lag[i];
        if(ix < 0) continue;
        tmp[iy] -= (hf.coef[i]*tmp[ix]);
      }
      y[iy] = tmp[iy];
    }
  }
  
  public static void helixDeconAdjoint(
      HelixFilter hf, float[] x, float[] y, String zero)
  { 
    zero(x);
    helixDeconAdjoint(hf, x, y);  
  }

  public static void helixDeconAdjoint(HelixFilter hf, float[] x, float[] y){ 
    // be aware that this take all the diagnals (The first filter coef) = 1 as granted
    // this method exclude the first 1
    //can be in place
    int nx = x.length;
    int ny = y.length;
    float[] tmp = new float[nx];
    for (int ix=nx-1, iy=0; ix>-1; ix--) {
      tmp[ix] = y[ix];
      for (int i=1; i<hf.nh; i++){
        iy = ix + hf.lag[i];
        if(iy > ny-1) continue;
        tmp[ix] -= (hf.coef[i]*tmp[iy]);
      }
      x[ix] = tmp[ix];
    }
  }
  
  public static int[] oneToND(int[] n, int i){
    // convert 1D linear index to ND indexes
    int dim = n.length;
    int[] ni = new int[dim];
    for(int d = 0; d < dim; d++){
      ni[d] = i % n[d];
      i /= n[d];
    }
    return ni;
  }

  public static int NToOneD(int[] n, int[] ni){
    // convert ND indexs to 1D linear index
    int dim = n.length;
    int i = ni[0]; 
    int nn = n[0];
    for(int d = 1; d < dim; d++){
      i += ni[d]*nn;
      nn *= n[d];
    }
    return i;
  }
  
//  public static float[] filterEstimate(HelixFilter h, float[] x, float[] y) {
//    //pad zeros to the end of input array x if nx < ny 
//    int nh  = h.lag.length;
//    int nx = x.length;
//    int ny = y.length;
//    float[] a = new float[nh];
//    if( nx < ny){
//      x = Arrays.copyOf(x, ny);
//    }
//    for(int i = 0; i < nh; i++){
//      for(int iy = h.lag[i], ix = 0; iy < ny; iy++){
//        if(h.mis[iy]) continue;
//        ix = iy - h.lag[i];
//        a[i] += x[ix]*y[iy];
//      }
//    }
//    return a;
//  }

//  public static float[] filterEstimate(
//      final HelixFilter h, final float[] x, final float[] y) 
//  {
//    //pad zeros to the end of input array x if nx < ny 
//    final int nh  = h.lag.length;
//    final int nx = x.length;
//    final int ny = y.length;
//    final float[] xx;
//    final float[] a = new float[nh];
//    if( nx < ny){
//      xx = Arrays.copyOf(x, ny);
//    }
//    else{
//      xx = x;
//    }
//    Parallel.loop(nh,new Parallel.LoopInt() {      
//      public void compute(int i) {
//        //for(int i = 0; i < nh; i++){
//        //final int j = i;
//        int lag = h.lag[i];
//        for(int iy = lag; iy < ny; iy++){
//          int ix;
//          if(!h.mis[iy]){
//            ix = iy - lag;
//            a[i] += xx[ix]*y[iy];
//          }
//        }
//      }
//    });
//    return a;
//  }
  
}
