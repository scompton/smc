package test;

public class HelixFilter {
  //# of helix filter coefficients, inputs, outputs, and dimensions
  int nh = 1; int nh0,nh1,nh2;
  int _nd, _nm, _dim;  
  int[] _nhf, _n;  
  int[] lag;
  float[] coef;
  boolean[] mis;
 
  
/*  enum fType{
    Laplacian,
    PEF
  };
  // need to know length of both data dimensions and filter 
  public HelixFilter(int[] n){
    this(n, fType.Laplacian);
  }*/
  public HelixFilter(){
    
  }
  
  public HelixFilter(int[] nhf, int[] n, int nd, int nm){
    _dim = n.length;
    _nhf = new int[_dim];
    for(int i = 0; i < _dim; i++){
      nh *= nhf[i];
      _nhf[i] = nhf[i];
    }
    nh -= _nhf[0]/2;
    lag = new int[nh];
    coef = new float[nh];
    _nd = nd;
    _nm = nm;
    _n = n;
    _nhf = nhf;
    helixLag(nhf,n);
  }
  
 public void helixLag(int[] nhf, int[] n){
   nh0 = nhf[0];
   int a1 = (nh0+1)/2;
   for(int i = 0; i < a1; i ++)
     lag[i] = i;     
   if(_dim == 1){
     for(int i = a1; i < nh0; i ++)
       lag[i] = i;
     return;
   }
   else{
     nh1 = nhf[1];
     int ind, ll;
     int l = -a1+nh0%2;
     for(int i = 0; i < nh1-1; i++){
       ind = i*nh0+a1;
       ll = (i+1)*n[0]+l;
       for(int j = 0; j < nh0; j++)
         lag[ind+j] = ll + j; 
     }
     if(_dim==2)     return;
     else{
       nh2 = nhf[2];
       int ind2, lll;
       int n12 = n[1]*n[0];
       int a2 = (nh1+1)/2;
       int a12 = (nh1-1)/2*nh0 + a1;
       l += (-a2 + nh1%2);
       for(int k = 0; k < nh2-1; k++){
         ind2 = k*n12+a12;
         lll = (k+1)*n12;
         for(int i = 0; i < nh1-1; i++){
           ind = ind2 + i*nh0;
           ll = (i+1)*n[0] + l + lll;
           for(int j = 0; j < nh0; j++)
             lag[ind+j] = ll + j;  
         }
       }
     }
   }
 }
       
//    if(type == fType.Laplacian){
//      helixLaplacian(n);
//    }
//    
  
//    else 
//      _filter = new PEF();

  
  // mask helix filter outputs;
  // for completeness, need to be nonzero even beyond the known data

  
  
  
//  
//  public void helixLaplacian (int[] n){
//    nh = 2*nd+1;
//    if(nh == 5)
//      coef = {1.8,-0.6,0.0,-0.2,-0.6};
//      lag[0] = n[0]/2;
//
//  }
//  
  /* Using Jon Clarebout's idea in his book 
   * Image Estimation By example: Geophysical Sounding Image Construction  Chapter 6
   */
  
  

      
  
  
//  public static void main(String[] args) {
//  int[] f = {5,3};
//  int[] d = {150,150};
//  HelixFilter hf = new HelixFilter(f,d,22500, 22500);
//  hf.helixLag(f,d);
//  hf.bound();
//
//  }

}
