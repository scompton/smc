package test;

import static edu.mines.jtk.util.ArrayMath.*;

import java.awt.Color;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.EigenTensors2;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Quantiler;

public class Util {

  public static void plot(float[][] S, String title, float min, float max){
    SimplePlot sp = new SimplePlot(); 
    PixelsView pv = sp.addPixels(S);
    sp.addColorBar();
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP); //X1DOWN_X2RIGHT
    pv.setColorModel(ColorMap.JET);
    sp.setTitle(title);
    ColorBar cb = sp.addColorBar();
    cb.setWidthMinimum(150);
    pv.setClips(min,max);
  }
  public static void plot(float[][] S, String title){
    SimplePlot sp = new SimplePlot(); 
    PixelsView pv = sp.addPixels(S);
    sp.addColorBar();
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP); //X1DOWN_X2RIGHT
    pv.setColorModel(ColorMap.JET);
    sp.setTitle(title);
    ColorBar cb = sp.addColorBar();
    cb.setWidthMinimum(150);
  }
 
  public static void plot(float[] S, int n1, int n2, String title){
    float[][] s = new float[n2][n1];
    Operations.convert(S,n1,n2,s);
    plot(s, title);
  }
  
  public static void getTensorEllipses(int n1, int n2, 
      int ns, EigenTensors2 et, float[][] image) {
    int nt = 51;
    int m1 = (n1-1)/ns;
    int m2 = (n2-1)/ns;
    int j1 = (n1-1-(m1-1)*ns)/2;
    int j2 = (n2-1-(m2-1)*ns)/2;
    int nm = m1*m2;
    float[][] sm = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1) {
  float[] s = et.getEigenvalues(i1,i2);
  sm[i2][i1] = s[1];
      }
    sm = copy(m1,m2,j1,j2,ns,ns,sm);
    float smq = Quantiler.estimate(0.98f,sm);
    double r = 0.45*ns/sqrt(smq);
    float[][] x1 = new float[nm][nt];
    float[][] x2 = new float[nm][nt];
    double dt = 2.0*PI/(nt-1);
    double ft = 0.0f;
    for (int i2=j2,k2=0,im=0; i2<n2 && k2<m2; i2+=ns,++k2)
      for (int i1=j1,k1=0; i1<n1 && k1<m1; i1+=ns,++k1,++im) {
        float[] u = et.getEigenvectorU(i1,i2);
        float[] s = et.getEigenvalues(i1,i2);
        double u1 = u[0], u2 = u[1];
        double du = s[0], dv = s[1];
        du = min(du,smq);
        dv = min(dv,smq);
        double a = r*sqrt(dv), b = r*sqrt(du);
        for (int it=0; it<nt; ++it) {
          double t = ft+it*dt;
          double cost = cos(t);
          double sint = sin(t);
          x1[im][it] = (float)(i1+b*cost*u1-a*sint*u2);
          x2[im][it] = (float)(i2+a*sint*u1+b*cost*u2);
        }
      }
    
    SimplePlot plot= new SimplePlot();
    plot.setTitle("Diffusion Tensors");
    PixelsView pv = plot.addPixels(image);
    ColorBar cb = plot.addColorBar();
    cb.setWidthMinimum(150);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.JET);
    PointsView _tensorsView = new PointsView(x2,x1);
    _tensorsView.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    _tensorsView.setLineColor(Color.YELLOW);
    plot.getPlotPanel().getTile(0, 0).addTiledView(_tensorsView);
  }
  
}
