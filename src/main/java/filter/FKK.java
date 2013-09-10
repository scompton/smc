package filter;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JProgressBar;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Fft;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.ZeroMask;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.DRectangle;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.Projector;
import edu.mines.jtk.mosaic.Tile;

import viewer.Viewer2D;
import viewer.Viewer3P;

import static edu.mines.jtk.util.ArrayMath.*;

public class FKK {

  public FKK(float[][][] f) {
    this(new Sampling(f[0][0].length),new Sampling(f[0].length),
        new Sampling(f.length),f);
  }

  public FKK(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    _f = f;
    _timeView = new Viewer3P(s1,s2,s3,_f);
    _timeView.getViewerFrame().setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    _timeView.setHeightMinimum(0,320);
    _timeView.setSize(WIDTH,HEIGHT);
    _timeView.show();
    goFft();
  }

  public void setZeroMask(ZeroMask zm) {
    _zm = zm;
  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    if (args.length!=4) {
      print("Usage: java fkk.FKK file path n1 n2 n3");
      System.exit(0);
    }
    int n1 = Integer.parseInt(args[1]);
    int n2 = Integer.parseInt(args[2]);
    int n3 = Integer.parseInt(args[3]);
    try {
      float[][][] f3 = new float[n3][n2][n1];
      ArrayInputStream ais = new ArrayInputStream(args[0]);
      ais.readFloats(f3);
      ais.close();
      new FKK(f3);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private 
  private final float[][][] _f;
  private PixelsView _fpv;
  private PointsView _pt00;
  private ZeroMask _zm;
  private Viewer3P _timeView;
  private static final double TWO_PI = 2.0*PI;
  private static final int WIDTH = 800;
  private static final int HEIGHT = 1000;

  private void goFft() {
    int n3 = _f.length;
    int n2 = _f[0].length;
    int n1 = _f[0][0].length;
    Fft fft = new Fft(n1,n2,n3);
    fft.setPadding1(n1);
    fft.setPadding2(n2);
    fft.setPadding3(n3);
    Sampling sf1 = fft.getFrequencySampling1();
    Sampling sf2 = fft.getFrequencySampling2();
    Sampling sf3 = fft.getFrequencySampling3();
    int nk2 = sf2.getCount();
    int nk3 = sf3.getCount();
    sf2 = new Sampling(nk2,sf2.getDelta(),-0.5);
    sf3 = new Sampling(nk3,sf3.getDelta(),-0.5);
    float[][][] forward = fft.applyForward(_f);
    print("n3="+forward.length+", n2="+forward[0].length+", n1="+forward[0][0].length);
//    print("nsf2="+nk2+", nsf1="+nk1);
    sortTraces(nk2,nk3,forward);
    Viewer3P v = new Viewer3P(sf1,sf2,sf3,cabs(forward));
    v.setColorModel1(ColorMap.JET);
    v.setPercentiles1(2.0f,99.9f);
    addPickListener(sf1,sf2,sf3,v,forward,fft);
    v.setLimits1(sf1.getFirst(),sf1.getLast());
    v.setLimits2(sf2.getFirst(),sf2.getLast());
    v.setLimits3(sf3.getFirst(),sf3.getLast());
    v.setLabel2("k2");
    v.setLabel3("k3");
    v.setSize(WIDTH,HEIGHT);
    v.setHeightElastic(0,50);
    v.setWidthElastic(0,100);
    v.setWidthElastic(1,100);
    v.show();
  }

  private void addPickListener(
      final Sampling sf1, final Sampling sf2, final Sampling sf3,
      final Viewer3P v, final float[][][] f, final Fft fft)
  {
    final int n2 = sf2.getCount();
    final int n3 = sf3.getCount();
    final double[] v1 = sf1.getValues();
    final double[] v2 = sf2.getValues();
    final double[] v3 = sf3.getValues();
    PlotPanelPixels3 ppp3 = v.getPlotPanelPixels3();
    final Tile tile00 = ppp3.getTile(0,0);
    tile00.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        Dimension dim = tile00.getSize();
        Projector hp = tile00.getHorizontalProjector();
        Projector vp = tile00.getVerticalProjector();
        DRectangle vr = tile00.getViewRectangle();
        int x = e.getX();
        int y = e.getY();
        double wx = hp.v(vr.width *((double)x/dim.width )+vr.x);
        double wy = vp.v(vr.height*((double)y/dim.height)+vr.y);
        int[] slices = v.getSlices();
        int i1 = slices[0];
        double r = sqrt(wx*wx+wy*wy);
        double m = sf1.getValue(i1)/r;
        float[][][] c = makeCircles(sf1,sf3,m);
        v.addPoints00(c,Color.WHITE,2.0f);
        print("applying filter...");
        float[][][] mask = getMask(v1,v2,v3,m);
        float[][][] filtered = mul(f,mask);
        v.addPixels(cabs(filtered));
        v.setColorModel2(ColorMap.JET);
        v.setPercentiles2(2.0f,99.9f);
//        sortTraces(n2,n3,filtered);
//        print("doing inverse fft...");
//        setData(fft.applyInverse(filtered));
        super.mouseClicked(e);
      }
    });
//    tile.addMouseListener(new MouseAdapter() {
//      @Override
//      public void mouseClicked(MouseEvent e) {
//        Dimension dim = tile.getSize();
//        Projector hp = tile.getHorizontalProjector();
//        Projector vp = tile.getVerticalProjector();
//        DRectangle vr = tile.getViewRectangle();
//        int x = e.getX();
//        int y = e.getY();
//        double wx = hp.v(vr.width *((double)x/dim.width )+vr.x);
//        double wy = vp.v(vr.height*((double)y/dim.height)+vr.y);
//        double m = wy/wx;
//        float[][] mask = getMask(sf1,sf2,m);
//        float[][] filtered = mul(f,mask);
//        pv.set(cabs(filtered));
//        for (int i=0; i<n2; i++)
//          p[i] = (float)abs(sf2.getValue(i)*m);
//        if (_pv==null) {
//          _pv = v.addPoints(p,k,"p");
//          _pv.setLineColor(Color.WHITE);
//        } else {
//          _pv.set(p,k);
//        }
//        float[][] fc = copy(filtered);
//        int h2 = n2/2;
//        copy(fc[0].length,h2,0,h2,fc,0, 0,filtered);
//        copy(fc[0].length,h2,0, 0,fc,0,h2,filtered);
//        setData(fft.applyInverse(filtered));
//        super.mouseClicked(e);
//      }
//    });
  }

  private float[][][] getMask(double[] v1, double[] v2, double[] v3, double m) {
    int n1 = v1.length;
    int n2 = v2.length;
    int n3 = v3.length;
    float[][][] mask = new float[n3][n2][n1*2];
    for (int i1=0; i1<n1; i1++) {
      int ic = i1*2;
      double r= v1[i1]/m;
      for (int i3=0; i3<n3; i3++) {
        double k3 = v3[i3];
        double k3k3 = k3*k3;
        for (int i2=0; i2<n2; i2++) {
          double k2 = v2[i2];
          float v = sqrt(k3k3+k2*k2)<r ? 1.0f : 0.0f;
          mask[i3][i2][ic  ] = v;
          mask[i3][i2][ic+1] = v;
        }
      }
    }
    return mask;
  }

  private float[][][] makeCircles(Sampling sf1, Sampling sf3, double m) {
    int n1 = sf1.getCount();
    double d3 = sf3.getDelta();
    float[][][] c = new float[n1][][];
    for (int i1=0; i1<n1; i1++) {
      double y = sf1.getValue(i1);
      double r = y/m;
      double rr = r*r;
      int hc = (int)ceil(2.0*DBL_PI*r/(d3*0.5))+1;
      int nc = hc*2-1;
      double dc = hc==1 ? 0.0: r*2.0/(hc-1);
      float[] xa = new float[nc];
      float[] ya = new float[nc];
      int i = 0;
      for (; i<hc; i++) {
        double k3 = -r+i*dc;
        double v = rr-k3*k3;
        v = v<0 ? 0.0 : v; 
        xa[i] = (float)sqrt(v);
        ya[i] = (float)k3;
      }
      for (int ih=hc-2; ih>=0; ih--, i++) {
        xa[i] = -xa[ih];
        ya[i] =  ya[ih];
      }
      c[i1] = new float[][]{xa,ya};
    }
    return c;
  }

  private void setData(float[][][] f) {
    if (_zm!=null) _zm.apply(0.0f,f);
    _timeView.addPixels(f);
  }

  private void sortTraces(int n2, int n3, float[][][] f) {
    float[][][] fc = copy(f);
    int h2 = n2/2;
    int h3 = n3/2;
    for (int i2=0; i2<n2; i2++) {
      int i2s = (h2+i2)%n2;
      for (int i3=0; i3<n3; i3++) {
        f[i3][i2] = fc[(h3+i3)%n3][i2s];
      }
    }
  }

  private static void print(String s) {
    System.out.println(s);
  }

}