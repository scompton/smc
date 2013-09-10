package filter;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Fft;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.ZeroMask;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;
import edu.mines.jtk.mosaic.DRectangle;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.Projector;
import edu.mines.jtk.mosaic.Tile;

import viewer.Viewer2D;
import viewer.ViewerFrame;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Simple interactive 2D F-K filtering.
 */
public class FK {

  /**
   * Constructs FK instance with default {@link Sampling}s.
   * @param f
   */
  public FK(float[][] f) {
    this(new Sampling(f[0].length),new Sampling(f.length),f);
  }

  /**
   * Constructs FK instance using provided {@link Sampling}s.
   * @param s1
   * @param s2
   * @param f
   */
  public FK(Sampling s1, Sampling s2, float[][] f) {
    _f = f;
    _timeView = new Viewer2D();
    _timeView.addPixels(s1,s2,f,"f");
    _timeView.getViewerFrame().setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    _timeView.setSize(WIDTH,HEIGHT);
    addWriteOption();
    _timeView.show();
    goFft();
  }

  /**
   * Set a {@link ZeroMask} to apply to the filtered data. 
   * @param zm
   */
  public void setZeroMask(ZeroMask zm) {
    _zm = zm;
  }

  /**
   * Runs FK filter for a file on disk.
   * @param args
   */
  public static void main(String[] args) {
    if (args.length!=5) {
      print("Usage: java fkk.FKK file path n1 n2 n3 i3");
      System.exit(0);
    }
    int n1 = Integer.parseInt(args[1]);
    int n2 = Integer.parseInt(args[2]);
    int n3 = Integer.parseInt(args[3]);
    int i3 = Integer.parseInt(args[4]);
    try {
      float[][][] f3 = new float[n3][n2][n1];
      ArrayInputStream ais = new ArrayInputStream(args[0]);
      ais.readFloats(f3);
      ais.close();
      new FK(f3[i3]);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private 
  private final float[][] _f;
  private float[][] _filtered = null;
  private PixelsView _fpv;
  private PointsView _pv;
  private ZeroMask _zm;
  private Viewer2D _timeView;
  private static final double TWO_PI = 2.0*PI;
  private static final int WIDTH = 800;
  private static final int HEIGHT = 1000;

  private void goFft() {
    int n2 = _f.length;
    int n1 = _f[0].length;
    Fft fft = new Fft(n1,n2);
    fft.setPadding1(n1);
    fft.setPadding2(n2);
    Sampling sf1 = fft.getFrequencySampling1();
    Sampling sf2 = fft.getFrequencySampling2();
    int nk = sf2.getCount();
    sf2 = new Sampling(nk,sf2.getDelta(),-0.5);
    float[][] forward = fft.applyForward(_f);
    sortTraces(nk,forward);
    Viewer2D v = new Viewer2D();
    PixelsView pv = v.addPixels(sf1,sf2,cabs(forward),"forward");
    pv.setColorModel(ColorMap.JET);
    pv.setPercentiles(2.0f,99.9f);
    Tile tile = v.getPlotPanel().getTile(0,0);
    addPickListener(sf1,sf2,tile,v,pv,forward,fft);
    v.setVLimits(0.0,0.5);
    v.setSize(WIDTH,HEIGHT);
    v.show();
  }

  private void addWriteOption() {
    final ViewerFrame vf = _timeView.getViewerFrame();
    Action saveFilteredAction = new AbstractAction("Save filtered data") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        if (_filtered==null) {
          JOptionPane.showMessageDialog(vf,
              "Cannot save, data has not been filtered.");
          return;
        }
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.showSaveDialog(_timeView.getViewerFrame());
        File file = fc.getSelectedFile();
        if (file!=null) {
          String filename = file.getAbsolutePath();
          if (!filename.toLowerCase().endsWith(".dat"))
            filename = filename+".dat";
          boolean success = true;
          try {
            ArrayOutputStream aos = new ArrayOutputStream(filename);
            aos.writeFloats(_filtered);
            aos.close();
          } catch (FileNotFoundException e) {
            success = false;
            e.printStackTrace();
          } catch (IOException e) {
            success = false;
            e.printStackTrace();
          }
          if (!success) {
            JOptionPane.showMessageDialog(vf,"Failed to save filtered data.");
          }
        }
      }
    };
    JMenuItem saveFiltered = new JMenuItem(saveFilteredAction);
    _timeView.getViewerFrame().addToMenu(saveFiltered);
  }

  private void addPickListener(
      final Sampling sf1, final Sampling sf2, final Tile tile, final Viewer2D v, 
      final PixelsView pv, final float[][] f, final Fft fft)
  {
    final int nk = sf2.getCount();
    final double[] kd = sf2.getValues();
    final float[] k = new float[nk];
    for (int i=0; i<nk; i++)
      k[i] = (float)kd[i];
    final float[] p = new float[nk];
    tile.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        Dimension dim = tile.getSize();
        Projector hp = tile.getHorizontalProjector();
        Projector vp = tile.getVerticalProjector();
        DRectangle vr = tile.getViewRectangle();
        int x = e.getX();
        int y = e.getY();
        double wx = hp.v(vr.width *((double)x/dim.width )+vr.x);
        double wy = vp.v(vr.height*((double)y/dim.height)+vr.y);
        double m = wy/wx;
        float[][] mask = getMask(sf1,sf2,m);
        float[][] filtered = mul(f,mask);
        pv.set(cabs(filtered));
        for (int i=0; i<nk; i++) p[i] = (float)abs(kd[i]*m);
        if (_pv==null) {
          _pv = v.addPoints(p,k,"p");
          _pv.setLineColor(Color.WHITE);
        } else {
          _pv.set(p,k);
        }
        sortTraces(nk,filtered);
        setData(fft.applyInverse(filtered));
        super.mouseClicked(e);
      }
    });
  }

  private float[][] getMask(Sampling sf1, Sampling sf2, double m) {
    int n1 = sf1.getCount();
    int n2 = sf2.getCount();
    double d = sf2.getDelta();
    double[] values =  sf2.getValues();
    float[][] mask = new float[n2][n1*2];
    for (int i1=0; i1<n1; i1++) {
      double k = abs(sf1.getValue(i1)/m);
      int ic = i1*2;
      float[] t = taper(k,d,values);
      for (int i2=0; i2<n2; i2++) {
        mask[i2][ic  ] = t[i2];
        mask[i2][ic+1] = t[i2];
      }
    }
    return mask;
  }

  private float[] taper(double k, double d, double[] values) {
    int n = values.length;
    float[] t = new float[n];
    double f = 0.5/(2*k+d);
    for (int i=0; i<n; i++) {
      if (values[i]<-k || values[i]>k) {
        t[i] = 0.0f;
        continue;
      }
      t[i] = (float)cos(TWO_PI*f*values[i]);
    }
    return t;
  }

  private void setData(float[][] f) {
    _filtered = f;
    if (_fpv==null) _fpv = _timeView.addPixels(_filtered,"filtered");
    if (_zm!=null) _zm.apply(0.0f,_filtered);
    _fpv.set(_filtered);
  }

  private void sortTraces(int n2, float[][] f) {
    int h2 = n2/2;
    float[][] fc = copy(f);
    for (int i2=0; i2<n2; i2++) {
      f[i2] = fc[(i2+h2)%n2];
    }
  }

  private static void print(String s) {
    System.out.println(s);
  }

}