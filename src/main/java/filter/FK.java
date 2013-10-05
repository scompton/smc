package filter;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

import viewer.*;

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
    _s1 = s1;
    _s2 = s2;
    _f = f;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Check.argument(n1==_f[0].length,"s1 consistent with array f");
    Check.argument(n2==_f.length,"s2 consistent with array f");
    _fft = new Fft(n1,n2);
    _fft.setPadding1(n1);
    _fft.setPadding2(n2);
    _sf1 = _fft.getFrequencySampling1();
    _sf2 = _fft.getFrequencySampling2();
    _sf2 = new Sampling(_sf2.getCount(),_sf2.getDelta(),-0.5);
    _forward = _fft.applyForward(_f);
    sortTraces(_forward);
  }

  public float[][] filter(double m) {
    float[][] mask = getMask(_sf1,_sf2,m);
    return mul(_forward,mask);
  }

  public float[][] inverse(float[][] filtered) {
    sortTraces(filtered);
    float[][] inverse = _fft.applyInverse(filtered);
    if (_zm!=null) _zm.apply(0.0f,inverse);
    return inverse;
  }

  public void designFilter() {
    _timeView = new Viewer2D();
    _timeView.addPixels(_s1,_s2,_f,"f");
    _timeView.getViewerFrame().setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    _timeView.setSize(WIDTH,HEIGHT);
    addWriteOption();
    _timeView.show();

    Viewer2D v = new Viewer2D();
    PixelsView pv = v.addPixels(_sf1,_sf2,cabs(_forward),"forward");
    pv.setColorModel(ColorMap.JET);
    pv.setPercentiles(2.0f,99.9f);
    Tile tile = v.getPlotPanel().getTile(0,0);
    addPickListener(tile,v,pv);
    v.setVLimits(0.0,0.5);
    v.setSize(WIDTH,HEIGHT);
    v.show();
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
  private final float[][] _f; // Input data
  private final Sampling _s1, _s2; // Data sampling
  private Sampling _sf1, _sf2; // Frequency sampling
  private float[][] _filtered = null;
  private float[][] _forward = null;
  private Fft _fft;
  private PixelsView _fpv;
  private PointsView _pv;
  private ZeroMask _zm;
  private Viewer2D _timeView;
  private static final double TWO_PI = 2.0*PI;
  private static final int WIDTH = 800;
  private static final int HEIGHT = 1000;

  private void goFft() {
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
      final Tile tile, final Viewer2D v, final PixelsView pv)
  {
    final int nk = _sf2.getCount();
    final double[] kd = _sf2.getValues();
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
        print("slope: "+m);
        float[][] filtered = filter(m);
        pv.set(cabs(filtered));
        for (int i=0; i<nk; i++) p[i] = (float)abs(kd[i]*m);
        if (_pv==null) {
          _pv = v.addPoints(p,k,"p");
          _pv.setLineColor(Color.WHITE);
        } else {
          _pv.set(p,k);
        }
        setData(inverse(filtered));
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

  private void sortTraces(float[][] f) {
    int n2 = _sf2.getCount();
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