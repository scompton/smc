package viewer;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.util.*;

/**
 * Wraps a PlotPanel with one Tile for 2D or 3D pixels and
 * includes convenient options for interactively editing the plot.
 * This class inlcludes methods to overlay points or pixels on
 * a base PixelsView.
 * </p>
 * For 3D arrays a slider allows panning through the volume.
 */
public class Viewer2D {

  /**
   * Constructs a pixels view of {@code f}.
   * @param f the pixel values.
   */
  public Viewer2D() {
    this((Orientation)null);
  }

  /**
   * Constructs a pixels view of {@code f}, with specified
   * orientation.
   * @param f the pixel values.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(Orientation orientation) {
    _orientation = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    _vf = new ViewerFrame(_orientation); // Make empty panel.
    _pp = _vf.getPlotPanel();
  }

  /**
   * Get the {@link ViewerFrame}.
   * @return the {@link ViewerFrame}.
   */
  public ViewerFrame getViewerFrame() {
    return _vf;
  }

  public PlotPanel getPlotPanel() {
    return _pp;
  }

  public Sampling getSampling1() {
    return _s1;
  }
  public Sampling getSampling2() {
    return _s2;
  }
  public Sampling getSampling3() {
    return _s3;
  }
  public double getHMin() {
    return _hMin;
  }
  public double getHMax() {
    return _hMax;
  }
  public double getVMin() {
    return _vMin;
  }
  public double getVMax() {
    return _vMax;
  }

  public PixelsView addPixels(double[][] f, String label) {
    Sampling s1 = new Sampling(f[0].length);
    Sampling s2 = new Sampling(f.length);
    return addPixels(s1,s2,f,label);
  }

  public PixelsView addPixels(
      Sampling s1, Sampling s2, double[][] f, String label)
  {
    return addPixels(s1,s2,convertToFloat(f),label);
  }

  public PixelsView addPixels(float[][] f, String label) {
    Sampling s1 = new Sampling(f[0].length);
    Sampling s2 = new Sampling(f.length);
    return addPixels(s1,s2,f,label);
  }

  public PixelsView addPixels(
      Sampling s1, Sampling s2, float[][] f, String label)
  {
    if (_s1==null && _s2==null) {
      _s1 = s1;
      _s2 = s2;
      if (_orientation==Orientation.X1DOWN_X2RIGHT) {
        _hMin = _s2.getFirst();
        _hMax = _s2.getLast();
        _vMin = _s1.getFirst();
        _vMax = _s1.getLast();
      } else {
        _hMin = _s1.getFirst();
        _hMax = _s1.getLast();
        _vMin = _s2.getFirst();
        _vMax = _s2.getLast();
      }
    } else {
      Check.argument(f.length==_s2.getCount(),
          "f.length is not consistent with sampling");
      Check.argument(f[0].length==_s1.getCount(),
          "f[0].length is not consistent with sampling");
    }
    PixelsView pv = _pp.addPixels(_s1,_s2,f);
    _vf.addOptions(new PixelsView[]{pv},label);
    _pv2Map.put(pv,f);
    updatePixels2();
    return pv;
  }

  public PixelsView addPixels(float[][][] f, String label) {
    Sampling s1 = new Sampling(f[0][0].length);
    Sampling s2 = new Sampling(f[0].length);
    Sampling s3 = new Sampling(f.length);
    return addPixels(s1,s2,s3,f,label);
  }

  public PixelsView addPixels(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, String label)
  {
    if (_s1==null && _s2==null && _s3==null) {
      _s1 = s1;
      _s2 = s2;
      _s3 = s3;
      if (_orientation==Orientation.X1DOWN_X2RIGHT) {
        _hMin = _s2.getFirst();
        _hMax = _s2.getLast();
        _vMin = _s1.getFirst();
        _vMax = _s1.getLast();
      } else {
        _hMin = _s1.getFirst();
        _hMax = _s1.getLast();
        _vMin = _s2.getFirst();
        _vMax = _s2.getLast();
      }
    } else {
      Check.argument(f.length==_s3.getCount(),
          "f.length is not consistent with sampling");
      Check.argument(f[0].length==_s2.getCount(),
          "f[0].length is not consistent with sampling");
      Check.argument(f[0][0].length==_s1.getCount(),
          "f[0][0].length is not consistent with sampling");
    }

    //al panel, displaying the middle frame.
    int n3 = _s3.getCount();
    _i3 = n3/2;
    //_title = String.valueOf(_i3);
    //_pp.setTitle(_title);
    PixelsView pv = _pp.addPixels(_s1,_s2,f[_i3]);
    Clips clips = new Clips(f);
    pv.setClips(clips.getClipMin(),clips.getClipMax());

    _vf.addOptions(new PixelsView[]{pv},label);
    _pv3Map.put(pv,f);
    updatePixels3();
    return pv;
  }

  /**
   * Adds a 2D array to be viewed as 1D slices.
   * @param p
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @param i2Label the label for the slice index displayed in the title,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[][] p, String label, String i2Label) {
    Sampling s1 = new Sampling(p[0].length);
    return addPoints(s1,p,label,i2Label);
  }

  /**
   * Adds a 2D array to be viewed as 1D slices.
   * @param s1 first(fast) dimension sampling.
   * @param p
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @param i2Label the label for the slice index displayed in the title,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(
      Sampling s1, float[][] p, String label, String i2Label)
  {
    if (_s1==null) {
      _s1 = s1;
      if (_orientation==Orientation.X1DOWN_X2RIGHT) {
        _vMin = _s1.getFirst();
        _vMax = _s1.getLast();
      } else {
        _hMin = _s1.getFirst();
        _hMax = _s1.getLast();
      }
    } else {
      Check.argument(p[0].length==_s1.getCount(),
          "p[0].length is not consistent with sampling");
    }
    int n2 = p.length;
    _i2 = n2/2;
    _i2Label = i2Label==null ? "" : i2Label;
    PointsView pt = _pp.addPoints(_s1,p[_i2]);
    _vf.addOptions(pt,label);
    _pt2Map.put(pt,p);
    updatePoints2();
    return pt;
  }

  /**
   * Adds a view of points (x1,x2) for a sampled function x2(x1).
   * @param x2 array of x2 coordinates.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[] x2, String label) {
    Check.argument(_s1.getCount()==x2.length,
        "x2.length is not consistend with sampling");
    PointsView pt = _pp.addPoints(_s1,x2);
    _vf.addOptions(pt,label);
    return pt;
  }

  /**
   * Adds a points view of the arrays x1 and x2 of point (x1,x2) coordinates.
   * @param x1 array of x1 coordinates.
   * @param x2 array of x2 coordinates.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[] x1, float[] x2, String label) {
    PointsView pt = _pp.addPoints(x1,x2);
    _vf.addOptions(pt,label);
    return pt;
  }

  /**
   * Adds a view of arrays of (x1,x2) coordinates for multiple plot segments.
   * The lengths of the specified arrays x1 and x2 must be equal.
   * @param x1 array of arrays of x1 coordinates.
   * @param x2 array of arrays of x2 coordinates.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[][] x1, float[][] x2, String label) {
    PointsView pt = _pp.addPoints(x1,x2);
    _vf.addOptions(pt,label);
    return pt;
  }

  /**
   * Adds a view of points (x1,x2) for a sampled function x2(x1),
   * specified for each slice of a 3D array.
   * @param x2 array of x2 coordinates for each x3.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints3(float[][] x2, String label) {
    Check.argument(_s3.getCount()==x2.length,
      "x2.length is not consistent with sampling");
    Check.argument(_s1.getCount()==x2[0].length,
      "x2[0].length is not consistent with sampling");
    PointsView pt = _pp.addPoints(_s1,x2[_i3]);
    _pt31DMap.put(pt,x2);
    updatePoints3();
    _vf.addOptions(pt,label);
    return pt;
  }

  /**
   * Adds a points view of the arrays x1 and x2 of point (x1,x2) coordinates,
   * specified for each slice of a 3D array.
   * @param x1 array of x1 coordinates.
   * @param x2 array of x2 coordinates.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints3(float[][] x1, float[][] x2, String label) {
    Check.argument(_s3.getCount()==x1.length,
      "x1.length is not consistent with sampling");
    Check.argument(_s3.getCount()==x2.length,
      "x2.length is not consistent with sampling");
    PointsView pt = _pp.addPoints(x2[_i3],x2[_i3]);
    _pt32DMap.put(pt,new float[][][]{x1,x2});
    updatePoints3();
    _vf.addOptions(pt,label);
    return pt;
  }

  /**
   * Adds a view of arrays of (x1,x2) coordinates for multiple plot segments,
   * specified for each slice of a 3D array.
   * The lengths of the specified arrays x1 and x2 must be equal.
   * @param x1 array of arrays of x1 coordinates for each slice of a 3D array.
   * @param x2 array of arrays of x2 coordinates for each slice of a 3D array.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints3(float[][][] x1, float[][][] x2, String label) {
    Check.argument(_s3.getCount()==x1.length,
      "x1.length is not consistent with sampling");
    Check.argument(_s3.getCount()==x2.length,
      "x2.length is not consistent with sampling");
    PointsView pt = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt33DMap.put(pt,new float[][][][]{x1,x2});
    updatePoints3();
    _vf.addOptions(pt,label);
    return pt;
  }

  public void addTensors(EigenTensors2 et2) {
    addTensors(et2,_s2.getCount()/10,Color.YELLOW,1.0f);
  }

  public void addTensors(EigenTensors2 et2, int ne, Color color, float width) {
    TensorsView tensorView = new TensorsView(_s1,_s2,et2);
    TensorsView.Orientation orientation;
    if (_orientation.equals(PlotPanel.Orientation.X1DOWN_X2RIGHT))
      orientation = TensorsView.Orientation.X1DOWN_X2RIGHT;
    else
      orientation = TensorsView.Orientation.X1RIGHT_X2UP;
    tensorView.setOrientation(orientation);
    tensorView.setEllipsesDisplayed(ne);
    tensorView.setLineColor(color);
    Sampling e1 = new Sampling(12,0.25,0.25);
    Sampling e2 = new Sampling(10,0.5,0.0);
    tensorView.setEllipsesDisplayed(e1,e2);
    tensorView.setLineWidth(width);
    tensorView.setScale(1);
    _pp.addTiledView(tensorView);
  }

  public void setTitle(String title) {
    _title = title;
    if (_i3!=Integer.MIN_VALUE)
      _pp.setTitle(_title+" "+_i3);
    else
      _pp.setTitle(_title);
  }

  public void setHLabel(String label) {
    _pp.setHLabel(label);
  }

  public void setVLabel(String label) {
    _pp.setVLabel(label);
  }

  public void setHLimits(double hmin, double hmax) {
    _hMin = hmin;
    _hMax = hmax;
    _pp.setHLimits(hmin,hmax);
  }

  public void setVLimits(double vmin, double vmax) {
    _vMin = vmin;
    _vMax = vmax;
    _pp.setVLimits(vmin,vmax);
  }

  public void setVFormat(String format) {
    _pp.setVFormat(format);
  }

  public void setHFormat(String format) {
    _pp.setHFormat(format);
  }

  public void setHInterval(double interval) {
    _pp.setHInterval(interval);
  }

  public void setVInterval(double interval) {
    _pp.setVInterval(interval);
  }

  public ColorBar addColorBar(String label) {
    return _pp.addColorBar(label);
  }

  public void setColorBarWidthMinimum(int widthMinimum) {
    _pp.setColorBarWidthMinimum(widthMinimum);
  }

  public void setSize(int width, int height) {
    _vf.setSize(width,height);
  }

  public void setFontSizeForPrint(double fontSize, double plotWidth) {
    _vf.setFontSizeForPrint(fontSize,plotWidth);
  }

  public void setFontSizeForSlide(double fracWidth, double fracHeight) {
    _vf.setFontSizeForSlide(fracWidth, fracHeight);
  }

  public void paintToPng(double dpi, double win, String fileName) {
    _vf.paintToPng(dpi,win,fileName);
  }

  public void show() {
    // Add LimitsFrame2D dialog to the options menu.
    JMenuItem changeLimits = new JMenuItem("Change Limits");
    final LimitsFrame2D lf2d = new LimitsFrame2D(this);
    changeLimits.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        lf2d.getFrame().setVisible(true);
      }
    });
    _vf.addToMenu(changeLimits);

    // Possibly add a slider to pan through data
    JSlider slider = null;
    if (_pt2.length>0) { // array of points views
      PointSliderListener sl = new PointSliderListener();
      int n2 = _pt2Map.get(_pt2[0]).length;
      slider = makeSlider(_i2,n2,sl);
    } else if (_pv3Map.size()>0) { // array of pixels views
      PixelSliderListener sl = new PixelSliderListener();
      int n3 = _s3.getCount();
      slider = makeSlider(_i3,n3,sl);
    }

    // add a bottom panel to the frame
    JPanel bottom = new JPanel();
    bottom.setLayout(new BorderLayout());
    if (slider!=null) {
      bottom.add(slider,BorderLayout.NORTH);
    }
    if (_pv2.length>0 || _pv3.length>0) { // add mouse tracking for pixels
      JTextArea text = addMouseTracker();
      bottom.add(text,BorderLayout.SOUTH);
    }
    _vf.add(bottom,BorderLayout.SOUTH);

    _vf.setVisible(true);
  }

  public static void main(String[] args) throws IOException {
    if (args.length != 4) {
      System.out.println("usage: java Viewer datasetPath n1 n2 n3");
      System.exit(0);
    }
    ArrayInputStream ais = new ArrayInputStream(args[0]);
    int n1 = Integer.parseInt(args[1]);
    int n2 = Integer.parseInt(args[2]);
    int n3 = Integer.parseInt(args[3]);
    float[][][] f = new float[n3][n2][n1];
    ais.readFloats(f);
    ais.close();
    Viewer2D v = new Viewer2D();
    v.addPixels(f,"f");
    v.show();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private ViewerFrame _vf;
  private PlotPanel _pp;

  private Map<PixelsView,float[][]> _pv2Map = new HashMap<>();
  PixelsView[] _pv2 = new PixelsView[0];

  private Map<PointsView,float[][]> _pt2Map = new HashMap<>();
  PointsView[] _pt2 = new PointsView[0];

  private Map<PixelsView,float[][][]> _pv3Map = new HashMap<>();
  PixelsView[] _pv3 = new PixelsView[0];

  // For 1D x1 and x2 arrays for each slice of the 3D data.
  private Map<PointsView,float[][]> _pt31DMap = new HashMap<>();
  PointsView[] _pt31D = new PointsView[0];

  // For 2D x1 and x2 arrays for each slice of the 3D data.
  private Map<PointsView,float[][][]> _pt32DMap = new HashMap<>();
  PointsView[] _pt32D = new PointsView[0];

  // For 3D x1 and x2 for each slice of the 3D data.
  private Map<PointsView,float[][][][]> _pt33DMap = new HashMap<>();
  PointsView[] _pt33D = new PointsView[0];

  private Orientation _orientation;
  private Sampling _s1 = null;
  private Sampling _s2 = null;
  private Sampling _s3 = null;
  private double _hMin;
  private double _hMax;
  private double _vMin;
  private double _vMax;
  private String _title = "";
  private String _i2Label = "";
  private int _i3 = Integer.MIN_VALUE;
  private int _i2 = Integer.MIN_VALUE;

  private JSlider makeSlider(int i, int n, ChangeListener cl) {
    DefaultBoundedRangeModel brm = new DefaultBoundedRangeModel(i,0,0,n-1);
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(n/10);
    slider.setMinorTickSpacing(n/150);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(cl);
    return slider;
  }

  /**
   * Adds a panel in the bottom of the ViewerFrame that
   * tracks the mouse location and displays the world
   * coordinates.
   */
  private JTextArea addMouseTracker() {
    final JTextArea text = new JTextArea("World coordinates:");
    text.setEditable(false);
    Mosaic mosaic = _pp.getMosaic();
    int nrows = mosaic.countRows();
    int ncols = mosaic.countColumns();
    for (int irow=0; irow<nrows; irow++) {
      for (int icol=0; icol<ncols; icol++) {
        final Tile tile = mosaic.getTile(irow,icol);
        tile.addMouseMotionListener(new MouseMotionAdapter() {
          @Override
          public void mouseMoved(MouseEvent e) {
            Dimension d = tile.getSize();
            Projector hp = tile.getHorizontalProjector();
            Projector vp = tile.getVerticalProjector();
            DRectangle vr = tile.getViewRectangle();
            int x = e.getX();
            int y = e.getY();
            double wx = hp.v(vr.width *((double)x/d.width)+vr.x);
            double wy = vp.v(vr.height*((double)y/d.height)+vr.y);
            int i1, i2;
            if (_orientation==Orientation.X1DOWN_X2RIGHT) {
              i1 = _s1.indexOfNearest(wy);
              i2 = _s2.indexOfNearest(wx);
            } else {
              i1 = _s1.indexOfNearest(wx);
              i2 = _s2.indexOfNearest(wy);
            }
            float v;
            //FIXME this isn't accurate if there are multiple pixels views.
            if (_pv2.length>0)
              v = (_pv2Map.get(_pv2[0]))[i2][i1];
            else
              v = (_pv3Map.get(_pv3[0]))[_i3][i2][i1];
            text.setText(String.format(
                "World coordinates: x=%6.3f, y=%6.3f, value=%6.8f",wx,wy,v));
            super.mouseMoved(e);
          }
        });
      }
    }
    return text;
  }

  private void updatePixels2() {
    _pv2 = _pv2Map.keySet().toArray(new PixelsView[0]);
  }

  private void updatePixels3() {
    _pv3 = _pv3Map.keySet().toArray(new PixelsView[0]);
  }

  private void updatePoints2() {
    _pt2 = _pt2Map.keySet().toArray(new PointsView[0]);
  }

  private void updatePoints3() {
    _pt31D = _pt31DMap.keySet().toArray(new PointsView[0]);
    _pt32D = _pt32DMap.keySet().toArray(new PointsView[0]);
    _pt33D = _pt33DMap.keySet().toArray(new PointsView[0]);
  }

  private static float[][] convertToFloat(double[][] a) {
    int n2 = a.length;
    float[][] b = new float[n2][];
    for (int i2=0; i2<n2; ++i2) {
      int n1 = a[i2].length;
      b[i2] = new float[n1];
      for (int i1=0; i1<n1; ++i1) {
        b[i2][i1] = (float)a[i2][i1];
      }
    }
    return b;
  }

  private class PointSliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i2 = source.getValue();
      for (PointsView pt : _pt2)
        pt.set(_s1,_pt2Map.get(pt)[i2]);
      _i2 = i2;
      _pp.removeTitle();
      _pp.setTitle(_title+" "+_i2Label+" "+_i2);
      _vf.repaint();
    }
  }

  private class PixelSliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i3 = source.getValue();
      for (PixelsView pv : _pv3)
        pv.set(_pv3Map.get(pv)[i3]);
      for (PointsView pt : _pt31D)
        pt.set(_s1,_pt31DMap.get(pt)[i3]);
      for (PointsView pt : _pt32D) {
        float[][][] x1x2 = _pt32DMap.get(pt);
        pt.set(x1x2[0][i3],x1x2[1][i3]);
      }
      for (PointsView pt : _pt33D) {
        float[][][][] x1x2 = _pt33DMap.get(pt);
        pt.set(x1x2[0][i3],x1x2[1][i3]);
      }
      _i3 = i3;
      _pp.removeTitle();
      _pp.setTitle(_title+" "+_i3);
      _vf.repaint();
    }
  }

}
