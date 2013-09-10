package warp;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;

//import org.eclipse.jface.dialogs.TitleAreaDialog;
//import org.eclipse.swt.SWT;
//import org.eclipse.swt.events.FocusAdapter;
//import org.eclipse.swt.events.FocusEvent;
//import org.eclipse.swt.events.MouseAdapter;
//import org.eclipse.swt.events.MouseEvent;
//import org.eclipse.swt.events.PaintEvent;
//import org.eclipse.swt.events.PaintListener;
//import org.eclipse.swt.events.SelectionAdapter;
//import org.eclipse.swt.events.SelectionEvent;
//import org.eclipse.swt.graphics.Color;
//import org.eclipse.swt.graphics.Point;
//import org.eclipse.swt.graphics.Rectangle;
//import org.eclipse.swt.layout.GridData;
//import org.eclipse.swt.layout.GridLayout;
//import org.eclipse.swt.widgets.Button;
//import org.eclipse.swt.widgets.Canvas;
//import org.eclipse.swt.widgets.Composite;
//import org.eclipse.swt.widgets.Control;
//import org.eclipse.swt.widgets.Display;
//import org.eclipse.swt.widgets.Group;
//import org.eclipse.swt.widgets.Layout;
//import org.eclipse.swt.widgets.Shell;

//import com.x4m.coreui.util.ErrorReporter;
//import com.x4m.coreui.util.TypedSimpleInput;
//import com.x4m.graphics2d.GridTiledView;
//import com.x4m.mosaic.Mosaic;
//import com.x4m.mosaic.Mosaic.BorderStyle;
//import com.x4m.mosaic.Projector;
//import com.x4m.mosaic.Tile;
//import com.x4m.mosaic.TileAxis;
//import com.x4m.util.parse.ParseException;
//import com.x4m.util.parse.StandardParsers;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.ArrayMath;

public class GammaConstraintPicker extends JFrame {

//  static final double GMIN_DEFAULT = 1.5; // Minimum gamma default
//  static final double GMAX_DEFAULT = 2.5; // Maximum gamma default
//
//  public GammaConstraintPicker(
//      Sampling s, Set<Pick> minPicks, Set<Pick> maxPicks)
//  {
//    _s = s;
//    double[] v = _s.getValues();
//    _ny = _s.getCount();
//    _y = new float[_ny];
//    for (int i=0; i<_ny; i++) {
//      _y[i] = (float)v[i];
//    }
//
//    _nx = (int)(MAX_DEFAULT-MIN_DEFAULT)+1;
//    _x = new float[_nx];
//    for (int i=0; i<_nx; i++) {
//      _x[i] = (float)MIN_DEFAULT+i;
//    }
//
//    if (minPicks==null) {
//      _gMinFPick = new Pick(_s.getFirst(),GMIN_DEFAULT);
//      _gMinLPick = new Pick(_s.getLast(),GMIN_DEFAULT);
//      _minGPicks.add(_gMinFPick);
//      _minGPicks.add(_gMinLPick);
//    } else {
//      _minGPicks = new HashSet<>(minPicks);
//      Pick[] picks = minPicks.toArray(new Pick[0]);
//      int nPicks = picks.length;
//      Arrays.sort(picks);
//      _gMinFPick = picks[0];
//      _gMinLPick = picks[nPicks-1];
//    }
//    
//    if (maxPicks==null) {
//      _gMaxFPick = new Pick(_s.getFirst(),GMAX_DEFAULT);
//      _gMaxLPick = new Pick(_s.getLast(),GMAX_DEFAULT);
//      _maxGPicks.add(_gMaxFPick);
//      _maxGPicks.add(_gMaxLPick);
//    } else {
//      _maxGPicks = new HashSet<>(maxPicks);
//      Pick[] picks = maxPicks.toArray(new Pick[0]);
//      int nPicks = picks.length;
//      Arrays.sort(picks);
//      _gMaxFPick = picks[0];
//      _gMaxLPick = picks[nPicks-1];
//    }
//    
//    _minF = interpolate(_minGPicks);
//    _maxF = interpolate(_maxGPicks);
//    correctMin();
//    correctMax();
//  }
//
//  public void setTime(double time) {
//    System.out.println("time"+time);
//    _timeLine = ArrayMath.filldouble(time,_ny);
//  }
//
//  public Set<Pick> getMinPicks() {
//    return _minGPicks;
//  }
//
//  public Set<Pick> getMaxPicks() {
//    return _maxGPicks;
//  }
//
//  public double[] getMinGamma() {
//    return _minF;
//  }
//
//  public double[] getMaxGamma() {
//    return _maxF;
//  }
//
//  public static void main(String[] args) {
//    Sampling s = new Sampling(1501,2,0.0);
//    GammaConstraintPicker gcp = new GammaConstraintPicker(s,null,null);
//    gcp.setTime(2000.0);
//    gcp.setVisible(true);
//    gcp.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//  }
//
////  @Override
////  protected void okPressed() {
////    super.okPressed();
////  }
////
////  @Override
////  protected void cancelPressed() {
////    _minF = ArrayMath.filldouble(GMIN_DEFAULT,_ny);
////    _maxF = ArrayMath.filldouble(GMAX_DEFAULT,_ny);
////    super.cancelPressed();
////  }
//
//  @Override
//  protected Control createDialogArea(Composite parent) {
//    setTitle("Gamma Constraint Picker");
//    setMessage("Pick constraints");
//
//    parent.setLayout(new GridLayout());
//    GridData gd = new GridData(GridData.FILL,GridData.FILL,true,true);
//
//    _container = new Composite(parent,SWT.NONE);
//    _container.setLayout(new GridLayout(1,false));
//    gd = new GridData(SWT.FILL,SWT.FILL,true,true);
//    _container.setLayoutData(gd);
//
//    ErrorReporter er = new ErrorReporter() {
//      @Override
//      public void setError(String arg0) {}
//      @Override
//      public void clearError() {}
//    };
//
//    Composite top = new Composite(_container,SWT.NONE);
//    top.setLayout(new GridLayout(2,false));
//    gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//    top.setLayoutData(gd);
//
//    // First and last gamma parameters
//    {
//      Group g = new Group(top,SWT.NONE);
//      g.setLayout(new GridLayout(3,false));
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      g.setLayoutData(gd);
//      g.setText("First and last gamma values");
//
//      Composite c1 = new Composite(g,SWT.NONE);
//      c1.setLayout(new GridLayout(1,false));
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      c1.setLayoutData(gd);
//
//      TypedSimpleInput<Double> gMinF = new TypedSimpleInput<>(c1,
//          StandardParsers.DOUBLE,"Minimum gamma first","",
//          Double.toString(_gMinFPick._gamma),er);
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      gMinF.setLayoutData(gd);
//      
//      TypedSimpleInput<Double> gMinL = new TypedSimpleInput<>(c1,
//          StandardParsers.DOUBLE,"Minimum gamma last","",
//          Double.toString(_gMinLPick._gamma),er);
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      gMinL.setLayoutData(gd);
//
//      Composite c2 = new Composite(g,SWT.NONE);
//      c2.setLayout(new GridLayout(1,false));
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      c2.setLayoutData(gd);
//
//      TypedSimpleInput<Double> gMaxF = new TypedSimpleInput<>(c2,
//          StandardParsers.DOUBLE,"Maximum gamma first","",
//          Double.toString(_gMaxFPick._gamma),er);
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      gMaxF.setLayoutData(gd);
//      
//      TypedSimpleInput<Double> gMaxL = new TypedSimpleInput<>(c2,
//          StandardParsers.DOUBLE,"Maximum gamma last","",
//          Double.toString(_gMaxLPick._gamma),er);
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      gMaxL.setLayoutData(gd);
//      
//      addGammaListener(gMinF,_gMinFPick, true,   0);
//      addGammaListener(gMinL,_gMinLPick, true,_ny-1);
//      addGammaListener(gMaxF,_gMaxFPick,false,   0);
//      addGammaListener(gMaxL,_gMaxLPick,false,_ny-1);
//    }
//
//    // Pick mode buttons
//    {
//      Group g = new Group(top,SWT.NONE);
//      g.setLayout(new GridLayout(1,false));
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      g.setLayoutData(gd);
//      g.setText("Picking mode");
//
//      Composite c3 = new Composite(g,SWT.NONE);
//      c3.setLayout(new GridLayout(1,false));
//      gd = new GridData(SWT.FILL,SWT.CENTER,true,false);
//      c3.setLayoutData(gd);
//      Button pickMin = new Button(c3,SWT.RADIO);
//      gd = new GridData(SWT.LEFT,SWT.CENTER,false,false);
//      pickMin.setLayoutData(gd);
//      pickMin.setText("Pick minimum");
//      pickMin.setSelection(false);
//      pickMin.addSelectionListener(new SelectionAdapter() {
//        @Override
//        public void widgetSelected(SelectionEvent e) {
//          _pickMax = false;
//        }
//      });
//      Button pickMax = new Button(c3,SWT.RADIO);
//      gd = new GridData(SWT.LEFT,SWT.CENTER,false,false);
//      pickMax.setLayoutData(gd);
//      pickMax.setText("Pick maximum");
//      pickMax.setSelection(true);
//      pickMax.addSelectionListener(new SelectionAdapter() {
//        @Override
//        public void widgetSelected(SelectionEvent e) {
//          _pickMax = true;
//        }
//      });
//    }
//
//    Composite mosaicArea = new Composite(_container,SWT.NONE);
//    gd = new GridData(SWT.FILL,SWT.FILL,true,true);
//    mosaicArea.setLayoutData(gd);
//    
//    // Make a special layout manager for our top container to position
//    // the mosaic child at the full size of the top container.
//    mosaicArea.setLayout(new Layout() {
//      protected Point computeSize(
//          Composite composite, int wHint, int hHint, boolean flushCache) 
//      {
//        return _mosaic.computeSize(wHint,hHint,flushCache);
//      }
//      protected void layout(Composite composite, boolean flushCache) {
//        Point size = composite.getSize();
//        _mosaic.setBounds(new Rectangle(0,0,size.x,size.y));
//      }
//    });
//    Set<Mosaic.AxesPlacement> axesPlacement = EnumSet.of(
//        Mosaic.AxesPlacement.TOP,
//        Mosaic.AxesPlacement.LEFT);
//    _mosaic = new Mosaic(mosaicArea,1,1,axesPlacement,BorderStyle.FLAT);
//    _mosaic.setWidthMinimum(0,10);
//    TileAxis yAxis = _mosaic.getTileAxisLeft(0);
//    yAxis.setLabel("Time (ms)");
//    TileAxis xAxis = _mosaic.getTileAxisTop(0);
//    xAxis.setLabel("Interval gamma");
//    Tile tile = _mosaic.getTile(0,0);
//    tile.setBestHorizontalProjector(new Projector(MIN_DEFAULT,MAX_DEFAULT));
//    tile.setBestVerticalProjector(new Projector(_s.getFirst(),_s.getLast()));
//    _hp = tile.getHorizontalProjector();
//    _vp = tile.getVerticalProjector();
//
//    GridTiledView tv = new GridTiledView(xAxis,yAxis);
//    tile.addTiledView(tv);
//    _canvas = tile.getCanvas();
//    _canvas.addMouseListener(new MouseAdapter() {
//      @Override
//      public void mouseDown(MouseEvent e) {
//        int x = e.x;
//        int y = e.y;
//
//        if (e.button==1) { // Add a new pick
//          double[] worldXY = getWorldXY(x,y);
//          double time = worldXY[1];
//          double gamma = worldXY[0];
//          int i = _s.indexOfNearest(time);
//          if (_pickMax) {
//            gamma = gamma<_minF[i] ? _minF[i] : gamma;
//            Pick pick = new Pick(time,gamma);
//            _maxGPicks.add(pick);
//          } else {
//            gamma = gamma>_maxF[i] ? _maxF[i] : gamma;
//            Pick pick = new Pick(time,gamma);
//            _minGPicks.add(pick);
//          }
//        }
//
//        if (e.button==3) { // Try to remove the nearest pick
//          if (_pickMax) {
//            Pick np = getNearestPick(x,y,_maxGPicks.toArray(new Pick[0]));
//            if (np!=_gMaxFPick && np!=_gMaxLPick)
//              _maxGPicks.remove(np);
//          } else {
//            Pick np = getNearestPick(x,y,_minGPicks.toArray(new Pick[0]));
//            if (np!=_gMinFPick && np!=_gMinLPick)
//              _minGPicks.remove(np);
//          }
//        }
//        _canvas.redraw();
//      }
//    });
//
//    _canvas.addPaintListener(new PaintListener() {
//      @Override
//      public void paintControl(PaintEvent event) {
//        paint(event);
//      }
//    });
//    _canvas.redraw();
//    return parent;
//  }
//
//  @Override
//  protected boolean isResizable() {
//    return true;
//  }
//
//  @Override
//  protected Point getInitialSize() {
//    return new Point(800,900);
//  }
//
//  //////////////////////////////////////////////////////////////////////////
//  // Private
//  private final Sampling _s;
//  private final int _nx, _ny;
//  private final float[] _x, _y;
//  private final Pick _gMinFPick, _gMinLPick, _gMaxFPick, _gMaxLPick;
//  private Set<Pick> _maxGPicks = new HashSet<>();
//  private Set<Pick> _minGPicks = new HashSet<>();
//  private double[] _minF, _maxF, _timeLine=null;
//  private boolean _pickMax = true;
//  private Projector _hp, _vp;
//
//  private Composite _container;
//  private Mosaic _mosaic;
//  private Canvas _canvas;
//  
//  private static final int D  = 10;  // Diameter
//  private static final int HD = D/2; // Half diameter
//  private static final double MIN_DEFAULT =  1.0; // Minimum allowed value
//  private static final double MAX_DEFAULT = 10.0; // Maximum allowed value
//
//  private double[] interpolate(Set<Pick> pickSet) {
//    Pick[] picks = pickSet.toArray(new Pick[0]);
//    Arrays.sort(picks);
//    int nPicks = picks.length;
//    float[] t = new float[nPicks];
//    float[] g = new float[nPicks];
//    for (int i=0; i<nPicks; i++) {
//      Pick pick = picks[i];
//      g[i] = (float)pick._gamma;
//      t[i] = (float)pick._time;
//    }
//    CubicInterpolator ci = new CubicInterpolator(t,g);
//    float[] xi = ci.interpolate(_y);
//    double[] xid = new double[_ny];
//    for (int i=0; i<_ny; i++)
//      xid[i] = xi[i]; 
//    return xid;
//  }
//
//  private void correctMin() {
//    for (int i=0; i<_ny; i++) {
//      if (_minF[i]>_maxF[i])
//        _minF[i] = _maxF[i];
//    }
//  }
//
//  private void correctMax() {
//    for (int i=0; i<_ny; i++) {
//      if (_maxF[i]<_minF[i])
//        _maxF[i] = _minF[i];
//    }
//  }
//
//  private void addGammaListener(
//      final TypedSimpleInput<Double> tsi, final Pick pick, final boolean min, 
//      final int index) 
//  {
//    tsi.getTextControl().addFocusListener(new FocusAdapter() {
//      @Override
//      public void focusLost(FocusEvent event) {
//        try {
//          double gamma = tsi.getInputValue().doubleValue();
//          if (min) {
//            gamma = gamma>_maxF[index] ? _maxF[index] : 
//                    gamma<MIN_DEFAULT  ? MIN_DEFAULT  : gamma;   
//            tsi.getTextControl().setText(Double.toString(gamma));
//          } else {
//            gamma = gamma<_minF[index] ? _minF[index] : 
//                    gamma>MAX_DEFAULT  ? MAX_DEFAULT  : gamma;   
//            tsi.getTextControl().setText(Double.toString(gamma));
//          }
//          pick._gamma = gamma;
//          _canvas.redraw();
//        } catch (ParseException e) {
//          
//        }
//        super.focusLost(event);
//      }
//    });
//  }
//
//  private double[] getWorldXY(int x, int y) {
//    Point p = _canvas.getSize();
//    int w = p.x;
//    int h = p.y;
//    double wx = _hp.v((double)x/w);
//    double wy = _vp.v((double)y/h);
//    return new double[]{wx,wy};
//  }
//
//  private Pick getNearestPick(int x, int y, Pick[] picks) {
//    Point point = _canvas.getSize();
//    int w = point.x;
//    int h = point.y;
//    float minD = Float.MAX_VALUE;
//    Pick nearestPick = null;
//    for (Pick p : picks) {
//      int xp = (int)(_hp.u(p._gamma)*w);
//      int yp = (int)(_vp.u(p._time)*h);
//      float a = x-xp;
//      float b = y-yp;
//      float d = ArrayMath.sqrt((a*a)+(b*b));
//      if (d<minD && d<10.0f) {
//        minD = d;
//        nearestPick = p;
//      }
//    }
//    return nearestPick;
//  }
//
//  private void paint(PaintEvent event) {
//    paintPoints(event);
//    event.gc.setLineWidth(3);
//    event.gc.setInterpolation(SWT.HIGH);
//    event.gc.setForeground(new Color(event.display,255,0,0));
//    _minF = interpolate(_minGPicks);
//    correctMin();
//    paintVerticalLine(event,_minF);
//    event.gc.setForeground(new Color(event.display,0,0,255));
//    _maxF = interpolate(_maxGPicks);
//    correctMax();
//    paintVerticalLine(event,_maxF);
//    if (_timeLine!=null) {
//      event.gc.setForeground(new Color(event.display,0,255,0));
//      paintHorizontalLine(event,_timeLine);
//    }
//  }
//
//  private void paintPoints(PaintEvent event) {
//    Point p = _canvas.getSize();
//    int w = p.x;
//    int h = p.y;
//    Pick[] minPicks = _minGPicks.toArray(new Pick[0]);
//    Pick[] maxPicks = _maxGPicks.toArray(new Pick[0]);
//    event.gc.setBackground(new Color(event.display,255,0,0));
//    for (Pick pick : minPicks) {
//      int x = (int)(_hp.u(pick._gamma)*w);
//      int y = (int)(_vp.u(pick._time)*h);
//      paintPoint(event,x,y);
//    }
//    event.gc.setBackground(new Color(event.display,0,0,255));
//    for (Pick pick : maxPicks) {
//      int x = (int)(_hp.u(pick._gamma)*w);
//      int y = (int)(_vp.u(pick._time)*h);
//      paintPoint(event,x,y);
//    }
//  }
//
//  private void paintPoint(PaintEvent event, int x, int y) {
//    event.gc.fillOval(x-HD,y-HD,D,D);
//  }
//
//  private void paintVerticalLine(PaintEvent event, double[] d) {
//    Point p = _canvas.getSize();
//    int w = p.x;
//    int h = p.y;
//    int[] pointArray = new int[_ny*2];
//    for (int i=0; i<_ny; i++) {
//      int j = i*2;
//      pointArray[j  ] = (int)(_hp.u( d[i])*w);
//      pointArray[j+1] = (int)(_vp.u(_y[i])*h);
//    }
//    event.gc.drawPolyline(pointArray);
//  }
//
//  private void paintHorizontalLine(PaintEvent event, double[] d) {
//    Point p = _canvas.getSize();
//    int w = p.x;
//    int h = p.y;
//    int[] pointArray = new int[_nx*2];
//    for (int i=0; i<_nx; i++) {
//      int j = i*2;
//      pointArray[j  ] = (int)(_hp.u(_x[i])*w);
//      pointArray[j+1] = (int)(_vp.u( d[i])*h);
//    }
//    event.gc.drawPolyline(pointArray);
//  }
//
//  public static class Pick implements Comparable<Pick> {
//
//    public Pick(double time, double gamma) {
//      _time = time;
//      _gamma = gamma;
//    }
//
//    @Override
//    public int compareTo(Pick p) {
//      if (this._time<p._time) return -1;
//      if (this._time>p._time) return  1;
//      return 0;
//    }
//
//    @Override
//    public boolean equals(Object obj) {
//      if (obj instanceof Pick) {
//        Pick p = (Pick)obj;
//        return p._time==this._time;
//      } else {
//        return false;
//      }
//    }
//
//    @Override
//    public int hashCode() {
//      int result = 17;
//      long l = Double.doubleToLongBits(_time);
//      result = 31*result+(int)(l^(l>>>32));
//      return result;
//    }
//
//    private double _time;
//    private double _gamma;
//  }

}
