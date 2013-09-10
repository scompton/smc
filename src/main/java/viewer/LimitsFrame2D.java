package viewer;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import edu.mines.jtk.mosaic.PlotPanelPixels3;

public class LimitsFrame2D extends JPanel implements PropertyChangeListener {
  
  /**
   * Constructs a JFrame with options to set the limits for the 
   * horizontal axis and vertical axis of a {@link PlotPanel}. The 
   * frame is constructed here, but is not made visible by default. 
   * Use the {@link #getFrame()} method to access the frame and make 
   * the it visible when required. 
   * @param v the Viewer3P instance that contains the 
   *  {@link PlotPanelPixels3} instance.
   */
  public LimitsFrame2D(Viewer2D v) {
    super(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    c.insets = new Insets(0,10,0,0);
    
    _v = v;
    _hMinVal = v.getHMin();
    _hMaxVal = v.getHMax();
    _vMinVal = v.getVMin();
    _vMaxVal = v.getVMax();
    
    NumberFormat limitsFormat = NumberFormat.getNumberInstance();
    limitsFormat.setGroupingUsed(false);
    
    JLabel k1Label = new JLabel("Horizontal Limits:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 0;
    add(k1Label,c);
    _hMin = new JFormattedTextField(limitsFormat);
    _hMin.setColumns(5);
    _hMin.setValue(_hMinVal);
    _hMin.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 1;
    c.gridy = 0;
    add(_hMin,c);
    _hMax = new JFormattedTextField(limitsFormat);
    _hMax.setColumns(5);
    _hMax.setValue(_hMaxVal);
    _hMax.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 2;
    c.gridy = 0;
    add(_hMax,c);
        
    JLabel k2Label = new JLabel("Vertical Limits:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 1;
    add(k2Label,c);
    _vMin = new JFormattedTextField(limitsFormat);
    _vMin.setColumns(5);
    _vMin.setValue(_vMinVal);
    _vMin.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;    
    c.gridx = 1;
    c.gridy = 1;
    add(_vMin,c);
    _vMax = new JFormattedTextField(limitsFormat);
    _vMax.setColumns(5);
    _vMax.setValue(_vMaxVal);
    _vMax.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 2;
    c.gridy = 1;
    add(_vMax,c);
    
    _frame = new JFrame("Enter Limits");
    _frame.add(this);
    _frame.setSize(300,110);
    _frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
  }
  
  /**
   * Returns the frame.
   * @return the frame.
   */
  public JFrame getFrame() {
    return _frame;
  }

  @Override
  public void propertyChange(PropertyChangeEvent e) {
    Object source = e.getSource();
    if (source==_hMin) {
      double hMinVal = ((Number)_hMin.getValue()).doubleValue();
      _v.setHLimits(hMinVal,_hMaxVal);
      _hMinVal = hMinVal;
    }
    if (source==_hMax) {
      double hMaxVal = ((Number)_hMax.getValue()).doubleValue();
      _v.setHLimits(_hMinVal,hMaxVal);
      _hMaxVal = hMaxVal;
    }
    if (source==_vMin) {
      double vMinVal = ((Number)_vMin.getValue()).doubleValue();
      _v.setVLimits(vMinVal,_vMaxVal);
      _vMinVal = vMinVal;
    }
    if (source==_vMax) {
      double vMaxVal = ((Number)_vMax.getValue()).doubleValue();
      _v.setVLimits(_vMinVal,vMaxVal);
      _vMaxVal = vMaxVal;
    }
  }
  
  public static void main(String[] args) {
//    LimitsFrame2D lf2d = 
//        new LimitsFrame2D(new Viewer2DNew(ArrayMath.zerofloat(1999,500,14)));
//    lf2d.getFrame().setVisible(true);
  }

  private final JFrame _frame;
  private final Viewer2D _v;
  private JFormattedTextField _hMin; // input for limits 1 min.
  private JFormattedTextField _hMax; // input for limits 1 max.
  private JFormattedTextField _vMin; // input for limits 2 min.
  private JFormattedTextField _vMax; // input for limits 2 max.
  private double _hMinVal;
  private double _hMaxVal;
  private double _vMinVal;
  private double _vMaxVal;
  private static final long serialVersionUID = 1L;
  
}
