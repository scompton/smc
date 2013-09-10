package viewer;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotPanelPixels3;

public class ClipFrame extends JPanel implements PropertyChangeListener {
  
  /**
   * Constructs a JFrame with options to change the minimum and
   * maximum clips. The frame is constructed here, but is not 
   * made visible by default. Use the {@link #getFrame()} method 
   * to make the frame visible when required. 
   * </p>
   * Clips are changed for all PixelsViews in the given array. For 
   * instance, if the plot panel contains multiple PixelsViews of
   * the same data, which is the case for a 
   * {@link PlotPanelPixels3}, this array should contain each of
   * those PixelsViews in order to update them simultaneously. 
   * @param pv the array of PixelsViews. 
   */
  public ClipFrame(PixelsView[] pv) {
    _pv = pv;
    NumberFormat clipFormat = NumberFormat.getNumberInstance();
    _min = new JFormattedTextField(clipFormat);
    _min.setColumns(10);
    _min.setValue(_minValue);
    _min.addPropertyChangeListener("value",this);
    _max = new JFormattedTextField(clipFormat);
    _max.setColumns(10);
    _max.setValue(_maxValue);
    _max.addPropertyChangeListener("value",this);
    JLabel minLabel = new JLabel("min:");
    JLabel maxLabel = new JLabel("max:");
    add(minLabel);
    add(_min);
    add(maxLabel);
    add(_max);
    
    _frame = new JFrame("Enter Clips");
    _frame.add(this);
    _frame.pack();
    _frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
  }

  /**
   * Returns the frame.
   * @return the frame.
   */
  public JFrame getFrame() {
    return _frame;
  }
  
  /**
   * Set the minimum and maximum values that the dialog should
   * display.
   * @param minValue
   * @param maxValue
   */
  public void setValues(float minValue, float maxValue) {
    _minValue = minValue;
    _maxValue = maxValue;
    _min.setValue(_minValue);
    _max.setValue(_maxValue);
  }

  @Override
  public void propertyChange(PropertyChangeEvent e) {
    Object source = e.getSource();
    if (source == _min) {
      float clipMin = ((Number)_min.getValue()).floatValue();
      for (PixelsView pv : _pv)
        pv.setClips(clipMin,_maxValue);
      _minValue = clipMin;
    }
    if (source == _max) {
      float clipMax = ((Number)_max.getValue()).floatValue();
      for (PixelsView pv : _pv)
        pv.setClips(_minValue,clipMax);
      _maxValue = clipMax;
    }
  }
  
  public static void main(String[] args) {
    new ClipFrame(null);
  }
  
  private JFrame _frame;
  private float _minValue = -1.0f;
  private float _maxValue =  1.0f;
  private PixelsView[] _pv;
  private JFormattedTextField _min;
  private JFormattedTextField _max;
  private static final long serialVersionUID = 1L;
  
}
