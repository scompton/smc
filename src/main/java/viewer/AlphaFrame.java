package viewer;

import java.awt.image.IndexColorModel;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PixelsView;

public class AlphaFrame extends JPanel implements ChangeListener {
  
  public AlphaFrame(PixelsView[] pv) {
    _pv = pv;
    DefaultBoundedRangeModel brm = new DefaultBoundedRangeModel(100,0,0,100);
    _slider = new JSlider(brm);
    _slider.setMajorTickSpacing(25);
    _slider.setMinorTickSpacing(5);
    _slider.setPaintLabels(true);
    _slider.setPaintTicks(true);
    _slider.addChangeListener(this);
    add(_slider);
    
    _frame = new JFrame("Transparency");
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
  
  public void setAlpha(float alpha) {
    _slider.setValue((int)(alpha*100));
  }
  
  @Override
  public void stateChanged(ChangeEvent e) {
    JSlider source = (JSlider)e.getSource();
    int perc = source.getValue();
    IndexColorModel icm = _pv[0].getColorModel();
    for (PixelsView pv : _pv)
      pv.setColorModel(ColorMap.setAlpha(icm,perc/100.0));
  }
  
  private JFrame _frame;
  private JSlider _slider;
  private PixelsView[] _pv;
  private static final long serialVersionUID = 1L;
  
}
