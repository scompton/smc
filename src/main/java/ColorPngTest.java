import java.awt.*;
import javax.swing.*;

import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import viewer.*;

public class ColorPngTest {

  public static void main(String[] args) {
    final float[] f = rampfloat(0.0f,1.0f,100);
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        tileBackground(f);
        //panelBackground(f);
        //frameBackground(f);
      }
    });
  }

  private static void tileBackground(float[] f) {
    Viewer2D v = new Viewer2D();
    v.addPixels(new float[100][100],"");
    v.addPoints(f,"");
    v.show();
    //SimplePlot sp = new SimplePlot();
    //sp.setTitle("tile background color");
    //sp.addPoints(f);
    //sp.addColorBar();
    //sp.getPlotPanel().getTile(0,0).setBackground(Color.LIGHT_GRAY); 
    //sp.paintToPng(360.0,3.0,"junk_tileBackground.png");
  }

  private static void panelBackground(float[] f) {
    SimplePlot sp = new SimplePlot();
    sp.setTitle("panel background color");
    sp.addPoints(f);
    sp.addColorBar();
    sp.getPlotPanel().setBackground(Color.LIGHT_GRAY); 
    sp.paintToPng(360.0,3.0,"junk_panelBackground.png");
  }

  private static void frameBackground(float[] f) {
    SimplePlot sp = new SimplePlot();
    sp.setTitle("frame background color");
    sp.addPoints(f);
    sp.setBackground(Color.LIGHT_GRAY); 
    sp.paintToPng(360.0,3.0,"junk_frameBackground.png");
  }

}

