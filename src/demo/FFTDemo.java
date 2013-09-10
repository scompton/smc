package demo;

import java.util.Random;

import javax.swing.SwingUtilities;

import edu.mines.jtk.dsp.Conv;
import edu.mines.jtk.dsp.Fft;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;

import static edu.mines.jtk.util.ArrayMath.*;

public class FFTDemo {
  
  public FFTDemo(final Sampling s, float fPeakHz) {
    // Convert to Cycles per Sample, and make a synthetic seismic trace.
    float fPeakCS = fPeakHz*(float)s.getDelta();
    final float[] f = makeSyntheticTrace(s.getCount(),fPeakCS);
    
    // Forward FFT and compute amplitude spectrum.
    final Fft fft = new Fft(s);
    fft.setPadding(2*s.getCount());
    final float[] g = fft.applyForward(f);
    final float[] amp = cabs(g);
    final float[] ampDb = mul(20.0f,log10(cabs(g)));
    
    // Inverse FFT
    final float[] h = fft.applyInverse(g);
    
    // Plot the results
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        SimplePlot spF = new SimplePlot();
        spF.addPoints(s,f);
        spF.setTitle("f");
        spF.setSize(WIDTH,HEIGHT);
        spF.setVisible(true);
        
        Sampling s1 = fft.getFrequencySampling1();
        SimplePlot spAmp = new SimplePlot();
        spAmp.addPoints(s1,amp);
        spAmp.setTitle("Amplitude");
        spAmp.setSize(WIDTH,HEIGHT);
        spAmp.setVisible(true);
        
        SimplePlot spAmpDb = new SimplePlot();
        spAmpDb.addPoints(s1,ampDb);
        spAmpDb.setTitle("Amplitude dB");
        spAmpDb.setSize(WIDTH,HEIGHT);
        spAmpDb.setVisible(true);
        
        SimplePlot spH = new SimplePlot();
        spH.addPoints(s,h);
        spH.setTitle("h");
        spH.setSize(WIDTH,HEIGHT);
        spH.setVisible(true);
      }
    });
  }

  public static void main(String[] args) {
    // Make sampling with 1001 samples, and a sampling rate of 2ms.
    Sampling s = new Sampling(1001,0.002,0.0);
    new FFTDemo(s,50.0f);
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Private 
  
  private static final int WIDTH = 1000;
  private static final int HEIGHT = 800;
  
  private static float[] makeSyntheticTrace(int n, float fpeak) {
    Random r = new Random(100);
    float[] f = pow(mul(sub(randfloat(r,n),0.5f),2.0f),7.0f);
    int ih = (int)(3.0f/fpeak);
    int nh = 1+2*ih;
    float[] h = zerofloat(nh);
    for (int jh=0; jh<nh; jh++) {
      float x = FLT_PI*fpeak*(jh-ih);
      h[jh] = (1.0f-2.0f*x*x)*exp(-x*x);
    }
    float[] g = zerofloat(n);
    Conv.conv(nh,-ih,h,n,0,f,n,0,g);
    return g;
  }

}
