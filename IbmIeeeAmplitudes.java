import java.util.logging.Logger;

import edu.mines.jtk.dsp.Fft;
import edu.mines.jtk.dsp.Real1;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.ArrayMath;
import edu.mines.jtk.util.Cdouble;


public class IbmIeeeAmplitudes {
	
	private static final Logger LOG = Logger.getLogger(IbmIeeeAmplitudes.class.getName());
	
	/**
	 * 
	 * @param f
	 * @param db
	 * @return
	 */
	public static Real1 computeAmplitudes(float[] f, boolean db) {
		Fft fft = new Fft(f);
		Sampling sk = fft.getFrequencySampling1();
		int nk = sk.getCount();
		float[] g = fft.applyForward(f);
		float[] mag = new float[nk];
		for (int kk=0,kr=0,ki=kr+1; kk<nk; kk++,kr+=2,ki+=2) {
			Cdouble c = new Cdouble(g[kr], g[ki]);
			mag[kk] = (float)c.abs();
		}
		if (db) {
			mag = ArrayMath.log10(mag);
			mag = ArrayMath.mul(10.0f, mag);
		}
		Real1 amplitudes = new Real1(sk, mag);
		return amplitudes;
	}
	
	/**
	 * Computes the ratio of a weighted average of high
	 * frequencies to a weighted average of low frequencies
	 * for the given array of amplitudes.
	 * @param a the array of amplitudes
	 * @return the ratio of high to low frequencies
	 */
	public static float computeRatio(float[] a) {
		int len = a.length;
		float hSum  = 0;
		float lSum  = 0;
		float hwSum = 0;
		float lwSum = 0;
		for (int i = 0; i < len; i++) {
			float hw = (float)i / (float)(len-1);
			hwSum += hw;
			hSum  += hw*a[i];
			
			float lw = (float)(len-1-i) / (float)(len-1);
			lwSum += lw;
			lSum  += lw*a[i];
		}
		float h = hSum/hwSum;
		float l = lSum/lwSum;
		float r = h/l;
		LOG.info("h="+h+", l="+l+", r="+r);
		return r;
	}
	
}
