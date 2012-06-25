package segy;

import java.util.logging.Logger;

import edu.mines.jtk.dsp.Fft;
import edu.mines.jtk.dsp.Real1;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.la.DMatrix;
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
	
	public static Sampling getSampling(float[] f) {
		return new Fft(f).getFrequencySampling1();
	}
	
	/**
	 * Computes the ratio of a weighted average of high
	 * frequencies to a weighted average of low frequencies
	 * for the given array of amplitudes.
	 * @param a the array of amplitudes
	 * @return the ratio of high to low frequencies
	 */
	public static float computeLinearRatio(float[] a, float[] highW, float[] lowW) {
		int len = a.length;
		float hSum  = 0;
		float lSum  = 0;
		float hwSum = 0;
		float lwSum = 0;
		for (int i = 0; i < len; i++) {
			float hw = (float)i / (float)(len-1);
			hwSum += hw;
			hSum  += hw*a[i];
			if (highW != null) {
				highW[i] = hw;
			}
			
			float lw = (float)(len-1-i) / (float)(len-1);
			lwSum += lw;
			lSum  += lw*a[i];
			if (lowW != null) {
				lowW[i] = lw;
			}
		}
		float h = hSum/hwSum;
		float l = lSum/lwSum;
		float r = h/l;
		LOG.info("h="+h+", l="+l+", r="+r);
		return r;
	}
	
	public static float computeSCRatio(float[] a, float[] highW, float[] lowW) {
		int len = a.length;
		float hSum  = 0;
		float lSum  = 0;
		float hwSum = 0;
		float lwSum = 0;
		for (int i = 0; i < len; i++) {
			float radians = (float)Math.toRadians((float)i / (float)(len-1) * 90f);
			float hw = (float)Math.sin(radians);
			hwSum += hw;
			hSum  += hw*a[i];
			if (highW != null) {
				highW[i] = hw;
			}
			
			float lw = (float)Math.cos(radians);
			lwSum += lw;
			lSum  += lw*a[i];
			if (lowW != null) {
				lowW[i] = lw;
			}
		}
		float h = hSum/hwSum;
		float l = lSum/lwSum;
		float r = h/l;
		LOG.info("h="+h+", l="+l+", r="+r);
		return r;
	}
	
	public static float computeRatio(float[] a, float p, float n,
			float[] highW, float[] lowW) {
		int len = a.length;
		float hSum  = 0;
		float lSum  = 0;
		float wSum = 0;
		for (int i = 0; i < len; i++) {
			float x = (float)i / (float)(len-1);
			float w = getWeight(x, p, n);
			wSum += w;
			hSum  += w*a[i];
			lSum  += w*a[(len-1)-i];
			if (highW != null) {
				highW[i] = w;
			}
			if (lowW != null) {
				lowW[len-1-i] = w;
			}
		}
		float h = hSum/wSum;
		float l = lSum/wSum;
		float r = h/l;
		LOG.info("h="+h+", l="+l+", r="+r);
		return r;
	}
	
	public static Real1[] getWeights(Sampling s, float p, float n) {
		int nk = s.getCount();
		float[] hwArray = new float[nk];
		float[] lwArray = new float[nk];
		for (int i = 0; i < nk; i++) {
			float x = (float)i / (float)(nk-1);
			float w = getWeight(x, p, n);
			hwArray[i] = w;
			lwArray[(nk-1)-i] = w;
		}
		Real1[] weights = {new Real1(s, hwArray), new Real1(s, lwArray)};
		return weights;
	}
	
	private static float getWeight(float x, float slope, float power) {
		float w;
		if (5*x < slope) {
			w = (float)((3*Math.pow(5*x/slope, 8) - 2*Math.pow(5*x/slope, 12))
					* Math.pow(x, power));
		} else {
			w = (float)Math.pow(x, power);
		}
		return w;
	}
}