package gbc;

import java.io.*;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A gridded horizon from Geokinetics Bradford County.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.28
 */
public class Horizon {

  public int ns; // number of samples
  public float[] x1,x2,x3; // sample coordinates
  public int nt; // number of triangles
  public int[] ia,ib,ic; // triangle indices

  /**
   * Reads a horizon from a binary file with the specified name.
   * @param fileName the file name.
   * @return the horizon.
   */
  public static Horizon readBinary(String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int ns = ais.readInt();
      float[] x1 = new float[ns]; ais.readFloats(x1);
      float[] x2 = new float[ns]; ais.readFloats(x2);
      float[] x3 = new float[ns]; ais.readFloats(x3);
      int nt = ais.readInt();
      int[] ia = new int[nt]; ais.readInts(ia);
      int[] ib = new int[nt]; ais.readInts(ib);
      int[] ic = new int[nt]; ais.readInts(ic);
      ais.close();
      return new Horizon(ns,x1,x2,x3,nt,ia,ib,ic);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a horizon from a text file with the specified name.
   * The file format is that used by Transform Software.
   * The returned horizon includes only those points that lie within
   * the specified sampling bounds.
   * @param fileName the file name.
   * @return the horizon.
   */
  public static Horizon readText(String fileName) {
    double nullValue = -999.99;
    int m = 500, n = 500; // (m,n) = max (inline,xline) dimensions
    int[][] index = fillint(-1,n,m);
    int ns = 0;
    FloatList x1List = new FloatList();
    FloatList x2List = new FloatList();
    FloatList x3List = new FloatList();
    try {
      BufferedReader br = new BufferedReader(new FileReader(fileName));
      String name = null;
      StringBuilder sb = null;
      String line = null;
      while ((line=br.readLine())!=null) {
        if (line.startsWith("#"))
          continue;
        String[] fields = line.split("\\s+");
        int i = Integer.parseInt(fields[0])-MIN_IL+1;
        int j = Integer.parseInt(fields[1])-MIN_XL+1;
        double xe = Double.parseDouble(fields[2]);
        double yn = Double.parseDouble(fields[3]);
        double zd = Double.parseDouble(fields[4]);
        assert index[i][j]<0:"no duplicate samples";
        if (xe==nullValue || yn==nullValue || zd==nullValue)
          continue;
	Coordinates.Map map = new Coordinates.Map(xe,yn);
	Coordinates.Nrm nrm = new Coordinates.Nrm(map);
	x1List.add((float)zd/1000.0f);
	x2List.add((float)nrm.x2);
	x3List.add((float)nrm.x3);
        index[i][j] = ns++;
      }
      br.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    float[] x1 = x1List.trim();
    float[] x2 = x2List.trim();
    float[] x3 = x3List.trim();
    IntList iaList = new IntList();
    IntList ibList = new IntList();
    IntList icList = new IntList();
    for (int i=1; i<m; ++i) {
      for (int j=1; j<n; ++j) {
        int imm = index[i-1][j-1];
        int imp = index[i-1][j  ];
        int ipm = index[i  ][j-1];
        int ipp = index[i  ][j  ];
        if (imm>=0 && imp>=0 && ipm>=0 && ipp>=0) {
          iaList.add(imp); ibList.add(imm); icList.add(ipm);
          iaList.add(imp); ibList.add(ipm); icList.add(ipp);
        } else if (  imp>=0 &&        imm>=0 &&        ipm>=0) {
          iaList.add(imp); ibList.add(imm); icList.add(ipm);
        } else if (  ipp>=0 &&        imp>=0 &&        imm>=0) {
          iaList.add(ipp); ibList.add(imp); icList.add(imm);
        } else if (  ipm>=0 &&        ipp>=0 &&        imp>=0) {
          iaList.add(ipm); ibList.add(ipp); icList.add(imp);
        } else if (  imm>=0 &&        ipm>=0 &&        ipp>=0) {
          iaList.add(imm); ibList.add(ipm); icList.add(ipp);
        }
      }
    }
    int[] ia = iaList.trim();
    int[] ib = ibList.trim();
    int[] ic = icList.trim();
    int nt = ia.length;
    return new Horizon(ns,x1,x2,x3,nt,ia,ib,ic);
  }

  /**
   * Clips this horizon to bounds implied by the specified samplings.
   * @param s2 sampling of CSM coordinate x2
   * @param s3 sampling of CSM coordinate x3
   */
  public void clip(Sampling s2, Sampling s3) {
    double f2 = s2.getFirst();
    double f3 = s3.getFirst();
    double l2 = s2.getLast();
    double l3 = s3.getLast();
    int[] index = fillint(-1,ns);
    FloatList x1List = new FloatList();
    FloatList x2List = new FloatList();
    FloatList x3List = new FloatList();
    int ms = 0;
    for (int is=0; is<ns; ++is) {
      if (f2<=x2[is] && x2[is]<=l2 &&
          f3<=x3[is] && x3[is]<=l3) {
        index[is] = ms++;
        x1List.add(x1[is]);
        x2List.add(x2[is]);
        x3List.add(x3[is]);
      }
    }
    IntList iaList = new IntList();
    IntList ibList = new IntList();
    IntList icList = new IntList();
    int mt = 0;
    for (int it=0; it<nt; ++it) {
      int ja = index[ia[it]];
      int jb = index[ib[it]];
      int jc = index[ic[it]];
      if (ja>=0 && jb>=0 && jc>=0) {
        ++mt;
        iaList.add(ja);
        ibList.add(jb);
        icList.add(jc);
      }
    }
    ns = ms;
    x1 = x1List.trim();
    x2 = x2List.trim();
    x3 = x3List.trim();
    nt = mt;
    ia = iaList.trim();
    ib = ibList.trim();
    ic = icList.trim();
  }

  /**
   * Writes this horizon to a binary file with the specified name.
   * @param fileName the file name.
   */
  public void writeBinary(String fileName) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeInt(ns);
      aos.writeFloats(x1);
      aos.writeFloats(x2);
      aos.writeFloats(x3);
      aos.writeInt(nt);
      aos.writeInts(ia);
      aos.writeInts(ib);
      aos.writeInts(ic);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Gets packed triplets (ia,ib,ic) of triangle indices.
   * @return array[3*nt] of triangle indices.
   */
  public int[] getIABC() {
    int[] iabc = new int[3*nt];
    for (int it=0,m=0; it<nt; ++it) {
      iabc[m++] = ia[it];
      iabc[m++] = ib[it];
      iabc[m++] = ic[it];
    }
    return iabc;
  }

  /**
   * Gets packed triplets (x3,x2,x1) of sample coordinates.
   * @return array[3*ns] of sample coordinates.
   */
  public float[] getX321() {
    float[] x321 = new float[3*ns];
    for (int is=0,m=0; is<ns; ++is) {
      x321[m++] = x3[is];
      x321[m++] = x2[is];
      x321[m++] = x1[is];
    }
    return x321;
  }

  /**
   * Gets packed triplets (x3,x2,x1) of sample coordinates for one x2
   * index. x3 and x2 coordinates are returned from given samplings.
   * This method checks that the horizon has exactly n2*n3 samples.
   * It is assumed that the horizon coordinates are nearly equivalent
   * to the Sampling coordinates, but this is not checked.
   * @return array[3*n2] of sample coordinates.
   */
  public float[] getI2X321(Sampling s2, Sampling s3, int i2) {
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    Check.argument(n2*n3==this.ns,"Number of input samplings (n2*n3) is "+
      "inconsistent with number of horizon samples");
    float x2 = (float)s2.getValue(i2);
    double d3 = s3.getDelta();
    double f3 = s3.getFirst();
    double x3 = f3;
    float[] x321 = new float[3*n3];
    for (int i3=0,ih=0,is=i2; i3<n3; i3++,is+=n2,x3=f3+i3*d3) {
      x321[ih++] = (float)x3;
      x321[ih++] = x2;
      x321[ih++] = this.x1[is];
    }
    return x321;
  }

  /**
   * Gets packed triplets (x3,x2,x1) of sample coordinates for one x3
   * index. x3 and x2 coordinates are returned from given samplings.
   * This method checks that the horizon has exactly n2*n3 samples.
   * It is assumed that the horizon coordinates are nearly equivalent
   * to the Sampling coordinates, but this is not checked.
   * @return array[3*n2] of sample coordinates.
   */
  public float[] getI3X321(Sampling s2, Sampling s3, int i3) {
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    Check.argument(n2*n3==this.ns,"Number of input samplings (n2*n3) is "+
      "inconsistent with number of horizon samples");
    float x3 = (float)s3.getValue(i3);
    double d2 = s2.getDelta();
    double f2 = s2.getFirst();
    double x2 = f2;
    float[] x321 = new float[3*n2];
    for (int i2=0,ih=0,is=i3*n2; i2<n2; i2++,x2=f2+i2*d2) {
      x321[ih++] = x3;
      x321[ih++] = (float)x2;
      x321[ih++] = this.x1[is++];
    }
    return x321;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  private Horizon(
    int ns, float[] x1, float[] x2, float[] x3,
    int nt, int[] ia, int[] ib, int[] ic) 
  {
    this.ns = ns; this.x1 = x1; this.x2 = x2; this.x3 = x3;
    this.nt = nt; this.ia = ia; this.ib = ib; this.ic = ic;
  }

  private static final int MIN_IL = 5947;
  private static final int MIN_XL = 5669;
}