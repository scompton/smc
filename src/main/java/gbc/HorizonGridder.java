package gbc;

import java.util.List;
import java.util.ArrayList;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;

public class HorizonGridder {

  public HorizonGridder(Sampling s2, Sampling s3, Horizon horizon) {
    _s2 = s2;
    _s3 = s3;
    _hns = horizon.ns;
    print("hns="+_hns+", hnt="+horizon.nt);
    _hx1 = horizon.x1;
    _hx2 = horizon.x2;
    _hx3 = horizon.x3;
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    _gh = new float[n3][n2];
    // getRow();

    // // Grid all the rows of horizon values (s2 coordinates).
    // double d2 = _s2.getDelta();
    // double f2 = _s2.getFirst();
    // double v2 = f2;
    // float[] hrow; // horizon row
    // float[] growX1 = new float[n2]; // gridded row
    // float[] growX3 = new float[n2];
    // List<float[][]> growList = new ArrayList<>();
    // int interpCount = 0;
    // int ongridCount = 0;
    // while ((hrow = getRow())!=null) {
    //   int ibs = 0; // binary search index
    //   int nr = hrow.length;
    //   for (int i2=0; i2<n2; i2++,v2=f2+i2*d2) {
    //     ibs = binarySearch(hrow,(float)v2,ibs);
    //     float x1,x3;
    //     if (ibs<0) { // linearly interpolate
    //       int is0 = -ibs-2;
    //       int is1 =  is0+1;
    //       if (is0<0) {
    //         int is = getSampleIndex(nr,0);
    //         x1 = _hx1[is]; x3 = _hx3[is];
    //         ibs = 0;
    //       } else if (is1>=nr) {
    //         int is = getSampleIndex(nr,nr-1);
    //         x1 = _hx1[is]; x3 = _hx3[is];
    //         ibs = nr-1;
    //       } else {
    //         is0 = getSampleIndex(nr,is0);
    //         is1 = getSampleIndex(nr,is1);
    //         float x20 = _hx2[is0];
    //         float x21 = _hx2[is1];
    //         float x10 = _hx1[is0];
    //         float x11 = _hx1[is1];
    //         float x30 = _hx3[is0];
    //         float x31 = _hx3[is1];
    //         float n = (float)v2-x20;
    //         float d = x21-x20;
    //         float s = n/d;
    //         x1 = x10+(x11-x10)*s;
    //         // print("x1="+x1);
    //         if (Float.isNaN(x1)) dump(hrow);
    //         assert !Float.isNaN(x1) : "hns="+_hns+", ibs="+ibs+", nr="+nr+
    //           ", is0="+is0+", is1="+is1+", x20="+x20+", x21="+x21+", x10="+
    //           x10+", x11="+x11+", v2="+v2;
    //         x3 = x30+(x31-x30)*s;
    //         ibs = is1;
    //       }
    //       interpCount++;
    //     } else {
    //       int is = getSampleIndex(nr,ibs);
    //       x1 = _hx1[is]; x3 = _hx3[is];
    //       ongridCount++;
    //     }
    //     // assert !Float.isNaN(x1);
    //     growX1[i2] = x1;
    //     growX3[i2] = x3;
    //   }
    //   growList.add(new float[][] {growX1,growX3});
    // }
    // print("interpolated="+interpCount+", ongridCount="+ongridCount);

    // // Grid all the columns of horizon values (s3 coordinates).
    // double d3 = _s3.getDelta();
    // double f3 = _s3.getFirst();
    // double v3 = f3;
    // int ngrow = growList.size();
    // float[] hcolX1 = new float[ngrow]; // horizon column x1 values
    // float[] hcolX3 = new float[ngrow]; // horizon column x3 values
    // for (int i2=0; i2<n2; i2++) {
    //   // fill hcol values
    //   for (int irow=0; irow<ngrow; irow++) {
    //     float[][] grow = growList.get(irow);
    //     hcolX1[irow] = grow[0][i2];
    //     hcolX3[irow] = grow[1][i2];
    //   }
    //   int ibs = 0; // binary search index
    //   for (int i3=0; i3<n3; i3++,v3=f3+i3*d3) {
    //     ibs = binarySearch(hcolX1,(float)v3,ibs);
    //     float x1;
    //     if (ibs<0) { // linearly interpolate
    //       int is0 = -ibs-2;
    //       int is1 =  is0+1;
    //       if (is0<0) {
    //         x1 = hcolX1[0];
    //         ibs = 0;
    //       } else if (is1>=ngrow) {
    //         x1 = hcolX1[ngrow-1];
    //         ibs = ngrow-1;
    //       } else {
    //         float x10 = hcolX1[is0];
    //         float x11 = hcolX1[is1];
    //         float x30 = hcolX3[is0];
    //         float x31 = hcolX3[is1];
    //         x1 = x10+(x11-x10)*(((float)v3-x30)/(x31-x30));
    //         ibs = is1;
    //       }
    //     } else {
    //       x1 = hcolX1[ibs];
    //     }
    //     _gh[i3][i2] = x1;
    //   }
    // }
  }


  // public float[] getI2(int i2) {

  // }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  private final Sampling _s2, _s3;
  private final float[][] _gh; // gridded horizon [n3][n2]
  private final int _hns; // number of horizon samples
  private final float[] _hx1, _hx2, _hx3; // horizon samples
  private int _is;

  // Get a row of horizon values
  private float[] getRow() {
    if (_is>=_hns) return null;
    float x3p = _hx3[_is];
    float x2p = _hx2[_is++];
    FloatList x2List = new FloatList();
    x2List.add(x2p);
    print("isFirst: x2v="+x2p+", x3v="+x3p);
    for (; _is<_hns; _is++) {
      float x2v = _hx2[_is];
      float x3v = _hx3[_is];
      if (x2v<x2p) break;
      print("is="+_is+", x2v="+x2v+", x3v="+x3v);
      x2List.add(x2v);
      x2p = x2v;
    }
    return x2List.trim();
  }

  private int getSampleIndex(int nr, int ir) {
    int is = _is-nr+ir;
    return is>=_hns ? _hns-1 : is;
  }

  private void testNaNs(float[] f) {
    int cNaN = 0;
    for (int i=0; i<f.length; i++) {
      if (Float.isNaN(f[i]))
        cNaN++;
    }
    print("Number of NaN values: "+cNaN);
  }

  private static void print(String s) {
    System.out.println(s);
  }

}