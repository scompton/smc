package utils;

import java.io.*;
import java.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

public class TextArrayReader2 {

  public TextArrayReader2(String fileName) {
    try {
      BufferedReader br = new BufferedReader(new FileReader(fileName)); 
      FloatList fl = new FloatList();
      String line = null;
      while ((line = br.readLine()) != null) {
        String[] sa = line.split("\\s+");
        int n = sa.length;
        float[] f1 = new float[n];
        for (int i=0; i<n; i++)
          f1[i] = Float.valueOf(sa[i]);
        _fl.add(f1);
      }
      br.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public float[][] getArray() {
    int n2 = _fl.size();
    float[][] f = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      f[i2] = _fl.get(i2);
    return f;
  }

  public static void main(String[] args) {
    TextArrayReader2 fr = new TextArrayReader2(args[0]);
    float[][] f = fr.getArray();
    dump(f);
  }

  private List<float[]> _fl = new ArrayList<>();
}
