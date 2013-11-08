package utils;

import java.io.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class TextArrayReader {

  public TextArrayReader(String fileName) {
    try {
      BufferedReader br = new BufferedReader(new FileReader(fileName)); 
      FloatList fl = new FloatList();
      String line = null;
      while ((line = br.readLine()) != null) {
        fl.add(Float.valueOf(line));
      }
      br.close();
      _f = fl.trim();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public float[] getArray() {
    return _f;
  }

  public static void main(String[] args) {
    TextArrayReader fr = new TextArrayReader(args[0]);
    float[] f = fr.getArray();
    dump(f);
  }

  private float[] _f;
}
