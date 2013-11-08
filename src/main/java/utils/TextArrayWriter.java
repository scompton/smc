package utils;

import java.io.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class TextArrayWriter {

  public TextArrayWriter(float[] f, String fileName) {
    try {
      BufferedWriter bw = new BufferedWriter(new FileWriter(fileName)); 
      int n = f.length;
      for (int i=0; i<n; i++) {
        String s = i!=n-1 ? String.valueOf(f[i])+"\n" : String.valueOf(f[i]);
        bw.write(s,0,s.length());
      }
      bw.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

}
