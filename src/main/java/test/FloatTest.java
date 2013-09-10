package test;

public class FloatTest {

  public static void main(String[] args) {
    float a = 1E29f;
//    float b = 1E31f;
    float b = Float.MAX_VALUE;
//    float b = 0.1f;
    float c = a+b;
    print("a="+a+", b="+b+", c=a+b="+c);
  }

  private static void print(String s) {
    System.out.println(s);
  }
}
