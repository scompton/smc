package warp;

public class FractionalShifts {
  
  public static void main(String[] args) {
    if (args.length==1)
      test1(Float.parseFloat(args[0]));
    if (args.length==2)
      test2(Float.parseFloat(args[0]),Integer.parseInt(args[1]));
    if (args.length==3)
      test3(Float.parseFloat(args[0]),Float.parseFloat(args[1]),
          Integer.parseInt(args[2]));
  }

  public static void test1(float vpvs) {
    int n1 = 1015;
    float scale = 2.0f/(vpvs+1.0f);
    int nl = n1-(int)(n1*scale);
    int n1max = n1-nl+1;
    System.out.println("n1="+n1+", n1max="+n1max+", nl="+nl);
    for (int i1=0; i1<n1max; ++i1) {
      int count = 0;
      for (int il=0, jl=i1; il<nl; ++il,++jl) {
        count++;
//        System.out.println("i1="+i1+", jl="+jl);
        if (jl>=n1)
          throw new RuntimeException("Error! Out of Bounds");
        if (i1==(n1max-1) && il==(nl-1))
          if (jl!=(n1-1))
            throw new RuntimeException("Error! Not at end");
      }
      if (count!=nl)
        throw new RuntimeException("Error! Too few lags. count="+count);
    }
  }
  
  public static void test2(float vpvs, int fr) {
    int n1 = 1001;
    float scale = 2.0f/(vpvs+1.0f);
    int nl = n1-(int)(n1*scale);
//    int nl = 4;
    int n1max = n1-nl+1;
    nl = (nl-1)*fr+1;
    int n1interp = (n1-1)*fr+1;
    System.out.println("n1="+n1+", n1max="+n1max+", n1interp="+n1interp+
        ", nl="+nl);
    for (int i1=0; i1<n1max; ++i1) {
      int count = 0;
      for (int il=0, jl=i1*fr; il<nl; ++il,++jl) {
        count++;
//        System.out.println("i1="+i1+", jl="+jl);
        if (jl>=n1interp)
          throw new RuntimeException("Error! Out of Bounds");
        if (i1==(n1max-1) && il==(nl-1))
          if (jl!=(n1interp-1))
            throw new RuntimeException("Error! Not at end");
      }
      if (count!=nl)
        throw new RuntimeException("Error! Too few lags. count="+count);
    }
  }
  
  public static void test3(float vpvs, float fsp, int fr) {
    int n1 = 100;
    float scale = 2.0f/(vpvs+1.0f);
    int nl1 = n1-(int)(n1*scale);
    int nlS = (int)(nl1*fsp);
//    int nlS = 2;
    int n1max = n1-nl1-nlS+2;
    nlS = (nlS-1)*fr+1;
    int n1interp = (n1-1)*fr+1;
    System.out.println("n1="+n1+", n1max="+n1max+", nl1="+nl1+", nlS="+nlS);
    for (int i1=0; i1<n1max; ++i1) {
      int count1 = 0;
      for (int il=0, jl=i1; il<nl1; ++il,++jl) {
        count1++;
//        System.out.println("i1="+i1+", jl="+jl);
        if (jl>=n1)
          throw new RuntimeException("Error! Out of Bounds");
        int countS = 0;
        for (int ilS=0,jS=jl*fr; ilS<nlS; ++ilS,++jS) {
          countS++;
//          System.out.println("i1="+i1+", il="+il+", jS="+jS);
          if (jS>=n1interp)
            throw new RuntimeException("Error! Out of Bounds");
          if (i1==(n1max-1) && il==(nl1-1) && ilS==(nlS-1))
            if (jS!=(n1interp-1))
              throw new RuntimeException("Error! Not at end. value="+jS);
        }
        if (countS!=nlS)
          throw new RuntimeException("Error! Too few lags. count="+countS);
      }
      if (count1!=nl1)
        throw new RuntimeException("Error! Too few lags. count="+count1);
    }
  }
}
