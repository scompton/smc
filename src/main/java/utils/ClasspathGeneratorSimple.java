package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class ClasspathGeneratorSimple {
  
  public ClasspathGeneratorSimple(String projectName) {
    File project = new File(HOME,projectName);
    if (project.exists() && project.isDirectory()) {
      addEntries(project);
    }
  }
  
  public static void main(String[] args) {
    if (args.length!=1) {
      print("Usage: java ClasspathGenerator packageName");
      System.exit(0);
    }
    String packageName = args[0];
    new ClasspathGeneratorSimple(packageName);
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // Private

  private static final String HOME = "/Users/scompton/svn/trunk2/";
  private static final String DEPLOY_ENTRY = "<classpathentry kind=\"con\" "+
      "path=\"org.eclipse.jdt.USER_LIBRARY/deploy\"/>";
  
  private List<String> addEntries(File project) {
    File cp = new File(project,".classpath");
    try {
      StringBuffer sb = new StringBuffer();
      sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      sb.append("<classpath>\n");
      sb.append("        <classpathentry kind=\"src\" path=\"src\"/>\n");
      sb.append("        <classpathentry kind=\"con\" path=\"org.eclipse.jdt.launching.JRE_CONTAINER\"/>\n");
      sb.append("        <classpathentry kind=\"con\" path=\"org.eclipse.pde.core.requiredPlugins\"/>\n");
      sb.append("        <classpathentry kind=\"output\" path=\"build/class\"/>\n");
      sb.append("        "+DEPLOY_ENTRY+"\n");  
      sb.append("</classpath>");
      BufferedWriter bw = new BufferedWriter(new FileWriter(cp));
      bw.write(sb.toString());
      bw.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return null;
  }
  
  private static void print(String s) {
    System.out.println(s);
  }
}
