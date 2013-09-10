package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ClasspathGenerator {
  
  public ClasspathGenerator(String projectName) {
    File project = new File(HOME,projectName);
    if (project.exists() && project.isDirectory()) {
      File manifest = new File(project.getAbsolutePath()+MANIFEST);
      if (manifest.exists()) {
        StringBuffer sb = new StringBuffer();
        try {
          BufferedReader br = new BufferedReader(new FileReader(manifest));
          String line;
          boolean collect = false;
          
          while ((line = br.readLine())!=null) {
            if (line.startsWith(START)) {
              collect = true;
            }
            
            if (line.startsWith(FINISH)) {
              collect = false;
            }
            
            if (collect)
              sb.append(line);
          }
          
          br.close();
        } catch (FileNotFoundException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        
//        List<String> names = getExistingNames(project);
        String collected = sb.toString();
        if (collected.isEmpty()) return;
        String bundle = collected.split(":")[1].trim().replaceAll("\\s","");
        boolean eclipse = false;
        if (bundle.contains(ECLIPSE))
          eclipse = true;
        String[] packages = bundle.split(",");
        List<String> cpEntries = new ArrayList<String>();
        if (eclipse) 
          cpEntries.add(ECLIPSE_ENTRY);
        cpEntries.add(JUNIT_ENTRY);
        for (String p : packages) {
          if (p.contains(X4M)) {
            p = p.split(X4M)[1];
            String[] entries = getClasspathEntry(p);
            if (entries==null) continue;
            for (String entry : entries) {
              if (entry==null) continue;
              cpEntries.add(entry);
            }  
          }
        }
        addEntries(project,cpEntries);
      }
    }
    
  }
  
  public static void main(String[] args) {
    if (args.length!=1) {
      print("Usage: java ClasspathGenerator packageName");
      System.exit(0);
    }
    String packageName = args[0];
    new ClasspathGenerator(packageName);
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // Private

  private static final String HOME = "/Users/scompton/svn/trunk2/";
  private static final String MANIFEST = "/META-INF/MANIFEST.MF";
  private static final String START = "Require-Bundle:";
  private static final String FINISH = "Export-Package:";
  private static final String ECLIPSE = "org.eclipse";
  private static final String X4M = "com.x4m.";
  private static final String SRC = "src";
  private static final String JAR = "jar";
  private static final String ECLIPSE_ENTRY = "<classpathentry kind=\"con\" "+
              "path=\"org.eclipse.jdt.USER_LIBRARY/Eclipse_Plugins_3.7.2\"/>";
  private static final String JUNIT_ENTRY = "<classpathentry kind=\"lib\" "+
              "path=\"../../tools/ant/apache-ant-1.7.0/lib/junit-4.4.jar\"/>";
  
  private String[] getClasspathEntry(String x4mName) {
    File x4mFile = new File(HOME,x4mName);
    String[] entries = null;
    if (x4mFile.exists() || x4mFile.isDirectory()) {
      
      boolean hasSrc = false;
      boolean hasJar = false;
      String[] fileList = x4mFile.list();
      for (String dir : fileList) {
        if (dir.contentEquals(SRC) && !x4mName.contains("jves")) {
          hasSrc = true;
        }
        if (dir.contentEquals(JAR)) {
          hasJar = true;
        }
      }
      
      
      if (hasSrc) {
        entries = new String[1];
        entries[0] = "<classpathentry combineaccessrules=\"false\" kind=\"src\" "+ 
                     "path=\"/"+x4mName+"\"/>";
      } else if (hasJar) {
        File[] jarFiles = new File(x4mFile,"jar").listFiles();
        entries = new String[jarFiles.length];
        for (int i=0; i<jarFiles.length; i++) {
          String name = jarFiles[i].getName();
          if (name.contains(".svn") || name.contains("idh.jar"))
            continue;
          entries[i] = "<classpathentry kind=\"lib\" "+ 
                       "path=\"../"+x4mName+"/jar/"+name+"\"/>";
        }
      } else if (x4mName.contains("jep")) {
        entries = new String[1];
        entries[0] = "<classpathentry kind=\"lib\" "+ 
                     "path=\"../"+x4mName+"/jep-java-3.4.jar\"/>";
      }
    }
    return entries;
  }
  
  private List<String> getExistingNames(File project) {
    File cp = new File(project,".classpath");
    try {
      BufferedReader br = new BufferedReader(new FileReader(cp));
      List<String> names = new ArrayList<String>();
      String line;
      while ((line = br.readLine())!=null) {
        if (line.contains("path=")) {
          String path = line.split("path=")[1];
          path = path.substring(path.indexOf("\"")+1,path.lastIndexOf("\""));
          String[] parts = path.split(File.separator);
          String name = parts[parts.length-1];
          names.add(name);
        }
      }
      br.close();
      return names;
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return null;
  }
  
  private List<String> addEntries(File project, List<String> entries) {
    File cp = new File(project,".classpath");
    try {
      StringBuffer sb = new StringBuffer();
      sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      sb.append("<classpath>\n");
      sb.append("        <classpathentry kind=\"src\" path=\"src\"/>\n");
      sb.append("        <classpathentry kind=\"con\" path=\"org.eclipse.jdt.launching.JRE_CONTAINER\"/>\n");
      sb.append("        <classpathentry kind=\"con\" path=\"org.eclipse.pde.core.requiredPlugins\"/>\n");
      sb.append("        <classpathentry kind=\"output\" path=\"build/class\"/>\n");
      for (String entry : entries) {
        sb.append("        "+entry+"\n");  
      }
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
