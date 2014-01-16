package utils;

import java.awt.*;
import javax.swing.*;

public class TestPanel {

  public static void main(String[] args) {
    JTextArea text = new JTextArea("test");
    JPanel panel = new JPanel();
    //panel.setBackground(Color.BLACK);
    panel.add(text);
    
    JFrame frame = new JFrame();
    frame.add(panel);
    frame.setBackground(Color.BLACK);
    frame.setSize(800,800);
    frame.setVisible(true);
  }

}
