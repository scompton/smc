package gbc;

import static java.lang.Math.*;

/**
 * Coordinate transforms for GBC data.
 * Coordinates in many data files are map coordinates (xe,yn,ze), where
 * xe is easting, yn is northing measured in km, and ze is time in ms.
 */
public class Coordinates {

  // Eastings and northings - (xe,yn) coordinates - of seismic rectangle.
  // Coordinates are measured in km.
  static final double XLL_KM=734.83856, YLL_KM=174.99298;//(i2,i3)=(0,0)
  static final double XLR_KM=731.3053 , YLR_KM=178.52492;//(i2,i3)=(n2-1,0)
  static final double XUR_KM=734.7187 , YUR_KM=181.93951;//(i2,i3)=(n2-1,n3-1)
  static final double XUL_KM=738.2519 , YUL_KM=178.4076 ;//(i2,i3)=(0,n3-1)

  // Rotation angle of normalized seismic rectangle. This is the angle between
  // the normalized seismic x2 (inline) axis and the map xe (easting) axis,
  // measured CCW in the map (xe,yn) coordinate system.
  // static final double PHID = atan2(YLR_KM-YLL_KM,XLR_KM-XLL_KM);
  static final double PHID = atan2(YLR_KM-YLL_KM,XLL_KM-XLR_KM);
  static final double COSD = cos(PHID);
  static final double SIND = sin(PHID);

   /**
   * Map coordinates xe (easting), yn (northing) in km.
   */
  public static class Map {
    public double xe,yn;
    public Map(double xe, double yn) {
      this.xe = xe;
      this.yn = yn;
    }
    public Map(Nrm nrm) {
      double x2 = nrm.x2;
      double x3 = nrm.x3;
      xe = COSD*x2-SIND*x3;
      yn = SIND*x2+COSD*x3;
      xe += XLL_KM;
      yn += YLL_KM;
    }
  }

  /**
   * Normalized coordinates x2 (inline) and x3 (crossline) in km.
   */
  public static class Nrm {
    public double x2,x3;
    public Nrm(double x2, double x3) {
      this.x2 = x2;
      this.x3 = x3;
    }
      public Nrm(Map map) {
      double xe = map.xe;
      double yn = map.yn;
      xe -= XLL_KM;
      yn -= YLL_KM;
      x2 = -COSD*xe+SIND*yn;
      x3 = SIND*xe+COSD*yn;
    }
  }

}