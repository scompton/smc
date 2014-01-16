package fault;

public class Ray {

  public Point p1, p2;

  public Ray(Point p1, Point p2) {
    this.p1 = p1;
    this.p2 = p2;
  }

  public boolean doIntersect(Ray b) {
    return getIntersection(b)!=null;
  }

  // For testing. Ignores the z coordinates and finds the intersection
  // assuming 2D rays.
  public Point getIntersection2(Ray b) {
    float ax21 = this.p2.x - this.p1.x;
    float ay21 = this.p2.y - this.p1.y;
    float bx21 = b.p2.x - b.p1.x;
    float by21 = b.p2.y - b.p1.y;
    //print("ax21="+ax21+",ay21="+ay21);
    //print("bx21="+bx21+",by21="+by21);
    float d = ax21*by21 - ay21*bx21;
    //print("d="+d);
    if (d==0.0f) return null;
    float bax = b.p1.x - this.p1.x;
    float bay = b.p1.y - this.p1.y;
    float n = bax*by21 - bay*bx21;
    float t = n/d;
    print("t="+t);
    float x = this.p1.x + t*(this.p2.x - this.p1.x);
    float y = this.p1.y + t*(this.p2.y - this.p1.y);
    return new Point(x,y,0.0f);
  }

  /**
   * Get the ray that defines the intersection between Ray @code{a} (this 
   * instance) and Ray @code{b}. In 3D two rays will generally not intersect at
   * a point. The shortest line segment that connects @code{a} and @code{b} is 
   * unique and can be considered their intersection in 3D.
   * @return the intersecting Ray.
   */
  public Ray getIntersection(Ray b) {
    float ax21 = this.p2.x - this.p1.x;
    float ay21 = this.p2.y - this.p1.y;
    float az21 = this.p2.z - this.p1.z;
    float bx21 = b.p2.x - b.p1.x;
    float by21 = b.p2.y - b.p1.y;
    float bz21 = b.p2.z - b.p1.z;
    float abx1 = this.p1.x - b.p1.x;
    float aby1 = this.p1.y - b.p1.y;
    float abz1 = this.p1.z - b.p1.z;
    float ab   = ax21*bx21 + ay21*by21 + az21*bz21;
    float a2   = ax21*ax21 + ay21*ay21 + az21*az21;
    float b2   = bx21*bx21 + by21*by21 + bz21*bz21;
    //print("ax21="+ax21+",ay21="+ay21+",az21="+az21);
    //print("bx21="+bx21+",by21="+by21+",bz21="+bz21);
    float ta = 
      ((abx1*bx21 + aby1*by21 + abz1*bz21)*ab -
       (abx1*ax21 + aby1*ay21 + abz1*az21)*b2) / ((a2*b2) - (ab*ab));
    float tb = (abx1*bx21 + aby1*by21 + abz1*bz21 + ta*ab)/b2;
    // TODO if ta or tb is less than 0 or greater than 1, the rays do not
    // intersect. Still, we can find where a and b would intersect if the
    // ray lengths were extended. Do we care where they would intersect, or
    // just whether or not the existing rays intersect?
    print("ta="+ta+", tb="+tb);
    float cx = this.p1.x + ta*(this.p2.x - this.p1.x);
    float cy = this.p1.y + ta*(this.p2.y - this.p1.y);
    float cz = this.p1.z + ta*(this.p2.z - this.p1.z);
    float dx =    b.p1.x + tb*(   b.p2.x -    b.p1.x);
    float dy =    b.p1.y + tb*(   b.p2.y -    b.p1.y);
    float dz =    b.p1.z + tb*(   b.p2.z -    b.p1.z);
    return new Ray(new Point(cx,cy,cz),new Point(dx,dy,dz));
  }

  /**
   * Get the squared distance between points that define this ray.
   * @return squared distance between points that define this ray.
   */
  public float getDistanceSqrd() {
    float dx = this.p2.x - this.p1.x;
    float dy = this.p2.y - this.p1.y;
    float dz = this.p2.z - this.p1.z;
    return dx*dx + dy*dy + dz*dz;
  }

  @Override
  public String toString() {
    return "p1="+this.p1+", p2="+this.p2;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  ////////////////////////////////////////////////////////////////////////////

  private static void print(String s) {
    System.out.println(s);
  }

}
