package warp;

/**
 * Reports percentage progress for dynamic warping to the command line.
 * @author Stefan Compton
 */
public class WarperProgress implements WarperWorkTracker {

  @Override
  public void setTotalWorkUnits(int totalWorkUnits) {
    _totalWorkUnits = (float)totalWorkUnits;
    _units = 0; // reset completed units

    Thread thread = new Thread(new Runnable() {
      @Override
      public void run() {
        while (_units!=_totalWorkUnits) {
          try {
            updateProgress();
            Thread.sleep(1000);
          } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
          }
        }
        println("\r100% complete");
      }
    });
    thread.start();
  }

  @Override
  public void worked(int units) {
    _units += units; // increment units completed
  }

  @Override
  public boolean isCanceled() {
    return false; // job cancelation not supported.
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private float _totalWorkUnits;
  private int _units; // units completed so far

  private void updateProgress() {
    int p = (int)(100.0f*(_units/_totalWorkUnits));
    print("\r"+p+"% complete");
  }

  private static void print(String s) {
    System.out.print(s);
  }

  private static void println(String s) {
    System.out.println(s);
  }

}
