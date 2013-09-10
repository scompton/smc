package warp;

/**
 * An interface used to monitor and report progress for warping. 
 */
public interface WarperWorkTracker {

  /**
   * Set the total number of work units that will be reported.
   * @param totalWorkUnits
   */
  public void setTotalWorkUnits(int totalWorkUnits);

  /**
   * Mark this many {@code units} as worked. This is not a cumulative total,
   * just an incremental number of work units that has been completed.
   * @param units
   */
  public void worked(int units);

  /**
   * Returns whether or not this task is in a canceled state.
   * @return {@code true} if the task is canceled, {@code false} otherwise.
   */
  public boolean isCanceled();

}
