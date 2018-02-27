#ifndef _SgeGridMonitor_H_
#define _SgeGridMonitor_H_

// Module Includes

// External Module Includes

// System Includes
#include <map>
#include <string>

namespace specnets
{
  class SgeGridMonitor
  {
  public:
    enum JobStatus
    {
      JS_UNKNOWN = 0,
      JS_DELETION = 1,
      JS_ERROR = 2,
      JS_HOLD = 3,
      JS_RUNNING = 4,
      JS_RESTARTED = 5,
      JS_SUSPENDED = 6,
      JS_WAITING = 7,
      JS_DONE = 8
    };

    SgeGridMonitor(void);
    ~SgeGridMonitor(void);

    JobStatus getJobStatus(const std::string & jobId);
    JobStatus getJobStatus(int jobId);

    bool refreshInfo(const std::string & sgePath);

    std::string submitJob(const std::string & sgePath,
                          const std::string & jobFilename,
                          const std::string & params,
                          const std::string & gridType);

  private:
    std::map<std::string, JobStatus> jobs;
  };

} //namespace specnets

#endif // _GridUtils_H_

