// author : Jane Tournois

#include <ctime>
#include <vector>
#include <fstream>
#include <iostream>

class Meshing_timers 
{
public:
  enum Chrono { REFINEMENT = 0,
                AABB_TREE_BUILD,
                KD_TREE_BUILD,
                FILL_QUEUE };

private:
  static const unsigned int NB_TIMERS = 4;

private:
  std::vector<time_t> timers;
  std::vector<bool> started;
  std::vector<time_t> buffers;
  
public:
  Meshing_timers() :
      timers(NB_TIMERS),
      started(NB_TIMERS),
      buffers(NB_TIMERS)
  {
    for(unsigned int i = 0; i < NB_TIMERS; i++)
      started[i] = false;
  };

public: // these operations should always take place in the following order
  // (pause and resume can be skipped or used several times but only together)
  
  void start(const Chrono& c)
  {
    timers[c] = 0;
    started[c] = true;
    buffers[c] = time(NULL);
  }
  
  void pause(const Chrono& c)
  {
    if(started[c])
      timers[c] += (time(NULL) - buffers[c]);
    else
      std::cout << "warning : timer paused before started.\n";
  }
  
  void resume(const Chrono& c)
  {
    if(started[c])
      buffers[c] = time(NULL);
    else
      start(c);
  }

  void end(const Chrono& c)
  {
    if(started[c])
      pause(c); 
    started[c] = false;
  }

  void end_all()
  {
    for(unsigned int i = 0; i < NB_TIMERS; i++)
      end(Chrono(i));
  }

public:
  time_t get(const Chrono& c)
  {
    return timers[c];
  }

  void report_timers()
  {
    std::ofstream fs("times.txt");
    fs << "Refinement :\t"      << timers[REFINEMENT] << "\n";
    fs << "AABB_tree build :\t" << timers[AABB_TREE_BUILD] << "\n";
    fs << "Kd_tree build :\t"   << timers[KD_TREE_BUILD] << "\n";
    fs << "Fill queue :\t"      << timers[FILL_QUEUE] << "\n";
    fs.close();
  }
};