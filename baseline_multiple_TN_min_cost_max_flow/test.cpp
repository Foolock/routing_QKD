#include "grid.hpp"
#include <iomanip>
#include <chrono>
// test driver
int main() {
  
  // run DIAG-4-4-4 setting
  std::vector<std::vector<int>> T{{6, 2},{3, 4}};

  int N = 9;

  double P = 0.75;

  double B = 0.85;

  double D = 0.02;
  Grid grid(T, N, P, B, D);

  int num_sample = 10;
  int num_sample_per_period = 1000;

  grid.display();

  Grid grid1(T, N, P, B, D); 

  auto beg = std::chrono::steady_clock::now();
  for(int i=0; i<num_sample; i++) {
    if(i > 0 && i % num_sample_per_period == 0) {
      // for a certain period, get prioritized edges
      grid1.getMaxFlow(grid.SS_global);
      grid1.getPriorityEdge();
    }
    std::cout << "dynamic round: " << i << "\n";
    
    grid1.stage1();
    grid1.stage2_global_dynamic();
    grid1.reset();
  }
  auto end = std::chrono::steady_clock::now();

  int num_key_dynamic = grid1.getMaxFlow(grid1.SS_global);

  auto beg2 = std::chrono::steady_clock::now();
  for(int i=0; i<num_sample; i++) {
    std::cout << "mcmf round: " << i << "\n";
    grid.stage1();
    grid.stage2_min_cost_max_flow();
    grid.reset(); 
  }
  auto end2 = std::chrono::steady_clock::now();

  int num_key_mcmf = grid.getMaxFlow(grid.SS_MCMF);
  
  std::cout << "key amount of mcmf routing = " << num_key_mcmf << "\n";
 
  std::cout << "key amount of dynamic global routing = " << num_key_dynamic << "\n";

  std::cout << "time of " << num_sample << " rounds of dynamic global routing is: " 
            << std::chrono::duration_cast<std::chrono::microseconds>(beg-end).count()
            << " ms\n";
 
  std::cout << "time of " << num_sample << " rounds mcmf routing is: " 
            << std::chrono::duration_cast<std::chrono::microseconds>(beg2-end2).count()
            << " ms\n";
  
 
  return 0;
}



