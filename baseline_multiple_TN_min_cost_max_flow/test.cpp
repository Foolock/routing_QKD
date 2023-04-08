#include "grid.hpp"
#include <iomanip>
// test driver
int main() {
  
  // run DIAG-4-4-4 setting
  std::vector<std::vector<int>> T{{5, 3},{4, 4},{3, 5}};

  int N = 9;

  double P = 0.75;

  double B = 0.85;

  double D = 0.02;
  Grid grid(T, N, P, B, D);

  int num_sample = 2000;
  int num_sample_per_period = 500;

  grid.display();

  grid.stage1();
  std::cout << "after stage 1: \n";
  grid.display();
//
//  std::cout << "total number of rounds: " << num_sample << "\n";
//  std::cout << "number of rounds in one period: " << num_sample_per_period << "\n";
// 
//  for(int i=0; i<num_sample; i++) {
//    if(i > 0 && i % num_sample_per_period == 0) {
//      // for a certain period, get prioritized edges
//      grid.getMaxFlow(grid.SS_global);
//      grid.getPriorityEdge();
//    }
//    std::cout << "dynamic round: " << i << "\n";
//    grid.stage1();
//    grid.stage2_global_dynamic();
//    grid.reset(); 
//  }
//
//  Grid grid1(T, N, P, B, D); 
//
//  for(int i=0; i<num_sample; i++) {
//    std::cout << "static round: " << i << "\n";
//    grid1.stage1();
//    grid1.stage2_global_static();
//    grid1.reset();
//  }
//
//  int num_key_dynamic = grid.getMaxFlow(grid.SS_global);
//  int num_key_static = grid1.getMaxFlow(grid1.SS_global);
//
//  std::cout << "key amount of dynamic local routing = " << num_key_dynamic << "\n";
//  std::cout << "key amount of static local routing = " << num_key_static << "\n";


  return 0;
}



