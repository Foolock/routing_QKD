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
 
  std::cout << "Initialization:\n"; 
  grid.display();

  grid.stage1();
  std::cout << "after stage1: \n";
  grid.display();

  grid.stage2_global();
  grid.stage2_local_IA();

    
  // get 100 simples: compare the performance of Global and Local routing
  int num_sample = 10;
  for(int i=0; i<num_sample; i++) {
    grid.stage1();
    grid.stage2_local_IA();
    grid.stage2_global();
    grid.reset(); 
  }
  std::cout << "after another " << num_sample + 1 << " samples\n";

  int num_key = grid.getMaxFlow(grid.SS_global);
 
  std::cout << "key amount of global routing = " << num_key << "\n";

  grid.displayNetworkGraph();

  grid.getPriorityEdge();
  return 0;
}


