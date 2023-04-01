#include "grid.hpp"
// test driver
int main() {
  
  // run DIAG-4-4-4 setting
  std::vector<std::vector<int>> T{{5, 3},{3, 5}};

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
  int num_sample = 100;
  for(int i=0; i<num_sample; i++) {
    grid.stage1();
    grid.stage2_local_IA();
    grid.stage2_global();
    grid.reset(); 
  }

  std::cout << "after another " << num_sample << " samples\n";

  std::cout << "Global routing: \n";
  std::cout << "the num of paths found between A and T1: ";
  std::cout << grid.SS_global[0][2].size() << "\n";

  std::cout << "the num of paths found between A and T2: ";
  std::cout << grid.SS_global[0][3].size() << "\n";

  std::cout << "the num of paths found between A and B: ";
  std::cout << grid.SS_global[0][1].size() << "\n";

  std::cout << "the num of paths found between T1 and T2: ";
  std::cout << grid.SS_global[2][3].size() << "\n";

  std::cout << "the num of paths found between T1 and B: ";
  std::cout << grid.SS_global[2][1].size() << "\n";

  std::cout << "the num of paths found between T2 and B: ";
  std::cout << grid.SS_global[3][1].size() << "\n";
  
  std::cout << "Local routing: \n";
  std::cout << "the num of paths found between A and T1: ";
  std::cout << grid.SS_local[0][2].size() << "\n";

  std::cout << "the num of paths found between A and T2: ";
  std::cout << grid.SS_local[0][3].size() << "\n";

  std::cout << "the num of paths found between A and B: ";
  std::cout << grid.SS_local[0][1].size() << "\n";

  std::cout << "the num of paths found between T1 and T2: ";
  std::cout << grid.SS_local[2][3].size() << "\n";

  std::cout << "the num of paths found between T1 and B: ";
  std::cout << grid.SS_local[2][1].size() << "\n";

  std::cout << "the num of paths found between T2 and B: ";
  std::cout << grid.SS_local[3][1].size() << "\n";
  

  
  std::cout << "key amount of global routing = " << grid.getMaxFlow(grid.SS_global) << "\n";


  return 0;
}


