#include "grid.hpp"

// test driver
int main() {
  
  // run DIAG-4-4-4 setting
  std::vector<std::vector<int>> T{{5, 3},{3, 5}};

  int N = 9;

  double P = 0.75;

  Grid grid(T, N, P);
 
  std::cout << "Initialization:\n"; 
  grid.display();

  grid.stage1();
  std::cout << "after stage1: \n";
  grid.display();

  grid.stage2_global();

  return 0;
}

