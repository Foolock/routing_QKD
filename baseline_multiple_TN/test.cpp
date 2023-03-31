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

  std::cout << "the paths found between A and T1: ";
  for(auto p : grid.SS[0][2]) {
    std::cout << p << " ";
  }
  std::cout << "\n";

  std::cout << "the paths found between A and T2: ";
  for(auto p : grid.SS[0][3]) {
    std::cout << p << " ";
  }
  std::cout << "\n";

  std::cout << "the paths found between A and B: ";
  for(auto p : grid.SS[0][1]) {
    std::cout << p << " ";
  }
  std::cout << "\n";

  std::cout << "the paths found between T1 and T2: ";
  for(auto p : grid.SS[2][3]) {
    std::cout << p << " ";
  }
  std::cout << "\n";


  std::cout << "the paths found between T1 and B: ";
  for(auto p : grid.SS[2][1]) {
    std::cout << p << " ";
  }
  std::cout << "\n";


  std::cout << "the paths found between T2 and B: ";
  for(auto p : grid.SS[3][1]) {
    std::cout << p << " ";
  }
  std::cout << "\n";


  return 0;
}


