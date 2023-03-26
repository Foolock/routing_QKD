#include "grid.hpp"


// test driver
int main() {
  
  std::vector<int> T = {2, 2};
  int N = 5;
  double P = 0.75;

  Grid grid(T, N, P);
  
  std::cout << "before stage 1: \n";
  grid.display();
  
  grid.stage1();
  std::cout << "after stage 1: \n";
  grid.display();

  grid.stage2_global();

  return 0;
}


