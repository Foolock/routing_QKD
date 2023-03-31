#include "grid.hpp"

// test driver
int main() {
  
  // run DIAG-4-4-4 setting
  std::vector<std::vector<int>> T{{5, 3},{3, 5}};

  int N = 10;

  double P = 0.75;

  Grid grid(T, N, P);
  

  return 0;
}

