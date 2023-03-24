#include "grid.hpp"


// test driver
int main() {
  
  std::vector<int> T = {2, 2};
  size_t N = 5;

  Grid grid(T, N);
  grid.display();

  return 0;
}


