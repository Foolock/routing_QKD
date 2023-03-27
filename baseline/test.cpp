#include "grid.hpp"


// test driver
int main() {
  
  std::vector<int> T = {2, 2};
  int N = 5;
  double P = 0.75;

  Grid grid(T, N, P);
  
  std::cout << "before stage 1: \n";
  grid.display();
//  
//  grid.stage1();
//  std::cout << "after stage 1: \n";
//  grid.display();
//
//  grid.stage2_global();
//
//  std::cout << "check SS : \n";
//  std::cout << "SSab: ";
//  for(auto& p : grid.SSab) {
//    std::cout << p << " ";
//  }
//  std::cout << "\n";

  // just for testing, remove the right edge of the node(2,2) 
  // and see if find2qubits_IA works 
//  grid.breakEdge(2, 2, 2);
  std::vector<int> available_q{0,1,3};
  grid.find2qubits_IA(2, 2, available_q);

  return 0;
}


