#include <vector>

struct Node {

//  // Node Index
//  std::vector<int> node_index = {0,0};

  // node role: 0 = router, 1 = Alice, 2 = Bob, 3 = TN 
  // At default, a node is a router
  int role = 0;

  // Distance to Alice, Bob and TN
  int Da = 0;
  int Db = 0;
  int Dt = 0;

  // 4 qubit memory
  Node* qubut_up = nullptr;
  Node* qubut_down = nullptr;
  Node* qubut_left = nullptr;
  Node* qubut_right = nullptr;
};
