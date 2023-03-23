#include <utility>

struct Node {

  // node index
  std::pair<int, int> index;
  
  // node role: 0 = Alice, 1 = Bob, 2 = TN, 3 = router 
  int role;

  // Distance to Alice, Bob and TN
  int Da;
  int Db;
  int Dt;

  // 4 qubit memory
  Node* qubut_up;
  Node* qubut_down;
  Node* qubut_left;
  Node* qubut_right;
};
