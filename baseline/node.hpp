#include <vector>

struct Qubit {

  // if a qubit is connected with an inter link, available = true
  bool available = false;

  // the index of qubit this qubit is linked to
  int to = -1;

};


struct Node {

  // node role: 0 = router, 1 = Alice, 2 = Bob, 3 = TN 
  // At default, a node is a router
  int role = 0;

  // Distance to Alice, Bob and TN
  int Da = 0;
  int Db = 0;
  int Dt = 0;

  // 4 qubit memory
  // 0 = uppper, 1 = left, 2 = right, 3 = bottom
  std::vector<Qubit> qubits = std::vector<Qubit>(4);

  // a vector of 4 direction
  // when direction[i] = true
  // it is indicatiing that direction i is connect with inter links
  // direction: 0 = uppper, 1 = left, 2 = right, 3 = bottom
  std::vector<bool> direction{std::vector<bool>(4, false)};

};

