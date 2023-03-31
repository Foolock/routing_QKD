#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <cmath>

/** 
 * @brief: calculate distance for each node to Alice, Bob, and multiple TN and assign
 *         it to the D vector of each node.
 *          
 *         for each node, the format of their distance to A, B, and multiple TNs
 *         is {Da, Db, Dt1, Dt2, Dt3, ...}
 *  
 */
void calDistanceToABTs(
  int grid_size,
  std::vector<int> A_index,
  std::vector<int> B_index,
  std::vector<std::vector<int>> T_indices,
  std::vector<std::vector<Node>>& temp_node_grid
) {

  // i -> row of grid
  // j -> column of grid
  for(int i=0; i<grid_size; i++) {
    for(int j=0; j<grid_size; j++) {

      // get distance: manhattan distance, and add it to node.D
      temp_node_grid[i][j].D.push_back(std::abs(A_index[0] - i) + std::abs(A_index[1] - j));
      temp_node_grid[i][j].D.push_back(std::abs(B_index[0] - i) + std::abs(B_index[1] - j));
     
      for(int t=0; t<T_indices.size(); t++) {
        temp_node_grid[i][j].D.push_back(std::abs(T_indices[t][0] - i) + std::abs(T_indices[t][1] - j));
      }
    }
  }
}



#endif
