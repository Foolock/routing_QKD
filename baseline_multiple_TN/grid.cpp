#include "grid.hpp"
#include "utility.hpp"

//
// @brief: (constructor) place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
Grid::Grid(std::vector<std::vector<int>>& TN_locations, int N, double P): 
  grid_size(N), 
  A_index{grid_size-1, 0},
  B_index{0, grid_size-1}, 
  T_indices(TN_locations), 
  P(P)
{
  // create a temporary object to store node_grid_per_round
  std::vector<std::vector<Node>> temp_node_grid(grid_size, std::vector<Node>(grid_size));   

  // place A, B in graph
  temp_node_grid[A_index[0]][A_index[1]].role = 1;
  temp_node_grid[B_index[0]][B_index[1]].role = 2;

  // place multiple Trusted Nodes in graph
  for(int i=0; i<T_indices.size(); i++) {
    temp_node_grid[T_indices[i][0]][T_indices[i][1]].role = 3 + i; // the role of TN ranges from 3 to ... 
  }

  // initialze distance to Alice, Bob and multiple TNs for each node
  calDistanceToABTs(grid_size, A_index, B_index, T_indices, temp_node_grid);


}
