#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "node.hpp"
#include "edge.hpp"
#include <cmath>
#include <queue>
#include <algorithm>

/**
 *
 *
 *
 *    ////////////////function declaration/////////////////////////
 *
 *
 *
 *
 *
 */

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
);

/**
 * @brief: helper: check if Da, Db, Dt is calculated correctly 
 */
std::vector<int> getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid);

/**
 * @brief: helper: transfer a node index from a integer form to a (x, y) form 
 *
 * input:
 *  node's index (int)
 * return(std::vector<int>):
 *  node's index : (x, y)
 */
std::vector<int> int2coordinate(int node_index, int grid_size);

/**
 * @brief(need more case test): helper: bfs, traverse the current node_grid_per_round. Find available path between s, t
 * https://www.geeksforgeeks.org/print-paths-given-source-destination-using-bfs/
 *
 * input:
 *  s, t: index(integer) of source and sink node. Notice: not the role of node
 * 
 */
std::vector<std::vector<int>> bfs(
    int s, int t,
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size
    );

// @brief: helper: check if a node is in the path vector
bool isInPath(int x, std::vector<int> path);


/**
 *
 *
 *
 *    ////////////////function definition/////////////////////////
 *
 *
 *
 *
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

std::vector<int> getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid) {
  std::vector<int> result; // result = {Da, Db, Dt1, Dt2, ...}
  result = temp_node_grid[x][y].D;
  return result;
}

std::vector<int> int2coordinate(int node_index, int grid_size) {
  std::vector<int> result(2); // result = {row, col}
  result[0] = node_index / grid_size;
  result[1] = node_index % grid_size;
  return result;
}

std::vector<std::vector<int>> bfs(
    int s, int t, 
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size
    ) {
  
  // operator for std::min_element: comparison function object  
  auto min_size = []( const auto &v1, const auto &v2 )
  {
      return std::size( v1 ) < std::size( v2 );
  };

  // a result vector to store the path we found
  // result[i] = a path found
  std::vector<std::vector<int>> result;

  // create a queue for bfs
  // which stores the paths
  std::queue<std::vector<int>> q;

  // path vector to store the current path
  std::vector<int> path;
  path.push_back(s);
  q.push(path);
  while(!q.empty()) {
    path = q.front();
    q.pop();

    auto it = std::min_element(std::begin(result), std::end(result), min_size);
    if(result.size()) { // "it" is nullptr at the beginning cuz nothing in result
      if(path.size() > it->size()) { // if the path is already longer than the result we have
                                     // no need to continue the loop
        continue;
      }
    } 

    int last = path[path.size() - 1];
    
    // if last vertex is the desired destination 
    // then store this path to our result path vector 
    if(last == t) {
      result.push_back(path);
    }

    // traverse to all the nodes connected to current node
    // if adjacent node is not visited, 
    // create a new path by copying the current path
    // and add this adjacent node into the new path
    // and push the new path to queue
    for(int i=0; i<edges_per_round[last].size(); i++) { // edges_per_round[last] = 4 cuz 4 neighbor(but some entries may be -1)
      if(edges_per_round[last][i].to != -1) {
        std::vector<int> coordinate = int2coordinate(edges_per_round[last][i].to, grid_size); 
        if(node_grid_per_round[coordinate[0]][coordinate[1]].role == 0 || edges_per_round[last][i].to == t) { 
          // if this index is some other
          // TN or user node, skip
          if(!isInPath(edges_per_round[last][i].to, path)) {
            std::vector<int> newpath = path;
            newpath.push_back(edges_per_round[last][i].to);
            q.push(newpath);
          }
        }
      }
    }
  }
  
  return result; 

}

bool isInPath(int x, std::vector<int> path)
{
  // x may be -1 for a non-existing adjacent node
  if(x != -1) {
    for(int i=0; i<path.size(); i++) {
      if(path[i] == x) {
        return true;
      }
    }
  }
  return false;
}


#endif
