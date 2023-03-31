#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "node.hpp"
#include "edge.hpp"
#include <cmath>
#include <queue>
#include <algorithm>
#include <limits.h>

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
 * @brief: helper: calculate how many qubits are available for each node
 *
 * input:
 *  x, y: index(coordinate) of a node in the node grid
 *
 * return(int):
 *  result: the number of qubits that are available for a node 
 *
 */
int numTakenQubits(int x, int y, std::vector<std::vector<Node>> node_grid_per_round);

/**
 * @brief: IA algorithm
 *
 * input:
 *  curr_r, curr_c: index of current node
 *  available_q: a vector storing the index of available qubits in qubits of each node
 *  this index should be the same as the index of direction (!= -1) in directions of each node
 *  so we can find neighbor node index through edges_per_round[curr][direction].to 
 *  with this available_q 
 *
 * return:
 *  2 qubit index indicating which 2 qubits to connect as intra link
 *
 */
std::vector<int> find2qubits_IA(int curr_r, int curr_c, std::vector<int> available_q, 
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size);

/**
 * @brief: helper: dfs, recurr from one sink(set as Alice or TN) until it meets a sink as another user node(Bob) or TN 
 *
 * input:
 *  s: index(integer) of source node
 *  target: index(integer) of the sink node
 *  curr_q: index([0,3]) of current qubit
 *
 * return(void):
 *  push the index of node to a path starting from s to TN or Bob or end node
 */
void dfs(int s, int target, int curr_q, std::vector<int>& path,
      std::vector<std::vector<Node>> node_grid_per_round,
      std::vector<std::vector<Edge>> edges_per_round,
      int grid_size,
      std::vector<std::vector<int>> users // users = index(in coordinate) of {A, T1, T2, ..., B}
    );

/**
 * @brief: stage 2 local routing: get the paths(stored in 2-D vector)
 *
 *
 *
 */
std::vector<std::vector<int>> getPathsDFS(
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size
    );



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

int numTakenQubits(int x, int y, std::vector<std::vector<Node>> node_grid_per_round) {

  int result = 0;
  for(int i=0; i<4; i++) {
    if(node_grid_per_round[x][y].qubits[i].available) {
      result ++;
    }
  }

  return result;

}

std::vector<int> find2qubits_IA(int curr_r, int curr_c, std::vector<int> available_q,
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size) {

  // available_q.size = num of neighbors which should be in [3,4]
  int num_neighbor = available_q.size();
  if(num_neighbor < 3 || num_neighbor > 4) {
    std::cerr << "error: number of neighbors is wrong.\n";
    std::exit(EXIT_FAILURE);
  }

  std::vector<int> result(2, -1);

  // a vector to store the index(integer) of neighbor
  std::vector<int> neighbor;

  // a vector to store the index(coordinate) of neighbor
  std::vector<std::vector<int>> neighbor_coor;

  // get the index(integer) of current node 
  int curr = curr_r*grid_size + curr_c; 

  // get the index(integer) of the 3 neighbor with inter link
  for(int i=0; i<num_neighbor; i++) {
    neighbor.push_back(edges_per_round[curr][available_q[i]].to);
  }

  // transform neighbors' index from integer to coordinate
  for(int i=0; i<num_neighbor; i++) {
    neighbor_coor.push_back(int2coordinate(neighbor[i], grid_size));
  }

  // for each neighbor pair, calculate Dij (distance between each Ti, Tj)
  // select the 2 neighbor with the minimum Dij
  std::vector<std::vector<int>> smallest_D_neighbors{{-1, -1}}; // make it a 2-D vector cuz there may be multiple smallest neighbor
  int smallest_D = INT_MAX;
  int temp_D = 0;
  int num_D = node_grid_per_round[0][0].D.size(); // get the total number of D of a node 
  for(int i=0; i<num_neighbor-1; i++) {
    for(int j=i+1; j<num_neighbor; j++) {
      // for each pair of neighbor
      for(int k1=0; k1<num_D; k1++) {
        for(int k2=0; k2<num_D; k2++) {
          if(k2 != k1) {
            temp_D = node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[k1] +
                      node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[k1];
            if(temp_D < smallest_D) {
              smallest_D = temp_D;
              smallest_D_neighbors[0] = {i, j};
            }
            else if(temp_D = smallest_D) {
              smallest_D_neighbors.push_back({i ,j});
            }
          }
        }
      }
    }
  }

  // now smallest_D_neighbors[i] stores 2 indices of the 2 neighbor nodes to construct intra link
  // we need to transfer it to 2 indices of the 2 qubits of current node to construct intra link
  // traverse the neighbor of current nodes, if find that 2 neighbor, then get the direction and 
  // get the qubit index in that direction
  for(int i=0; i<smallest_D_neighbors.size(); i++) {

    for(int direction=0; direction<4; direction++) {
      if(edges_per_round[curr][direction].to == smallest_D_neighbors[i][0]) {
        smallest_D_neighbors[i][0] = direction;
      } 
      else if(edges_per_round[curr][direction].to == smallest_D_neighbors[i][1]) {
        smallest_D_neighbors[i][1] = direction;
      }
    }

  }

  // get the first result first
  result = smallest_D_neighbors[0];

  // then try to get the first pair in smallest_D_neighbors that can make vertical or horizontal link(if there is one)
  for(int i=0; i<smallest_D_neighbors.size(); i++) {
    // if the qubit index add up to 3, then it is vertical or horizontal
    if(smallest_D_neighbors[i][0] + smallest_D_neighbors[i][1] == 3) {
      result = smallest_D_neighbors[i];
    }

  }

  // check if result is legit
  for(int i=0; i<result.size(); i++) {
    if(result[i] < 0 || result[i] > 3) {
      std::cerr << "error: result of find2qubits_IA() is wrong.\n";
      std::exit(EXIT_FAILURE);
    }
  }


  return result;

}

/**
 * @brief: helper: dfs, recurr from one sink(set as Alice or TN) until it meets a sink as another user node(Bob) or TN 
 */
void dfs(int s, int target, int curr_q, std::vector<int>& path,
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size
    ) {

  /*
   * from the source node(s), keep recurrsion until it reach target or an end node.
   * push back every node index to the path cuz there won't be joint path 
   * end condition: 1. either Node.role = target 
   *                2. or an end node is reach -> no intra link reachable 
   *                                           -> the current qubit.to = -1 
   *
   * continueous condition: 1. from the qubit and the inter link, get the next node 
   *                        2. traverse the qubits next_q of the next node, 
   *                           if next_q != 3 - curr_q, keep dfs(next, next_q)
   */

  // push back the current node index into path   
  path.push_back(s);
  
  // end condition check
  std::vector<int> s_coor = int2coordinate(s, grid_size);
  if(s == target  
        || node_grid_per_round[s_coor[0]][s_coor[1]].qubits[curr_q].to == -1) { 
    // if this node is our target or the qubit inside this node has no intra link
    // return 
    return;
  }
  // continueous condition check
  else { 
 
    // go through the intra link
    int next_q = node_grid_per_round[s_coor[0]][s_coor[1]].qubits[curr_q].to;
   
    // get the index of next node which the current qubit is connected to through inter link 
    // by edges_per_round[curr_node][direction] = next, here direction = next_q
    // (both are 0 = uppper, 1 = left, 2 = right, 3 = bottom)
    int next = edges_per_round[s][next_q].to;
    
    // before next dfs, next_q = 3 - next_q cuz for neighbor node's current qubit
    // in the next dfs, it is in reverse direction to next_q in this iteration
    next_q = 3 - next_q;
 
    // recurr
    dfs(next, target, next_q, path,
        node_grid_per_round,
        edges_per_round,
        grid_size
        );
  }

}

/**
 * @brief: stage 2 local routing: get the paths(stored in 2-D vector)
 */
std::vector<std::vector<int>> getPathsDFS(
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size,
    std::vector<std::vector<int>> users // users = index(in coordinate) of {A, T1, T2, ..., B} 
    ) {

  // a 2-d vector result to store all the paths
  std::vector<std::vector<int>> result;

  // a path vector to store the path we found
  std::vector<int> path;

  // end node of path(index in coordinate)
  int end_node = -1;

  // for each pair in users, use bfs to get available paths
  for(int i=0; i<users.size() - 1; i++) {
    for(int j=i+1; j<users.size(); j++) {
      int s = users[i][0]*grid_size + users[i][1];   
      int t = users[j][0]*grid_size + users[j][1];
      // for each pair of s, t, the dfs starts at the neighbor nodes of s 
      // and the qubit of that neighbor node connected through inter link 
      for(int direction=0; direction<4; direction++) { // traverse s's directions to find its neighbor 
        int neighbor_s = edges_per_round[s][direction].to;
        path.push_back(s); // before dfs, put A as the starting node
        int neighbor_q = 3 - direction; // This neighbor's qubit connected with s through inter link is in the reverse direction
        dfs(neighbor_s, t, neighbor_q, path,
          node_grid_per_round,
          edges_per_round,
          grid_size
          );
        if(end_node == s) {
          result.push_back(path);
        }
        path.clear();
      }
    }
  }

  return result;

}













#endif
