#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "node.hpp"
#include "edge.hpp"
#include <cmath>
#include <queue>
#include <algorithm>
#include <limits.h>
#include <climits>
#include <queue>
#include <cstring>
#include <utility>
#include <random>

/**
 *
 *
 * helper from ortool
 *
 *
 */
#include <cstdint>
#include <vector>

#include "ortools/graph/min_cost_flow.h"
namespace operations_research {
// MinCostFlow simple interface example.
std::pair<int, std::vector<std::vector<int>>> SimpleMinCostFlowProgram(
    const std::vector<int>& start_nodes,
    const std::vector<int>& end_nodes,
    const std::vector<int>& capacities,
    const std::vector<int>& unit_costs,
    const std::vector<int>& supplies
    ) {
  std::pair<int, std::vector<std::vector<int>>> results; // result = {status, paths_2D}

  // Instantiate a SimpleMinCostFlow solver.
  SimpleMinCostFlow min_cost_flow;

  // Add each arc.
  for (int i = 0; i < start_nodes.size(); ++i) {
    int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(
        start_nodes[i], end_nodes[i], capacities[i], unit_costs[i]);
    if (arc != i) LOG(FATAL) << "Internal error";
  }

  // Add node supplies.
  for (int i = 0; i < supplies.size(); ++i) {
    min_cost_flow.SetNodeSupply(i, supplies[i]);
  }

  // a 2-D vector to store edges from min cost max flow solver
  std::vector<std::vector<int>> edges_MCMF;

  // Find the min cost flow.
  int status = min_cost_flow.Solve();

  if (status == MinCostFlow::OPTIMAL) {
//    LOG(INFO) << "Minimum cost flow: " << min_cost_flow.OptimalCost();
//    LOG(INFO) << "";
//    LOG(INFO) << " Arc   Flow / Capacity  Cost";
    for (std::size_t i = 0; i < min_cost_flow.NumArcs(); ++i) {
      int cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i);
      if(cost != 0) {
//        LOG(INFO) << min_cost_flow.Tail(i) << " -> " << min_cost_flow.Head(i)
//                  << "  " << min_cost_flow.Flow(i) << "  / "
//                  << min_cost_flow.Capacity(i) << "       " << cost;
        edges_MCMF.push_back({min_cost_flow.Tail(i), min_cost_flow.Head(i)}); 
      }
    }
  } else {
//    LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
//              << status;
  }
  results = std::make_pair(status, edges_MCMF);
  return results;
}

}  // namespace operations_research



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
 * @brief: IA algorithm - dynamic version 
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
std::vector<int> find2qubits_IA_dynamic(int curr_r, int curr_c, std::vector<int> available_q, 
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size,
    Grid* grid
    );

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
    int grid_size,
    std::vector<std::vector<int>> users, // users = index(in coordinate) of {A, T1, T2, ..., B} 
    Grid* grid
    );

/**
 * returns true if there is a path from s to t
 * aslo fill parent[] to store the path
*/
bool bfs_maxflow(int V, std::vector<std::vector<int>> rGraph, int s, int t, int parent[]);

int fordFulkerson(int V, std::vector<std::vector<int>> graph, int s, int t);

/**
 * SPFA algorithm return the shortest path(store in int path[]) 
 */
void shortestPathFaster(const std::vector<std::vector<int>>& graph, int S, int V, int path[]);

/**
 * @brief: function used to stage2 global routing dynamic
 *         choose which path to contruct in global_paths 
 *         according to _prioritizedEdges
 *
 */
std::pair<std::vector<int>, int> pathToConstruct(Grid* grid, std::vector<std::vector<int>> global_paths);

/**
 * @brief: segment a 1-D vector storing the paths according to some router nodes(TN)
 */
std::vector<std::vector<int>> segmentPath(const std::vector<int>& path, const std::vector<int> TNs);

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

    // if the length of the current path is already exceed the distance of half of the perimeter of the grid while 
    // still having no results
    // we don't need to continue actually, notice that the stuck-in while loop problem is actually it is taking
    // a long detour, so we set this criteria to avoid long unecessary detour
    if(path.size() > 2*grid_size - 1) {
      break;
    }

    // this criteria can help us leave the while loop early if we have some results
    auto it = std::min_element(std::begin(result), std::end(result), min_size);
    if(result.size() > 0) { // "it" is nullptr at the beginning cuz nothing in result
      if(path.size() > it->size()) { // if the path is already longer than the result we have
                                     // no need to continue the loop
        break; // if you use continue, there could be a potential issue!!
               // let's for last round you pop path1, then you get path2 by path1 and go
               // to this round, but path2 now has longer length than shortest path,
               // use continue will skip this round, but then in one of the future rounds 
               // you will pop path1 again and then pop path2, thus stuck in the loop 
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
  std::vector<std::vector<int>> smallest_D_neighbors; // make it a 2-D vector cuz there may be multiple smallest neighbor
  int smallest_D = INT_MAX; // smallest distance among all pairs
  int smallest_D_each_pair = INT_MAX; // smallest distance between each pair
  int temp_D = 0;
  int num_D = node_grid_per_round[0][0].D.size(); // get the total number of D of a node 
  for(int i=0; i<num_neighbor-1; i++) {
    for(int j=i+1; j<num_neighbor; j++) {
      // for each pair of neighbor, get their smallest distance
      for(int k1=0; k1<num_D; k1++) {
        for(int k2=0; k2<num_D; k2++) {
          if(k2 != k1) {
            temp_D = node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[k1] +
                      node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[k2];
            if(temp_D < smallest_D_each_pair) {
              smallest_D_each_pair = temp_D;
            }
          }
        }
      }
      
      // get the smallest distance among all pairs of users
      if(smallest_D_each_pair < smallest_D) {
        smallest_D = smallest_D_each_pair;
        smallest_D_neighbors.clear();
        smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
      }
      else if(smallest_D_each_pair == smallest_D){
        smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
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

    if(next == -1) { // if you reach a end node, next is -1, then you should not pass it to next iteration
      return;
    }

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
    std::vector<std::vector<int>> users, // users = index(in coordinate) of {A, T1, T2, ..., B} 
    Grid* grid
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
        if(neighbor_s != -1) {
          path.push_back(s); // before dfs, put A as the starting node
          int neighbor_q = 3 - direction; // This neighbor's qubit connected with s through inter link is in the reverse direction
          dfs(neighbor_s, t, neighbor_q, path,
            node_grid_per_round,
            edges_per_round,
            grid_size
            );
          end_node = path.back();
          if(end_node == t) {
            result.push_back(path);
          }
          path.clear();
        }
      }
    }
  }
  
  return result;

}

/**
 * returns true if there is a path from s to t
 * aslo fill parent[] to store the path
*/
bool bfs_maxflow(int V, std::vector<std::vector<int>> rGraph, int s, int t, int parent[])
{
  // create a visited array and mark all vertices as not visited
  bool visited[V];
  std::memset(visited, 0, sizeof(visited)); // make all entries of visited 0

  // create a queue, enqueue source vertex and mark source vertex 
  // as visited
  std::queue<int> q;
  q.push(s);
  visited[s] = true;
  parent[s] = -1;

  // standard bfs loop
  while(!q.empty()) {
    int u = q.front();
    q.pop();

    for (int v=0; v<V; v++) {

      // traverse the graph to find u's unvisited adjacent nodes 
      if (visited[v] == false && rGraph[u][v] >0) {
        // if we find a connection to the sink node
        // then there is no need to be in bfs anymore 
        // we just need to set its parent and return true
        if (v == t) {
          parent[v] = u;
          return true;
        }
        q.push(v);
        parent[v] = u;
        visited[v] = true;
      }
    }
  }

  // after while loop we didn't reach t in bfs from s,
  // then return false 
  return false;
}

int fordFulkerson(int V, std::vector<std::vector<int>> graph, int s, int t)
{
  int u, v;

  // create a residual graph and fill the residual graph 
  // with given capacities in the original graph as
  // residual capacities. 
  std::vector<std::vector<int>> rGraph(V, std::vector<int>(V));
  for (u = 0; u < V; u ++) {
    for (v = 0; v < V; v ++) {
      // rGraph[u][v] indicates residual capacity of edge
      // from u to v, if rGraph[u][v] = 0, then there is
      // no edge
      rGraph[u][v] = graph[u][v];
    }
  }

  // parent array to store the path
  // which will be filled by bfs
  int parent[V];

  // initialize max_flow as 0
  int max_flow = 0;

  // augment the flow while there is path from source to sink
  while (bfs_maxflow(V, rGraph, s, t, parent)) {
    // if there is an augmenting path

    // find minimum residual capacity of the edges along the
    // path filled by bfs. i.e., find the maximum flow through
    // the path found.
    int path_flow = INT_MAX;
    for (v = t; v != s; v = parent[v]) {
      u = parent[v];
      path_flow = std::min(path_flow, rGraph[u][v]);
    }

    // update residual capacities of the edges and 
    // reverse edges along the path
    for (v = t; v != s; v = parent[v]) {
      u = parent[v];
      // becuase augmenting path can have 2 forms 
      // one is through non full forward edges
      // two is through non empty backward edges
      // my guess: the backward edges is a type of 
      // abstraction to reduce the flow of some forward edges 
      // to achieve bigger overall flow
      rGraph[u][v] -= path_flow;
      rGraph[v][u] += path_flow;
    }

    // add path flow to overall flow
    max_flow += path_flow;
  }

  // return the overall flow
  return max_flow;
}

// SPFA
void shortestPathFaster(const std::vector<std::vector<int>>& graph, int S, int V, int path[])
{
  // create array d to store shortest distance
  int d[V];

  // create boolean array to check if vertex 
  // is present in queue or not
  // bool inQueue[V+1] = {false}; // won't work, c++ won't allow initializer with variable length of array
  bool inQueue[V];
  for(int i=0; i<V; i++) {
    inQueue[i] = false;
  }
// Initialize the distance from source to other vertices as INT_MAX
  for (int i=0; i<V; i++) {
    d[i] = INT_MAX;
  }

  // Source distance to itself is 0
  d[S] = 0;

  std::queue<int> q;
  q.push(S);
  inQueue[S] = true;
  path[S] = -1;

  while (!q.empty()) {

    // take the front vertex from queue
    int u = q.front();
    q.pop();
    inQueue[u] = false;

    // relax all adjacent edges of vertex taken from the queue
    for (int i=0; i<graph[u].size(); i++) {

      if(graph[u][i] != 0) {
        // get vertex and weight of edge
        int v = i;
        int weight = graph[u][i];

        if(d[v] > d[u] + weight) {
          d[v] = d[u] + weight;
          path[v] = u;

          // check if vertex v is in Queue or not
          // if not, add it to queue
          if (!inQueue[v]) {
            // this is where SPFA differs from bellman ford,
            // for each round,
            // it only traverses the vertices that has been updated last round
            // this is what queue is used for
            q.push(v);
            inQueue[v] = true;
          }
        }
      }
    }
  }
}

/**
 * @brief: function used to stage2 global routing dynamic
 *         choose which path to contruct in global_paths 
 *         according to _prioritizedEdges
 *
 */
std::pair<std::vector<int>, int> pathToConstruct(Grid* grid, std::vector<std::vector<int>> global_paths) {
  // the path to construct(named as shortest in static global routing)
  std::vector<int> shortest;
  int shortest_index;
  int min_length;
        
  /*      
   * you need to randomly shuffle the global_paths before you find users pair to prioritized
   * to avoid over-prioritized one user pair. But let's check out the performance first.
   */     
          
  // if there are user pairs to prioritized,
  // find the corresponding paths in global_paths 
  // to construct first
  bool prioritized {false};
  if(grid->_priorityEdges.size() > 0) {
    // first, travserse the global_paths to find if there is any match in prioritized user pair
    for (int i=0; i<global_paths.size(); i++) {
      int user0 = global_paths[i][0];
      int user1 = global_paths[i][global_paths[i].size() - 1];
      std::vector<int> user_pair = {user0, user1};
      // if found in grid->_priorityEdges, get the corresponding path and construct it  
      if(std::find(grid->_priorityEdges.begin(), grid->_priorityEdges.end(), user_pair) != grid->_priorityEdges.end()) {
        shortest = global_paths[i];
        prioritized = true;
      }
    }
  }

  if(!prioritized) {
    // if no prioritized user pair
    // get the shortest path in global_paths 
    shortest = global_paths[0];
    shortest_index = 0;
    min_length = global_paths[0].size();
    for (int i = 1; i < global_paths.size(); i++) {
      if (global_paths[i].size() < min_length) {
          min_length = global_paths[i].size();
          shortest = global_paths[i];
          shortest_index = i;
        }
    }
  }
  std::pair<std::vector<int>, int> result = std::make_pair(shortest, shortest_index);

  return result;
}

std::vector<int> find2qubits_IA_dynamic(int curr_r, int curr_c, std::vector<int> available_q,
    std::vector<std::vector<Node>> node_grid_per_round,
    std::vector<std::vector<Edge>> edges_per_round,
    int grid_size,
    Grid* grid
    ) {

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


  /*
   * for dynamic local routing, instead of calculating the distance 
   * between each Ti, Tj, only calculate the distance Ti, Tj pair 
   * within the priority edges 
   */
  std::vector<std::vector<int>> smallest_D_neighbors; // make it a 2-D vector cuz there may be multiple smallest neighbor
  int smallest_D; // smallest distance among all pairs
  int smallest_D_each_pair; // smallest distance between each pair
  int temp_D;
  int num_D; // get the total number of D of a node 
  if(grid->_priorityEdges.size() > 0) {
    
    smallest_D = INT_MAX; // smallest distance among all pairs
    smallest_D_each_pair = INT_MAX; // smallest distance between each pair   
    temp_D = 0;
    num_D = node_grid_per_round[0][0].D.size(); // get the total number of D of a node(every node has the same amount) 

    // priorityUsers[i] is a pair of priority user's indexs in node.D
    std::vector<std::vector<int>> priorityUsers;  
    for(int i=0; i<grid->_priorityEdges.size(); i++) {
   
      // transfer the node indices in _priorityEdges into node role representation
      std::vector<int> user1 = int2coordinate(grid->_priorityEdges[i][0], grid_size);
      std::vector<int> user2 = int2coordinate(grid->_priorityEdges[i][1], grid_size);
      int role1 = node_grid_per_round[user1[0]][user1[1]].role; 
      int role2 = node_grid_per_round[user2[0]][user2[1]].role; 

      // traverse the D of a node to get the index (in D, shown as user role) of the user node pair to prioritized  
      int index1;
      int index2;
      for(int i=0; i<num_D; i++) {
        if(node_grid_per_round[0][0].D[i] == role1) {
          index1 = i;	
	}
	if(node_grid_per_round[0][0].D[i] == role2) {
          index2 = i;	
	}
      }
     
      priorityUsers.push_back({index1, index2});
    }

    /*
     * same as dynamic global routing, everytime before we use priorityUsers, 
     * we need to shuffle it, but let's check out the performance first to
     * see if it is worth it.
     * performance sucks...
     */
    // obtain a random seed
    std::random_device rd;
    
    // create a random number engine using the seed
    std::mt19937 eng(rd());
    
    // shuffle the vector using the random number engine
    std::shuffle(priorityUsers.begin(), priorityUsers.end(), eng); 

    // after getting the prioritity users, select the 2 neighbors that can give the shortest distance between
    // priority users  
    for(int i=0; i<num_neighbor-1; i++) {
      for(int j=i+1; j<num_neighbor; j++) {
        // for each pair of neighbor, get their smallest distance
        temp_D = std::min(
          node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[priorityUsers[0][1]] +
            node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[priorityUsers[0][0]],
          node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[priorityUsers[0][0]] +
            node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[priorityUsers[0][1]]
        );
        if(temp_D < smallest_D_each_pair) {
          smallest_D_each_pair = temp_D;
        }
        // get the smallest distance among all pairs of users
        if(smallest_D_each_pair < smallest_D) {
          smallest_D = smallest_D_each_pair;
          smallest_D_neighbors.clear();
          smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
        }
        else if(smallest_D_each_pair == smallest_D){
          smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
        }
      }
    }
  }
  else {

  // for each neighbor pair, calculate Dij (distance between each Ti, Tj)
  // select the 2 neighbor with the minimum Dij
  smallest_D = INT_MAX; // smallest distance among all pairs
  smallest_D_each_pair = INT_MAX; // smallest distance between each pair
  temp_D = 0;
  num_D = node_grid_per_round[0][0].D.size(); // get the total number of D of a node(every node has the same amount) 
  for(int i=0; i<num_neighbor-1; i++) {
    for(int j=i+1; j<num_neighbor; j++) {
      // for each pair of neighbor, get their smallest distance
      for(int k1=0; k1<num_D; k1++) {
        for(int k2=0; k2<num_D; k2++) {
          if(k2 != k1) {
            temp_D = node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[k1] +
                      node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[k2];
            if(temp_D < smallest_D_each_pair) {
              smallest_D_each_pair = temp_D;
            }
          }
        }
      }
      
      // get the smallest distance among all pairs of users
      if(smallest_D_each_pair < smallest_D) {
        smallest_D = smallest_D_each_pair;
        smallest_D_neighbors.clear();
        smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
      }
      else if(smallest_D_each_pair == smallest_D){
        smallest_D_neighbors.push_back({neighbor[i], neighbor[j]});
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
 * @brief: segment a 1-D vector storing the paths according to one router nodes(TN)
 */
std::vector<std::vector<int>> segmentPath(const std::vector<int>& path, const std::vector<int> TNs) {
  
  std::vector<std::vector<int>> results; 

  bool found_TN {false};

  // for each TN in TNs, traverse the path to see if it exists, if so 
  // mark down their indices in the path
  std::vector<int> TN_indices_in_path;
  for(int i=0; i<path.size(); i++) {
    for(int j=0; j<TNs.size(); j++) {
      if(TNs[j] == path[i]) {
        TN_indices_in_path.push_back(i);
        found_TN = true;
      } 
    }
  }
  
  if(!found_TN) {
    results.push_back(path);
  }
  else {
    // sort TN_indices
    std::sort(TN_indices_in_path.begin(), TN_indices_in_path.end());

    int start = 0;
    std::vector<int> segment;
    for (int index : TN_indices_in_path) {
        if(start != 0) {
          segment = std::vector<int>(path.begin() + start - 1, path.begin() + index + 1);
        }
        else {
          segment = std::vector<int>(path.begin() + start, path.begin() + index + 1);
        }
        results.push_back(segment);
        start = index + 1;
    }
    segment = std::vector<int>(path.begin() + start - 1, path.end());
    results.push_back(segment);
  }
  

  return results;
}

#endif

