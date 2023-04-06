#include "grid.hpp"
#include "utility.hpp"
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <algorithm>
#include <utility>

/** 
 * @brief: place Alice(1), Bob(2), TNs(3,4,...) and initialize Da, Db, Dt1, Dt2, ... for each node
 */
Grid::Grid(std::vector<std::vector<int>>& TN_locations, int N, double P, double B, double D): 
  grid_size(N), 
  T_indices(TN_locations), 
  P(P),
  B(B),
  D(D)
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

  // initialize user node list
  users.push_back(A_index);
  for(int i=0; i<T_indices.size(); i++) {
    users.push_back(T_indices[i]);
  }
  users.push_back(B_index); // B should be the last cuz it is always the sink

  // initialize shared status buffer size
  // for SS[i][j], the 1st dimension of SS stands for the Ti,
  // the 2nd dimension of SS stands for Tj
  // SS[i][j] is a vector that stores the lengths of all the paths between Ti and Tj
  SS_global.resize(users.size());
  for(auto& row : SS_global) {
    row.resize(users.size());
  }
  
  SS_local.resize(users.size());
  for(auto& row : SS_local) {
    row.resize(users.size());
  }

  // also initialize size of RK, SK
  RK.resize(users.size());
  for(auto& row : RK) {
    row.resize(users.size());
  }

  SK.resize(users.size());
  for(auto& row : SK) {
    row.resize(users.size());
  }

  // check if distance results are correct
//  std::vector<int> result;
//  for(int i=0; i<N; i++) {
//    for(int j=0; j<N; j++) {
//      result = getDistanceToABT(i, j, temp_node_grid);
//      std::cout << "for node #" << (i*N+j) << ", Distance = {";
//      for(int k=0; k<result.size(); k++) {
//        std::cout << result[k] << " ";
//      }
//      std::cout << "}\n";
//    }
//  }

  // initialize edges_per_round: except for side node, each node is connected to its 
  // upper, bottom, left and right neighbor
  // create a temporary object to store edges_per_round
  std::vector<std::vector<Edge>> temp_edges(grid_size*grid_size, std::vector<Edge>(4));
  for(int row=0; row<grid_size; row++) {
    for(int col=0; col<grid_size; col++) {
      // for a curr_node, add its neighbor's index to the 2nd dimension of edges_per_round[][]
      if (row > 0) addEdge(row, col, (row-1), col, 0, temp_node_grid, temp_edges); // upper neighbor
      if (col > 0) addEdge(row, col, row, (col-1), 1, temp_node_grid, temp_edges); // left neighbor
      if (col < N-1) addEdge(row, col, row, (col+1), 2, temp_node_grid, temp_edges); // right neighbor
      if (row < N-1) addEdge(row, col, (row+1), col, 3, temp_node_grid, temp_edges); // bottom neighbor   
    }
  }

  // assign temporary objects to node_grid_per_round and edges_per_round
  node_grid_per_round = temp_node_grid;
  node_grid = node_grid_per_round;
  edges_per_round = temp_edges;
  edges = edges_per_round;


}

/**
 * @brief: display the node grid
 */
void Grid::display() {

  std::cout << "this is the current grid(with edges_per_round): \n\n";
  for(int i=0; i<grid_size; i++) {
    for(int j=0; j<grid_size; j++) {
      if(node_grid[i][j].role == 1) {
        std::cout << std::setfill('0') << "\033[32m" << std::setw(3) << i*grid_size+j << "\033[0m";
      
      }
      else if(node_grid[i][j].role == 2) {
        std::cout << std::setfill('0') << "\033[32m" << std::setw(3) << i*grid_size+j << "\033[0m";
      }
      else if(node_grid[i][j].role == 0) {
        std::cout << std::setfill('0') << std::setw(3) << i*grid_size+j;
      }
      else {
        std::cout << std::setfill('0') << "\033[31m" << std::setw(3) << i*grid_size+j << "\033[0m";
      }
      if(node_grid_per_round[i][j].direction[2]) {std::cout << "--";}
      else {std::cout << "  ";}
    }
    std::cout << "\n";
    // before it print next row of grid
    // check if current row node have some connect with next row
    for(int j=0; j<grid_size; j++) {
      if(node_grid_per_round[i][j].direction[3]) {std::cout << " |   ";}
      else {std::cout << "     ";}
    }
    std::cout << "\n";
  }

}

/**
 * @brief: helper: a function to show the qubits and intra link status of a node
 */
void Grid::displayNodeStatus() {
  

  std::cout << "checking node's qubit and intra links status:\n";

  int check = 0;
 
  std::cout << "do you need to check? yes(1), no(0)\n";

  std::cin >> check;

  while(check) {
    
    // get index 
    int index;
    std::cout << "input index: ";
    std::cin >> index;
    std::cout << "for node " << index  << ":\n";
   
    std::vector<int> coordinate = int2coordinate(index, grid_size);

    int row = coordinate[0];
    int col = coordinate[1];

    // print node's qubits status
    std::cout << "qubits status: \n";
    for (int i=0; i<4; i++) {
      if(node_grid_per_round[row][col].qubits[i].available) {
        std::cout << "y ";
      }
      else {
        std::cout << "n ";
      }
    } 
    std::cout << "\n";

    // print node's intra links status
    std::cout << "intra links status: \n";
    for (int i=0; i<4; i++) {
      std::cout << i << "(" << node_grid_per_round[row][col].qubits[i].to << ") "; 
    }
    std::cout << "\n";


    std::cout << "check(1), stop(0)\n";

    std::cin >> check;

 }
}

/**
 * @brief: display constructed networkflow Graph
 */
void Grid::displayNetworkGraph() {
  std::cout << "print networkGraph:\n";
  for(int i=0; i<_networkGraph.size(); i++) {
    if(i == 0) {
      std::cout << "     A ";
    }
    else if(i == _networkGraph.size() - 1) {
      std::cout << " B ";
    }
    else {
      std::cout << " " << "T" << i;
    }
  }
  std::cout << "\n";
  
  for(int i=0; i<_networkGraph.size(); i++) {
    if(i == 0) {
      std::cout << "A ";
    }
    else if(i == _networkGraph.size() - 1) {
      std::cout << "B ";
    }
    else {
      std::cout << "T" << i;
    }
    std::cout<< "{ ";
    for(int j=0; j<_networkGraph[i].size(); j++) {
      std::cout << std::setfill(' ') << std::setw(2) << _networkGraph[i][j] << " ";
    }
    std::cout << " }\n";
  }
}

/**
 * @brief: reset edges_per_round[] and node_grid_per_round[] as its original copy
 *          also reset _priorityEdges
 */
void Grid::reset() {
  edges_per_round = edges;
  node_grid_per_round = node_grid;
}

/** 
 * @brief: helper: add edge when initialize adjacent list edges_per_round
 */
void Grid::addEdge(
  int cur_row, int cur_col,
  int neighbor_row, int neighbor_col,
  int direction,
  std::vector<std::vector<Node>>& temp_node_grid,
  std::vector<std::vector<Edge>>& temp_edges)
{
  temp_edges[cur_row * grid_size + cur_col][direction].to = neighbor_row * grid_size + neighbor_col;
  temp_edges[neighbor_row * grid_size + neighbor_col][3 - direction].to = cur_row * grid_size + cur_col;
  // if neighbor is upper of curr,
  // then curr is bottom of neighbor

  // also update direction info for nodes in temp_node_grid
  temp_node_grid[cur_row][cur_col].direction[direction] = true;
  temp_node_grid[neighbor_row][neighbor_col].direction[3-direction] = true;

  // also update qubits info for nodes in temp_node_grid
  temp_node_grid[cur_row][cur_col].qubits[direction].available = true;
  temp_node_grid[neighbor_row][neighbor_col].qubits[3-direction].available = true;

}

/**
 * @brief: break edge: break edge in adjacent list edges_per_round
 */
void Grid::breakEdge(int curr_row, int curr_col, int direction) {
  
  // to break a edge 
  // we have to break it both for curr and neighbor
  // so before the edge is broken, i.e., the index 
  // of neighbor is set to be -1, we need to get it first
  int neighbor = edges_per_round[curr_row * grid_size + curr_col][direction].to;
  if(neighbor != -1) {
    // if neighbor == -1, it means it has no neighbor at the beginning(side node) 
    edges_per_round[curr_row * grid_size + curr_col][direction].to = -1;
    edges_per_round[neighbor][3 - direction].to = -1;

    // also we need to update the direction info of nodes in node_grid_per_round 
    std::vector<int> neighbor_coordinate = int2coordinate(neighbor, grid_size);
    node_grid_per_round[curr_row][curr_col].direction[direction] = false;
    node_grid_per_round[neighbor_coordinate[0]][neighbor_coordinate[1]].direction[3-direction] = false;
    
    // also we need to update the qubits info for nodes in node_grid_per_round
    node_grid_per_round[curr_row][curr_col].qubits[direction].available = false;
    node_grid_per_round[neighbor_coordinate[0]][neighbor_coordinate[1]].qubits[3-direction].available = false;
  }
}

/**
 * @brief: add intra link 
 */
void Grid::addIntraEdge(int x, int y, int q1, int q2) {
  
  node_grid_per_round[x][y].qubits[q1].to = q2;
  node_grid_per_round[x][y].qubits[q2].to = q1;

}

/**
 * @brief: break intra link
 */
void Grid::breakIntraEdge(int x, int y) {
  if(node_grid_per_round[x][y].role == 0) {
    std::cerr << "error: you cannot break the intra link of a router.\n";
    std::exit(EXIT_FAILURE);
  }
  for(int q=0; q<4; q++) {
    node_grid_per_round[x][y].qubits[q].to = -1;
  }
}

/**
 * @brief: stage 1: intialize inter link with a success rate = P
 */
void Grid::stage1() {

  std::srand(time(NULL)); // seed the random number generator
  
  for(int row=0; row<grid_size; row++) {
    for(int col=0; col<grid_size; col++) {
      // to break the edge, I only traverse the right(direction = 2)
      // and the bottom(direction = 3) edges_per_round of each node, once the 
      // edge is broken, breakEdge() will break the other edge for 
      // that neighbor 
      for(int direction = 2; direction<4; direction++) {
        double rand_num = ((double) rand() / RAND_MAX); // generate a random number between 0 and 1
        if(rand_num > P) {
          breakEdge(row, col, direction); 
        }
      }
    }
  }

}

// @brief: stage 2: (global routing) create intra link with a success rate = R
void Grid::stage2_global_static() {
  
  /*
   * first, we need to find the shortest path (in hops) between any pair of nodes in {Alice, Bob and TNs}  
   */

  // create a global path pool including all the paths.
  std::vector<std::vector<int>> global_paths;

  // for each pair in users, use bfs to get available paths
  for(int i=0; i<users.size() - 1; i++) {
    for(int j=i+1; j<users.size(); j++) {
      int s = users[i][0]*grid_size + users[i][1];   
      int t = users[j][0]*grid_size + users[j][1];
      std::vector<std::vector<int>> path_st = bfs(s, t, node_grid_per_round, edges_per_round, grid_size);
      global_paths.reserve(global_paths.size() + path_st.size());
      global_paths.insert(global_paths.end(), path_st.begin(), path_st.end());
    }
  }
  
  /*
   * second, if there are edges(user pair) to prioritized, constructed the path 
   * among prioritized user pair first, 
   * if no, recurrsively find the shortest path and delete it from the grid
   * by marking them as visited until there is no available path
   */
  while(global_paths.size()) {
 
    // if there are user pairs to prioritized,
    // find the corresponding paths in global_paths 
    // to construct first
    if(_priorityEdges.size() > 0) {
        
    }

    // get the shortest path in global_paths 
    std::vector<int> shortest = global_paths[0];
    int shortest_index = 0;
    int min_length = global_paths[0].size();
    for (int i = 1; i < global_paths.size(); i++) {
      if (global_paths[i].size() < min_length) {
          min_length = global_paths[i].size();
          shortest = global_paths[i];
          shortest_index = i;
        }
    }

    // get the role of source and sink node in the shortest path
    int head = shortest[0];
    int tail = shortest[shortest.size() - 1];
    std::vector<int> head_coor = int2coordinate(head, grid_size);
    std::vector<int> tail_coor = int2coordinate(tail, grid_size);
    int head_role = node_grid_per_round[head_coor[0]][head_coor[1]].role; 
    int tail_role = node_grid_per_round[tail_coor[0]][tail_coor[1]].role; 

    // before we mark edges_per_round along the shortest path visited
    // we need to traverse the shortest path to see if there is any edge 
    // already marked as visited(in the last iteration)
    for(int i=0; i<shortest.size(); i++) {
      for(int j=0; j<4; j++) { // 0 <= direction <= 3
        if(edges_per_round[shortest[i]][j].to == shortest[i+1]) {
        if(edges_per_round[shortest[i]][j].visited == true) {
          // if this path has visited edges
          // i.e., we have disjoint paths
          // skip it
          goto next_shortest;
        }
        }
      }
    }

    // mark the edges_per_round along the shortest path visited
    // the path stores node's index(integer), so we need to find the edge first
    for(int i=0; i<shortest.size() - 1; i++) {
      // edges_per_round[shortest[i]][direction].visited = true
      // we have shortest[i], what is the direction?
      // it is the one when edges_per_round[shortest[i]][direction].to = shortest[i+1] 
      for(int j=0; j<4; j++) { // 0 <= direction <= 3
        if(edges_per_round[shortest[i]][j].to == shortest[i+1]) {
          edges_per_round[shortest[i]][j].visited = true;
          edges_per_round[shortest[i+1]][3-j].visited = true; // we need to mark for its neighbor too
        }  
      }
    }

    // print this path before we erase it
    /*
    std::cout << "found shortest path: \n";
    for(int i=0; i<shortest.size(); i++) {
      std::cout << shortest[i] << " -- ";    
    }
    std::cout << "\n\n";
    */
    
    // put shortest path length into the corresponding SS before we erase it
    // according to the role of source and sink node
    // role: A(1), B(2), T1(3), T2(4), ....
    // To fit into SS index, role need to -1
    // also shortest.size() - 1 cuz it counts the number of nodes
    SS_global[head_role-1][tail_role-1].push_back(shortest.size()-1);

    next_shortest:
      // delete this path from global path
      global_paths.erase(global_paths.begin() + shortest_index);

  }
}

// @brief: dynamic version of stage 2: (global routing) create intra link with a success rate = R
void Grid::stage2_global_dynamic() {
  
  /*
   * first, we need to find the shortest path (in hops) between any pair of nodes in {Alice, Bob and TNs}  
   */

  // create a global path pool including all the paths.
  std::vector<std::vector<int>> global_paths;

  // for each pair in users, use bfs to get available paths
  for(int i=0; i<users.size() - 1; i++) {
    for(int j=i+1; j<users.size(); j++) {
      int s = users[i][0]*grid_size + users[i][1];   
      int t = users[j][0]*grid_size + users[j][1];
      std::vector<std::vector<int>> path_st = bfs(s, t, node_grid_per_round, edges_per_round, grid_size);
      global_paths.reserve(global_paths.size() + path_st.size());
      global_paths.insert(global_paths.end(), path_st.begin(), path_st.end());
    }
  }
  
  /*
   * second, if there are edges(user pair) to prioritized, constructed the path 
   * among prioritized user pair first, 
   * if no, recurrsively find the shortest path and delete it from the grid
   * by marking them as visited until there is no available path
   */
  while(global_paths.size()) {

    // the path to construct(named as shortest in static global routing)
    std::vector<int> shortest;
    int shortest_index;
    std::pair<std::vector<int>, int> result;
    result = pathToConstruct(this, global_paths); 

    shortest = result.first;
    shortest_index = result.second;

    // get the role of source and sink node in the shortest path
    int head = shortest[0];
    int tail = shortest[shortest.size() - 1];
    std::vector<int> head_coor = int2coordinate(head, grid_size);
    std::vector<int> tail_coor = int2coordinate(tail, grid_size);
    int head_role = node_grid_per_round[head_coor[0]][head_coor[1]].role; 
    int tail_role = node_grid_per_round[tail_coor[0]][tail_coor[1]].role; 

    // before we mark edges_per_round along the shortest path visited
    // we need to traverse the shortest path to see if there is any edge 
    // already marked as visited(in the last iteration)
    for(int i=0; i<shortest.size(); i++) {
      for(int j=0; j<4; j++) { // 0 <= direction <= 3
        if(edges_per_round[shortest[i]][j].to == shortest[i+1]) {
        if(edges_per_round[shortest[i]][j].visited == true) {
          // if this path has visited edges
          // i.e., we have disjoint paths
          // skip it
          goto next_shortest;
        }
        }
      }
    }

    // mark the edges_per_round along the shortest path visited
    // the path stores node's index(integer), so we need to find the edge first
    for(int i=0; i<shortest.size() - 1; i++) {
      // edges_per_round[shortest[i]][direction].visited = true
      // we have shortest[i], what is the direction?
      // it is the one when edges_per_round[shortest[i]][direction].to = shortest[i+1] 
      for(int j=0; j<4; j++) { // 0 <= direction <= 3
        if(edges_per_round[shortest[i]][j].to == shortest[i+1]) {
          edges_per_round[shortest[i]][j].visited = true;
          edges_per_round[shortest[i+1]][3-j].visited = true; // we need to mark for its neighbor too
        }  
      }
    }

    // print this path before we erase it
    /*
    std::cout << "found shortest path: \n";
    for(int i=0; i<shortest.size(); i++) {
      std::cout << shortest[i] << " -- ";    
    }
    std::cout << "\n\n";
    */
    
    // put shortest path length into the corresponding SS before we erase it
    // according to the role of source and sink node
    // role: A(1), B(2), T1(3), T2(4), ....
    // To fit into SS index, role need to -1
    // also shortest.size() - 1 cuz it counts the number of nodes
    SS_global[head_role-1][tail_role-1].push_back(shortest.size()-1);

    next_shortest:
      // delete this path from global path
      global_paths.erase(global_paths.begin() + shortest_index);

  }
}


/**
 * @brief: stage 2: (local routing: IA algorithm) create intra link with a success rate = R
 * according to the node grid from stage 1
 */
void Grid::stage2_local_IA() {

  /*
   * first, we need to traverse the graph and calculate how many qubits are available for each node
   */
  
  // create a vector to store which qubits are available in qubits of the node
  std::vector<int> available_q;

  std::vector<int> connect_q;

  for(int row=0; row<grid_size; row++) {
    for(int col=0; col<grid_size; col++) {
      // for each node in the node_grid

      // empty available qubits vector first(it may have result from last itertaion
      available_q.clear();

      // traverse the qubits of node get available qubits
      for(int i=0; i<4; i++) {
        if(node_grid_per_round[row][col].qubits[i].available) {
          available_q.push_back(i);
        }
      }

      // if num == 2, connect the 2 available qubits 
      if(available_q.size() == 2) {
        // connect intra link between these 2 qubits
        addIntraEdge(row, col, available_q[0], available_q[1]);
      }
      // if num == 3, connect the 2 qubits with min D
      else if(available_q.size() == 3) {
        // get the 2 qubits that needs to connect
        connect_q = find2qubits_IA(row, col, available_q, node_grid_per_round, edges_per_round, grid_size);
           
        // connect intra link between these 2 qubits
        addIntraEdge(row, col, connect_q[0], connect_q[1]);
      }
      else if(available_q.size() == 4) {
        // get the 2 qubits that needs to connect
        connect_q = find2qubits_IA(row, col, available_q, node_grid_per_round, edges_per_round, grid_size);

        // connect intra link between these 2 qubits
        addIntraEdge(row, col, connect_q[0], connect_q[1]);
      
        // also connect the 2 remaining qubits
        std::vector<int> remain_q;
        std::sort(available_q.begin(), available_q.end());
        std::sort(connect_q.begin(), connect_q.end());
        remain_q.reserve(available_q.size());
        std::set_difference(available_q.begin(), available_q.end(), connect_q.begin(), connect_q.end(), std::back_inserter(remain_q));

        addIntraEdge(row, col, remain_q[0], remain_q[1]);
        
      }
    }
  }

  /*
   * second, find paths through the connect inter and intra links
   */

  // before you find paths
  // you need to break all the intra links of A, B, Ts in case in dfs you will stuck in loop
  breakIntraEdge(A_index[0], A_index[1]);
  breakIntraEdge(B_index[0], B_index[1]);
  for(int i=0; i<T_indices.size(); i++) {
    breakIntraEdge(T_indices[i][0], T_indices[i][1]);    
  }

  // to find the path, we can use dfs cuz the links have been fixed, i.e., we cannot make choices
  std::vector<std::vector<int>> paths = getPathsDFS(
      node_grid_per_round,
      edges_per_round,
      grid_size,
      users, // users = index(in coordinate) of {A, T1, T2, ..., B} 
      this
      ); 

  /*
   * finally, push back the path length to the corresponding SS
   */
  for(int i=0; i<paths.size(); i++) {
    if(paths[i].size() != 0) {
      
      // print the path
      /*
      std::cout << "found paths: \n";
      for(int j=0; j<paths[i].size(); j++) {
        std::cout << paths[i][j] << " -- "; 
      }
      std::cout << "\n\n";
      */

      int head = paths[i][0];
      int tail = paths[i][paths[i].size() - 1];
      std::vector<int> head_coor = int2coordinate(head, grid_size);
      std::vector<int> tail_coor = int2coordinate(tail, grid_size);
      int head_role = node_grid_per_round[head_coor[0]][head_coor[1]].role;
      int tail_role = node_grid_per_round[tail_coor[0]][tail_coor[1]].role; 
      // role: A(1), B(2), T1(3), T2(4), ....
      // To fit into SS index, role need to -1
      // also shortest.size() - 1 cuz it counts the number of nodes
      SS_local[head_role-1][tail_role-1].push_back(paths[i].size()-1);
    }
  }
}


/**
 * @brief: construct network flow graph and get the maximum flow value
 */
int Grid::getMaxFlow(std::vector<std::vector<std::vector<int>>> SS) {

  int result;

  /* 
   * first, transform SS to RK to SK
   */
  
  // for a path of size N, the success rate of it being able to transform to 1 bit in RK is B^(N-1)
  // for a path of size N, the success rate of it being able to transform to 1 bit in SK from RK is (1-D)^(N)
  std::srand(time(NULL)); // seed the random number generator
  double S1 = 0.0; // success rate to add a bit to RK, S1 = B^(N-1)
  double S2 = 0.0; // success rate to add a bit from RK to SK, S2 = (1-D)^(N)
  for(int i=0; i<SS.size(); i++) {
    for(int j=0; j<SS[i].size(); j++) {
      // for each path between Ti, Tj
      for(int p=0; p<SS[i][j].size(); p++) {
        
        double rand_num = ((double) rand() / RAND_MAX); // generate a random number between 0 and 1
       
        S1 = std::pow(B, SS[i][j][p] - 1);
        S2 = std::pow(1-D, SS[i][j][p]); 
        if(rand_num < S1) {
          RK[i][j] ++;
          if(rand_num < S2) {
            SK[i][j] ++;
          }
        }
      }
    }
  }

  /*
   * second, construct a network flow graph from SK
   */
  std::vector<std::vector<int>> networkGraph(users.size(), std::vector<int>(users.size()));
  for(int i=0; i<users.size(); i++) {
    for(int j=0; j<users.size(); j++) {
      networkGraph[i][j] = SK[i][j];
    }
  }

  // store in class member _networkGraph for future use
  _networkGraph = networkGraph;

  /*
   * finally, get maximum flow result  
   */
  result = fordFulkerson(users.size(), networkGraph, 0, 1);

  _key_num = result;

  return result;
}

/**
 * @brief: get a set of user pair(Ti, Tj) to prioritize from the network flow graph constructed in getMaxFlow()
 */
void Grid::getPriorityEdge() {

  // before you get _priorityEdges, you need to reset it 
  _priorityEdges.clear();

  // put row B at the end of the networkGraph 
  for(int i=0; i<_networkGraph.size(); i++) {
    std::rotate(_networkGraph[i].begin() + 1, _networkGraph[i].begin() + 2, _networkGraph[i].end());
  }
  std::rotate(_networkGraph.begin() + 1, _networkGraph.begin() + 2, _networkGraph.end());

  // display networkGraph
  displayNetworkGraph();  

  std::vector<std::vector<int>> _priorityEdges;
  
  /**
   * if _key_num = 0, just return???
   */
  if(_key_num == 0) {
    return;
  }

  /*
   * first, find a edge in _networkGraph with the maximum capacity
   * if multiple, find the last one
   */
  int max_capacity = 0;
  std::vector<int> max_edge(2, -1); // edge with max capacity 
  for(int i=0; i<users.size(); i++) {
    for(int j=0; j<users.size(); j++) {
      if(_networkGraph[i][j] >= max_capacity) {
        max_capacity = _networkGraph[i][j];
        max_edge = {i, j}; // j > i, cuz in the graph all edges goes in A->B direction
      }  
    }
  }

  std::cout << "emax in networkGraph = {" << max_edge[0] << ", " << max_edge[1] << "}\n";

  /*
   * second, find shortest path(cost = capacity) in _networkGraph
   * from Alice to Bob that contains max_edge
   * i.e., find 2 shortest paths, one from A to max_edge[0], one
   * from max_edge[1] to Bob
   */
  // use sfpa algorithm to find shortest path, here path[i] is the parent of i
  int path[users.size()];
  std::fill(path, path+users.size(), -1);
  // find the shortest path from A to max_edge[0] 
  shortestPathFaster(_networkGraph, 0, users.size(), path);
  // get the edges along A to max_edge[0]
  int curr_node = max_edge[0];
  if(curr_node != 0) { // if max_edge[0] = 0(A), it means A is the start node of the max edge
                         // then no need to find edges from A
    while(1) {
      if(path[curr_node] == -1) { // path[i] = -1 means node path[i] has no connection to node i
        break;
      }
      _priorityEdges.push_back({path[curr_node], curr_node});
      if(path[curr_node] == 0) {
        break;
      }
      else{
        curr_node = path[curr_node]; 
      }
    }
  }

  // reset the path array
  std::fill(path, path+users.size(), -1);
  // find the shortest path from max_edge[1] to B 
  shortestPathFaster(_networkGraph, max_edge[1], users.size(), path);
  // get the edges along max_edge[1] to B
  curr_node = users.size()-1;
  if(curr_node != max_edge[1]) { // if max_edge[1] is 1(B), it means B is the end node of the max edge
                         // then no need to find edges from max_edge[1] to B
    while(1) {
      if(path[curr_node] == -1) {
        break;
      }
      _priorityEdges.push_back({path[curr_node], curr_node});
      if(path[curr_node] == max_edge[1]) {
        break;
      }
      else{
        curr_node = path[curr_node]; 
      }
    }
  }

  // transfer _priorityEdges index to its actual user node index
  if(_priorityEdges.size() > 0) {
    for(int i=0; i<_priorityEdges.size(); i++) {
      if(_priorityEdges[i][0] == 0) { // when = 0, it is A
        _priorityEdges[i][0] = A_index[0]*grid_size + A_index[1]; 
      }
      else { // otherwise, it is a T, no way it will be B 
        std::vector<int> T = T_indices[_priorityEdges[i][0] - 1]; 
        _priorityEdges[i][0] = T[0]*grid_size + T[1];   
      }
      if(_priorityEdges[i][1] == users.size()-1) { // when = users.size() - 1, it is B
        _priorityEdges[i][1] = B_index[0]*grid_size + B_index[1];
      }
      else { // otherwise, it is a T, no way it will be A
        std::vector<int> T = T_indices[_priorityEdges[i][1] - 1];
        _priorityEdges[i][1] = T[0]*grid_size + T[1];
      }

      // legit check
      if(_priorityEdges[i][0] == users.size() - 1 || _priorityEdges[i][1] == 0) {
        std::cerr << "error: _prioritized edges node index wrong.\n";
        std::exit(EXIT_FAILURE);
      }
    }
  }

  std::cout << "the prioritize edges are: \n";
  for(int i=0; i<_priorityEdges.size(); i++) {
    std::cout << "{" << _priorityEdges[i][0] << ", " << _priorityEdges[i][1] << "}\n"; 
  }

}










