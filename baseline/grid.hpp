#include "node.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <limits.h>

// @brief: inter edge 
struct Edge {

  int to = -1;

  /*
   * inter edge is implemented in this way in case I need to add some properties later
   */

  // used in stage 2 global routing 
  // if an edge is visited, visited = true 
  bool visited = false;
};

// @brief: Grid class including operations in stage 1 to stage 5
class Grid { 
  public:
    
    /** 
     * @brief: place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
     * Assume only one TN
     *
     * input:
     *  TN_location: location of TN
     *  N: size of the grid
     *  P: success rate of bell state transmission in fiber channel(edge)
     * return(void): 
     *  initialize grid_size, node_grid_per_round, A_index, B_index, C_index, edges_per_round(inter edge adjacent list)
     *  and P(success rate of fiber channel transmission
     */
    Grid(std::vector<int> TN_location, int N, double P);

    /** 
     * @brief: calculate distance for each node to Alice, Bob, and TN and assign  
     *
     * return(void):
     *  assign Da, Db, Dt for each node in the node grid
     */
    void calDistanceToABT(std::vector<std::vector<Node>>& temp_node_grid);

    /**
     * @brief: helper: check if Da, Db, Dt is calculated correctly 
     *
     * input:
     *  x: x index of node
     *  y: y index of node
     * return(std::vector<int>): 
     *  a vector of {Da, Db, Dt}
     */
    std::vector<int> getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid); 

    /**
     * @brief: display the node grid
     *
     */
    void display();

    /** 
     * @brief: helper: add edge when initialize adjacent list edges_per_round
     *
     * input:
     *  curr_row, curr_col: current node's row and col index
     *  neighbor_row, neighbor_col: neighbor node's row and col index 
     *  direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void):
     *  added edges_per_round are stored in temp_edges, also node that are connected with edges_per_round are updated(its direction vector)
     *  in the temp_node_grid
     */
    void addEdge(
      int cur_row, int cur_col,
      int neighbor_row, int neighbor_col,
      int direction,
      std::vector<std::vector<Node>>& temp_node_grid,
      std::vector<std::vector<Edge>>& temp_edges);

    /**
     * @brief: break edge: break edge in adjacent list edges_per_round
     *
     * curr_row, curr_col: index of current node
     * direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void): 
     *  change direction in edges_per_round to -1 for broken nodes, change direction of nodes in node_grid_per_round to false
     */
    void breakEdge(int curr_row, int curr_col, int direction);

    /**
     * @brief: helper: transfer a node index from a integer form to a (x, y) form 
     *
     * input:
     *  node's index (int)
     * return(std::vector<int>):
     *  node's index : (x, y)
     */
    std::vector<int> int2coordinate(int node_index);

    /**
     * @brief: stage 1: intialize inter link with a success rate = P
     *  i.e., break edges_per_round with a rate = 1-P
     */
    void stage1();

    /**
     * @brief: stage 2: (global routing) create intra link with a success rate = R
     * according to the node grid from stage 1
     *
     * global routing search for shortest path for each pair A->B, A->T, T->B
     * then construct the intra link along that path. then record it to SS. then delete it
     * 
     * it will keep traversing path until there is no path available between all user.
     */
    void stage2_global();

    /**
     * @brief(need more case test): helper: bfs, traverse the current node_grid_per_round. Find available path between s, t
     * https://www.geeksforgeeks.org/print-paths-given-source-destination-using-bfs/
     *
     * input:
     *  s, t: index(integer) of source and sink node. Notice: not the role of node
     * 
     */
    std::vector<std::vector<int>> bfs(int s, int t);
 
    /**
     * @brief: helper: check if a node is in the path vector
     *
     * input:
     *  x: node index(integer)
     *  path: a path vector
     * return(bool):
     *  true: when x is in the path 
     *  false: when x is not presented in the path
     */
    bool isInPath(int x, std::vector<int> path);
 
    /**
     * @brief: reset edges_per_round[] and node_grid_per_round[] as its original copy
     *
     */
    void reset();

    void stage2_local_IA(); 

    int numTakenQubits(int x, int y); 

    void addIntraEdge(int x, int y, int q1, int q2); 

    std::vector<int> find2qubits_IA(int curr_r, int curr_c, std::vector<int> available_q);

//  private:
    // node grid
    // format: node_grid_per_round[row][col] stands for the node in row-1 row and col-1 col
    std::vector<std::vector<Node>> node_grid; // original node_grid
    std::vector<std::vector<Node>> node_grid_per_round; // node_grid copy per round

    // grid size
    int grid_size;

    // location of A, B, T
    // At default, Alice is placed at left bottom corner
    // Bob is placed at right upper corner
    std::vector<int> A_index = {grid_size-1, 0};
    std::vector<int> B_index = {0, grid_size-1};
    std::vector<int> T_index;

    // adjacent list to represent edges_per_round in the grid
    // format: edges_per_round[grid_size*grid_size][4]
    //         1st dimension: node index (integer format)
    //         2nd dimension: Edge.to = adjacent node index(integer format) in the order of upper, left, right, bottom
    // all entries are initailzed = -1 to avoid conflict with the first node(index = 0)
    std::vector<std::vector<Edge>> edges; // original initialized edges
    std::vector<std::vector<Edge>> edges_per_round; // edges copy per round

    // success rate of bell state transmission in fiber channel(edge)
    double P;

    // raw key pool (implemented as counter)
    int RKab = 0;
    int RKat = 0;
    int RKtb = 0;

    // shared state buffer
    // SSab means a shared state buffer between Alice and Bob
    // An entry in SS stores the length of a path
    // e.g., SSab[0] = 8, the first path between Alice and Bob has a length of 8
    std::vector<int> SSab;
    std::vector<int> SSat;
    std::vector<int> SStb;
};


// @brief: (constructor) place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
Grid::Grid(std::vector<int> TN_location, int N, double P): 
  grid_size(N), 
  A_index{grid_size-1, 0},
  B_index{0, grid_size-1}, 
  T_index(TN_location), 
  P(P) 
{
  // create a temporary object to store node_grid_per_round
  std::vector<std::vector<Node>> temp_node_grid(grid_size, std::vector<Node>(grid_size));      

  // place A, B, T in graph
  temp_node_grid[A_index[0]][A_index[1]].role = 1;
  temp_node_grid[B_index[0]][B_index[1]].role = 2;
  temp_node_grid[T_index[0]][T_index[1]].role = 3;

  // initialize distance to Alice, Bob and TN for each node
  calDistanceToABT(temp_node_grid);

  // check if distance results are correct
//  std::vector<int> result;
//  for(int i=0; i<N; i++) {
//    for(int j=0; j<N; j++) {
//      result = getDistanceToABT(i, j, temp_node_grid);
//      std::cout << "for node #" << (i*N+j) << ", Distance = {"
//                << result[0] << ", " << result[1] << ", " << result[2]
//                << "}\n";
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

// @brief: calculate distance for each node to Alice, Bob, and TN and assign  
void Grid::calDistanceToABT(
  std::vector<std::vector<Node>>& temp_node_grid
) {
  int Da;
  int Db;
  int Dt;

  // i -> row of grid
  // j -> column of grid
  for(int i=0; i<grid_size; i++) {
    for(int j=0; j<grid_size; j++) {

      // get distance: manhattan distance
      Da = std::abs(A_index[0] - i) + std::abs(A_index[1] - j);
      Db = std::abs(B_index[0] - i) + std::abs(B_index[1] - j);
      Dt = std::abs(T_index[0] - i) + std::abs(T_index[1] - j);

      // assign distance to the corresponding node
      temp_node_grid[i][j].D = {Da, Db, Dt};
    }
  }
}

// @brief: helper: check if Da, Db, Dt is calculated correctly 
std::vector<int> Grid::getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid) {
  std::vector<int> result(3); // result = {Da, Db, Dt}
  result = temp_node_grid[x][y].D;
  return result;
}

// @brief: display the node grid
void Grid::display() {

//  std::cout << "this is the current grid: \n\n";
//  for(int i=0; i<grid_size; i++) {
//    for(int j=0; j<grid_size; j++) {
//      std::cout << node_grid_per_round[i][j].role << " ";
//    }
//    std::cout << "\n\n";
//  }
//
  std::cout << "this is the current grid(with edges_per_round): \n\n";
  for(int i=0; i<grid_size; i++) {
    for(int j=0; j<grid_size; j++) {
      std::cout << node_grid_per_round[i][j].role;
      if(node_grid_per_round[i][j].direction[2]) {std::cout << "--";}
      else {std::cout << "  ";}
    }
    std::cout << "\n";
    // before it print next row of grid
    // check if current row node have some connect with next row
    for(int j=0; j<grid_size; j++) {
      if(node_grid_per_round[i][j].direction[3]) {std::cout << "|  ";}
      else {std::cout << "   ";}
    }
    std::cout << "\n";
  }

}


// @brief: helper: add edge when initialize adjacent list edges_per_round
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


// @brief: break edge: break edge in adjacent list edges_per_round
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
    std::vector<int> neighbor_coordinate = int2coordinate(neighbor);
    node_grid_per_round[curr_row][curr_col].direction[direction] = false;
    node_grid_per_round[neighbor_coordinate[0]][neighbor_coordinate[1]].direction[3-direction] = false;
    
    // also we need to update the qubits info for nodes in node_grid_per_round
    node_grid_per_round[curr_row][curr_col].qubits[direction].available = false;
    node_grid_per_round[neighbor_coordinate[0]][neighbor_coordinate[1]].qubits[3-direction].available = false;
  }
}

// @brief: helper: transfer a node index from a integer form to a (x, y) form 
std::vector<int> Grid::int2coordinate(int node_index) {
  std::vector<int> result(2); // result = {row, col}
  result[0] = node_index / grid_size;
  result[1] = node_index % grid_size;
  return result;
}

// @brief: stage 1: intialize inter link with a success rate = P
// i.e., break edges_per_round with a rate = 1-P
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
void Grid::stage2_global() {
  
  /*
   * first, we need to find the shortest path (in hops) between any pair of nodes in {Alice, Bob and TN}  
   */

  // get available paths from A to B and put it into path pool ab, paths_ab may be null  
  int s = A_index[0]*grid_size + A_index[1];
  int t = B_index[0]*grid_size + B_index[1];
  std::vector<std::vector<int>> paths_ab = bfs(s,t);

  // get available paths from A to T and put it into path pool at, paths_at may be null   
  s = A_index[0]*grid_size + A_index[1];
  t = T_index[0]*grid_size + T_index[1];
  std::vector<std::vector<int>> paths_at = bfs(s,t);

  // get available paths from T to B and put it into path pool tb, paths_tb may be null
  s = T_index[0]*grid_size + T_index[1];
  t = B_index[0]*grid_size + B_index[1];
  std::vector<std::vector<int>> paths_tb = bfs(s,t);

  // create our path pool for global routing
  std::vector<std::vector<int>> global_paths;
  global_paths.reserve(paths_ab.size() + paths_at.size() + paths_tb.size());
  global_paths.insert(global_paths.end(), paths_ab.begin(), paths_ab.end());
  global_paths.insert(global_paths.end(), paths_at.begin(), paths_at.end());
  global_paths.insert(global_paths.end(), paths_tb.begin(), paths_tb.end());

  /*
   * second, recurrsively find the shortest path and delete it from the grid
   * by marking them as visited until there is no available path
   */
  int iteration = 0;
  while(global_paths.size()) {

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
//          std::cout << "found joint path: \n";
//          std::cout << "curr: " << shortest[i] << " next: " << edges_per_round[shortest[i]][j].to << "\n";
//          for(int i=0; i<shortest.size(); i++) {
//            std::cout << shortest[i] << " -- ";
//          }
//          std::cout << "\n\n";
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
    std::cout << "found shortest path: \n";
    for(int i=0; i<shortest.size(); i++) {
      std::cout << shortest[i] << " -- ";    
    }
    std::cout << "\n\n";

    // put shortest path into the corresponding SS before we erase it
    if(shortest[0] == A_index[0]*grid_size+A_index[1] &&
        shortest[shortest.size()-1] == B_index[0]*grid_size+B_index[1]) {SSab.push_back(shortest.size());}
    else if(shortest[0] == A_index[0]*grid_size+A_index[1] &&
        shortest[shortest.size()-1] == T_index[0]*grid_size+T_index[1]) {SSat.push_back(shortest.size());}
    else if(shortest[0] == T_index[0]*grid_size+T_index[1] &&
        shortest[shortest.size()-1] == B_index[0]*grid_size+B_index[1]) {SStb.push_back(shortest.size());}
    else {
      std::cerr << "stage 2 error: shortest path source and sink incorrect.\n";
      std::cerr << "source :" << shortest[0] << " sink: " << shortest[shortest.size()-1] << "\n";
      std::exit(EXIT_FAILURE);
    }

  next_shortest:
    // delete this path from global path
    global_paths.erase(global_paths.begin() + shortest_index);
    iteration ++;
  }


}

// @brief(need more case test): helper: bfs, traverse the current node_grid_per_round. Find available path between s, t
std::vector<std::vector<int>> Grid::bfs(int s, int t) {

  // helper 
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
      std::vector<int> coordinate = int2coordinate(edges_per_round[last][i].to); 
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
  
  return result; 

}

// @brief: helper: check if a node is in the path vector
bool Grid::isInPath(int x, std::vector<int> path)
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

// @brief: reset edges_per_round[] and node_grid_per_round[] as its original copy
void Grid::reset() {
  edges_per_round = edges;
  node_grid_per_round = node_grid;
} 

/**
 * @brief: stage 2: (local routing: IA algorithm) create intra link with a success rate = R
 * according to the node grid from stage 1
 * 
 * IA local routing traverse the node_grid to check their qubits
 * if only 1 qubits available, does nothing 
 * if 2 qubits available, connect this 2 qubits
 * if 3 qubits available, connect 2 qubits with IA (compare Dab, Dat, Dtb for 2 options, connect the 2 qubits 
 * with the min D among {Dab, Dat, Dtb for option 1 and Dab, Dat, Dtb for option 2}, IA means when there 
 * is tie, prefer vertical or horizontal intra link
 * if 4 qubits available, connect the 2 qubits first with IA, then connect the rest 2 
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
        connect_q = find2qubits_IA(row, col, available_q);
           
        // connect intra link between these 2 qubits
        addIntraEdge(row, col, connect_q[0], connect_q[1]);
      }
      else if(available_q.size() == 4) {
        // get the 2 qubits that needs to connect
        connect_q = find2qubits_IA(row, col, available_q);

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

  // to find the path, we can use dfs cuz the links have been fixed, i.e., we cannot make choices
  

}

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
int Grid::numTakenQubits(int x, int y) {

  int result = 0;
  for(int i=0; i<4; i++) {
    if(node_grid_per_round[x][y].qubits[i].available) {
      result ++; 
    }
  }

  return result;

}

/**
 * @brief: add intra link 
 * 
 * input:
 *  x, y: index(coordinate) of a node in the node grid
 *  q1, q2: index of 2 qubits to construct intra link.
 *  0 = uppper, 1 = left, 2 = right, 3 = bottom
 */
void Grid::addIntraEdge(int x, int y, int q1, int q2) {
  
  node_grid_per_round[x][y].qubits[q1].to = q2;
  node_grid_per_round[x][y].qubits[q2].to = q1;

}

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
std::vector<int> Grid::find2qubits_IA(int curr_r, int curr_c, std::vector<int> available_q) {
 
  std::vector<int> result(2, -1);
  
  // a vector to store the index(integer) of neighbor
  std::vector<int> neighbor;

  // a vector to store the index(coordinate) of neighbor
  std::vector<std::vector<int>> neighbor_coor;

  // a 2-D vector to store Dab, Dat, Dtb of each pair of 3~4 neighbors
  /* 
   * first dimension: stores D between {01, 02, 03, 10, 12, 13, 20, 21, 23, 30, 31, 32},
   * 01 means every Dab, Dat, Dtb combinations(totally 12) between neighbor 0 and neighbor 1
   *
   * second dimension: a vector to store Dab, Dat, Dtb of 2 neighbors(there should be 12 combinations,
   * as Dab12 may not be equal to Dab21
   * 
   * this vector is initialized with all entries = -1 for testing
   */
  std::vector<std::vector<int>> D(12, std::vector<int>(12, -1));

  // temp result for D calculation
  int temp_D;

  // number of neighbors
  int num_neighbor = available_q.size();
  
  if(num_neighbor == 3) {
    
    // when there are 3 neighbors, the first dimension of 2-d vector D only has
    // 6 entris: {01, 02, 10, 12, 20, 21}, so remove half of its first dimension,
    // same for second dimension
    D.erase(D.begin(), D.begin() + D.size()/2);
    for(int i=0; i<D.size(); i++) {
      D[i].erase(D[i].begin(), D[i].begin() + D[i].size()/2);
    }

  }

  // get the index(integer) of current node 
  int curr = curr_r*grid_size + curr_c; 
  
  // get the index(integer) of the 3 neighbor with inter link
  for(int i=0; i<num_neighbor; i++) {
    neighbor.push_back(edges_per_round[curr][available_q[i]].to);
  }

  // transform neighbors' index from integer to coordinate
  for(int i=0; i<num_neighbor; i++) {
    neighbor_coor.push_back(int2coordinate(neighbor[i]));
  }

  // calculate Dab, Dat, Dtb and add them to the corresponding location of D
  int m = 0; // first dimension of D
  int n = 0; // second dimension of D
  for(int i=0; i<num_neighbor; i++) {
    for(int j=0; j<num_neighbor; j++) {
      if(j != i) {
        n = 0;
        // calculate all Dab, Dat, Dtb between neighbor i and neighbor j
        for(int k1=0; k1<3; k1++) { // k1 means D[k1] = {Da, Db, Dt}
          for(int k2=0; k2<3; k2++) { // k2 means D[k2] = {Da, Db, Dt}
            if(k2 != k1) {
              // e.g. when i=0, j=1, k1=0, k2 = 1, it is calculating Dab01

              // calculate Dab, Dat, Dtb, 
              temp_D = node_grid_per_round[neighbor_coor[i][0]][neighbor_coor[i][1]].D[k1] + 
                         node_grid_per_round[neighbor_coor[j][0]][neighbor_coor[j][1]].D[k2];
           
              // put it to corresponding location of first dimension of D 
              D[m][n] = temp_D;
              n ++;
            }
          }
        }
        m ++; 
      }        
    }
  }

  // traverse D to make sure it is correct (by checking if there is any -1)
  for(int i=0; i<D.size(); i++) {
    for(int j=0; j<D[i].size(); j++) {
      if(D[i][j] == -1) {
        std::cerr << "error: Dab, Dat, Dtb calculation mistakes!\n";
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // traverse D to get the 2 qubits with the least D within D[i][j] 
  int minVal = INT_MAX;
  std::vector<std::vector<int>> minIndex; // {minRow, minCol}
  int minRow = -1;
  int minCol = -1;
  for(int i=0; i<D.size(); i++) {
    for(int j=0; j<D[i].size(); j++) {
      if(D[i][j] < minVal) {
        minVal = D[i][j];
        minIndex.push_back({i, j});
      }
      else if(D[i][j] == minVal) {
        minIndex.push_back({i, j});
      }
    }
  }

  // create a temporary result cuz there may be ties
  std::vector<std::vector<int>> temp_result;
  for(int i=0; i<minIndex.size(); i++) {
    // the idea of my IA implementation here is to get a set of temp results
    // traverse this result set, select the first pair that can make vertical
    // or horizontal intra link, i.e., their qubit index add up as 3
    // cuz: 0 = uppper, 1 = left, 2 = right, 3 = bottom

    minRow = minIndex[i][0];
    minCol = minIndex[i][1];

    if(num_neighbor == 3) {
      // the following mapping is based on {01, 02, 10, 12, 20, 21} 
      if(minRow == 0) {temp_result.push_back({neighbor[0], neighbor[1]});}
      else if(minRow == 1) {temp_result.push_back({neighbor[0], neighbor[2]});}
      else if(minRow == 2) {temp_result.push_back({neighbor[1], neighbor[0]});}
      else if(minRow == 3) {temp_result.push_back({neighbor[1], neighbor[2]});}
      else if(minRow == 4) {temp_result.push_back({neighbor[2], neighbor[0]});}
      else if(minRow == 5) {temp_result.push_back({neighbor[2], neighbor[1]});}
      else {
        std::cerr << "error: cannot find neighbor nodes to connect with minRow.\n";
        std::exit(EXIT_FAILURE);
      }
    }
    else if(num_neighbor == 4) {
      // the following mapping is based on {01, 02, 03, 10, 12, 13, 20, 21, 23, 30, 31, 32} 
      if(minRow == 0) {temp_result.push_back({neighbor[0], neighbor[1]});}
      else if(minRow == 1) {temp_result.push_back({neighbor[0], neighbor[2]});}
      else if(minRow == 2) {temp_result.push_back({neighbor[0], neighbor[3]});}
      else if(minRow == 3) {temp_result.push_back({neighbor[1], neighbor[0]});}
      else if(minRow == 4) {temp_result.push_back({neighbor[1], neighbor[2]});}
      else if(minRow == 5) {temp_result.push_back({neighbor[1], neighbor[3]});}
      else if(minRow == 6) {temp_result.push_back({neighbor[2], neighbor[0]});}
      else if(minRow == 7) {temp_result.push_back({neighbor[2], neighbor[1]});}
      else if(minRow == 8) {temp_result.push_back({neighbor[2], neighbor[3]});}
      else if(minRow == 9) {temp_result.push_back({neighbor[3], neighbor[0]});}
      else if(minRow == 10) {temp_result.push_back({neighbor[3], neighbor[1]});}
      else if(minRow == 11) {temp_result.push_back({neighbor[3], neighbor[2]});}
      else {
        std::cerr << "error: cannot find neighbor nodes to connect with minRow.\n";
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // now temp_result[i] stores 2 indices of the 2 neighbor nodes to construct intra link
  // we need to transfer it to 2 indices of the 2 qubits of current node to construct intra link
  // traverse the neighbor of current nodes, if find that 2 neighbor, then get the direction and 
  // get the qubit index in that direction
  for(int i=0; i<temp_result.size(); i++) {

    for(int direction=0; direction<4; direction++) {
      if(edges_per_round[curr][direction].to == temp_result[i][0]) {
        temp_result[i][0] = direction;
      } 
      else if(edges_per_round[curr][direction].to == temp_result[i][1]) {
        temp_result[i][1] = direction;
      }
    }

  }

  // now get the first pair in temp_result that can make vertical or horizontal link
  for(int i=0; i<temp_result.size(); i++) {
    // if the qubit index add up to 3, then it is vertical or horizontal
    if(temp_result[i][0] + temp_result[i][1] == 3) {
      result = temp_result[i];
    }
    // notice that there may not be vertical or horizontal link sometimes, then just select the last one
    else {
      result = temp_result[i];
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
 *
 *
 *
 *
 */












































