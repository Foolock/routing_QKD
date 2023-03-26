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

// @brief: inter edge 
struct Edge {

  int to = -1;

  /*
   * inter edge is implemented in this way in case I need to add some properties later
   */

  // used in stage 2 global routing 
  // if an edge is visited, visited = true 
  bool visited;
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

    std::vector<std::vector<int>> bfs(int s, int t);
  
    bool isInPath(int x, std::vector<int> path);
  
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
      temp_node_grid[i][j].Da = Da;
      temp_node_grid[i][j].Db = Db;
      temp_node_grid[i][j].Dt = Dt;
    }
  }
}

// @brief: helper: check if Da, Db, Dt is calculated correctly 
std::vector<int> Grid::getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid) {
  std::vector<int> result(3); // result = {Da, Db, Dt}
  result[0] = temp_node_grid[x][y].Da;
  result[1] = temp_node_grid[x][y].Db;
  result[2] = temp_node_grid[x][y].Dt;
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

/**
 * @brief: stage 2: (global routing) create intra link with a success rate = R
 * according to the node grid from stage 1
 *
 * global routing search for shortest path for each pair A->B, A->T, T->B
 * then construct the intra link along that path. then record it to SS. then delete it
 * 
 * it will keep traversing path until there is no path available between all user.
 */
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
  while(global_paths.size()) {

    // print size of global paths
    std::cout << "the size of global paths = " << global_paths.size() << "\n";

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

    // mark the edges_per_round along the shortest path visited
    // the path stores node's index(integer), so we need to find the edge first
    for(int i=0; i<shortest.size() - 1; i++) {
      // edges_per_round[shortest[i]][direction].visited = true
      // we have shortest[i], what is the direction?
      // it is the one when edges_per_round[shortest[i]][direction].to = shortest[i+1] 
      for(int j=0; j<4; j++) { // 0 <= direction <= 3
        if(edges_per_round[shortest[i]][j].to == shortest[i+1]) {
//          if(edges_per_round[shortest[i]][j].visited == true) {goto next_shortest;}; // if this path has visited edges
//                                                                                     // i.e., we have disjoint paths
//                                                                                     // skip it
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
  }


}

/**
 * @brief(need more case test): helper: bfs, traverse the current node_grid_per_round. Find available path between s, t
 * https://www.geeksforgeeks.org/print-paths-given-source-destination-using-bfs/
 *
 * input:
 *  s, t: index(integer) of source and sink node. i.e., role of node. (1 = Alice, 3=TN, 2=Bob)
 * 
 */
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


















