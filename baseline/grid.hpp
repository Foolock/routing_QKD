#include "node.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

// @brief: inter edge 
struct Edge {

  int to = -1;

  /*
   * inter edge is implemented in this way in case I need to add some properties later
   */
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
     *  initialize grid_size, node_grid, A_index, B_index, C_index, edges(inter edge adjacent list)
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
     * @brief: helper: add edge when initialize adjacent list edges
     *
     * input:
     *  curr_row, curr_col: current node's row and col index
     *  neighbor_row, neighbor_col: neighbor node's row and col index 
     *  direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void):
     *  added edges are stored in temp_edges, also node that are connected with edges are updated(its direction vector)
     *  in the temp_node_grid
     */
    void addEdge(
      int cur_row, int cur_col,
      int neighbor_row, int neighbor_col,
      int direction,
      std::vector<std::vector<Node>>& temp_node_grid,
      std::vector<std::vector<Edge>>& temp_edges);

    /**
     * @brief: break edge: break edge in adjacent list edges
     *
     * curr_row, curr_col: index of current node
     * direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void): 
     *  change direction in edges to -1 for broken nodes, change direction of nodes in node_grid to false
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

    // stage 1: intialize inter link with a success rate = P
    // i.e., break edges with a rate = 1-P
    void stage1();

  private:
    // node grid
    std::vector<std::vector<Node>> node_grid;

    // grid size
    int grid_size;

    // location of A, B, T
    // At default, Alice is placed at left bottom corner
    // Bob is placed at right upper corner
    std::vector<int> A_index = {grid_size-1, 0};
    std::vector<int> B_index = {0, grid_size-1};
    std::vector<int> T_index;

    // adjacent list to represent edges in the grid
    // format: edges[grid_size*grid_size][4]
    //         1st dimension: node index
    //         2nd dimension: adjacent node index in the order of upper, left, right, bottom
    // all entries are initailzed = -1 to avoid conflict with the first node(index = 0)
    std::vector<std::vector<Edge>> edges;

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
  // create a temporary object to store node_grid
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

  // initialize edges: except for side node, each node is connected to its 
  // upper, bottom, left and right neighbor
  // create a temporary object to store edges
  std::vector<std::vector<Edge>> temp_edges(grid_size*grid_size, std::vector<Edge>(4));
  for(int row=0; row<grid_size; row++) {
    for(int col=0; col<grid_size; col++) {
      // for a curr_node, add its neighbor's index to the 2nd dimension of edges[][]
      if (row > 0) addEdge(row, col, (row-1), col, 0, temp_node_grid, temp_edges); // upper neighbor
      if (col > 0) addEdge(row, col, row, (col-1), 1, temp_node_grid, temp_edges); // left neighbor
      if (col < N-1) addEdge(row, col, row, (col+1), 2, temp_node_grid, temp_edges); // right neighbor
      if (row < N-1) addEdge(row, col, (row+1), col, 3, temp_node_grid, temp_edges); // bottom neighbor   
    }
  }

  // assign temporary objects to node_grid and edges
  node_grid = temp_node_grid;
  edges = temp_edges;
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
//      std::cout << node_grid[i][j].role << " ";
//    }
//    std::cout << "\n\n";
//  }
//

  std::cout << "this is the current grid(with edges): \n\n";
  for(int i=0; i<grid_size; i++) {
    for(int j=0; j<grid_size; j++) {
      std::cout << node_grid[i][j].role;
      if(node_grid[i][j].direction[2]) {std::cout << "--";}
      else {std::cout << "  ";}
    }
    std::cout << "\n";
    // before it print next row of grid
    // check if current row node have some connect with next row
    for(int j=0; j<grid_size; j++) {
      if(node_grid[i][j].direction[3]) {std::cout << "|  ";}
      else {std::cout << "   ";}
    }
    std::cout << "\n";
  }

}


// @brief: helper: add edge when initialize adjacent list edges
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


// @brief: break edge: break edge in adjacent list edges
void Grid::breakEdge(int curr_row, int curr_col, int direction) {
  
  // to break a edge 
  // we have to break it both for curr and neighbor
  // so before the edge is broken, i.e., the index 
  // of neighbor is set to be -1, we need to get it first
  int neighbor = edges[curr_row * grid_size + curr_col][direction].to;
  if(neighbor != -1) {
  // if neighbor == -1, it means it has no neighbor at the beginning(side node) 
  edges[curr_row * grid_size + curr_col][direction].to = -1;
  edges[neighbor][3 - direction].to = -1;

  // also we need to update the direction info of nodes in node_grid 
  std::vector<int> neighbor_coordinate = int2coordinate(neighbor);
  node_grid[curr_row][curr_col].direction[direction] = false;
  node_grid[neighbor_coordinate[0]][neighbor_coordinate[1]].direction[3-direction] = false;
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
// i.e., break edges with a rate = 1-P
void Grid::stage1() {

  std::srand(time(NULL)); // seed the random number generator
  
  for(int row=0; row<grid_size; row++) {
    for(int col=0; col<grid_size; col++) {
      // to break the edge, I only traverse the right(direction = 2)
      // and the bottom(direction = 3) edges of each node, once the 
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


