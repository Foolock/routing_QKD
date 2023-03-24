#include "node.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

class Grid {
  public:

    // constructor: place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
    /*
     * Assume only one TN
     * TN_location: location of TN
     * N: size of the grid
     * P: success rate of bell state transmission in fiber channel(edge)
     * return(void):  
     */
    Grid(std::vector<int> TN_location, int N, double P): grid_size(N), A_index{grid_size-1, 0}, 
      B_index{0, grid_size-1}, T_index(TN_location), P(P) {
    
      // create a temporary object to store node_grid
      std::vector<std::vector<Node>> temp_node_grid(grid_size, std::vector<Node>(grid_size));

      temp_node_grid[A_index[0]][A_index[1]].role = 1; 
      temp_node_grid[B_index[0]][B_index[1]].role = 2; 

      // initialize TN location 
      temp_node_grid[T_index[0]][T_index[1]].role = 3; 

      // initialize distance to Alice, Bob and TN for each node
      calDistanceToABT(temp_node_grid);
     
      // check if distance results are correct
//      std::vector<int> result;
//      for(int i=0; i<N; i++) {
//        for(int j=0; j<N; j++) {
//          result = getDistanceToABT(i, j, temp_node_grid);
//          std::cout << "for node #" << (i*N+j) << ", Distance = {"
//                    << result[0] << ", " << result[1] << ", " << result[2]
//                    << "}\n";
//        }
//      }
//
      // assign temporary node grid to node grid
      node_grid = temp_node_grid;

      // initialize edges: except for side node, each node is connected to its 
      // upper, bottom, left and right neighbor
      // create a temporary object to store edges
      std::vector<std::vector<int>> temp_edges(grid_size*grid_size, std::vector<int>(4, -1));
      for(int i=0; i<grid_size; i++) {
        for(int j=0; j<grid_size; j++) {
          int curr_node = i * grid_size + j;
          // for a curr_node, add its neighbor's index to the 2nd dimension of edges[][]
          if (i > 0) addEdge(curr_node, (i-1)*N + j, 0, temp_edges); // upper neighbor
          if (j > 0) addEdge(curr_node, i*N + (j-1), 1, temp_edges); // left neighbor
          if (j < N-1) addEdge(curr_node, i*N + (j+1), 2, temp_edges); // right neighbor
          if (i < N-1) addEdge(curr_node, (i+1)*N + j, 3, temp_edges); // bottom neighbor
        }
      }
      edges = temp_edges;

    }

    // calculate distance for each node to Alice, Bob, and TN and assign  
    /*
     * A_index : location of Alice in the grid
     * B_index : location of Bob in the grid
     * T_index : location of TN in the grid
     * return(void): 
     */
    void calDistanceToABT(
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

    // helper: check if Da, Db, Dt is calculated correctly 
    /*
     * x: x index of node
     * y: y index of node
     * return(std::vector<int>): a vector of {Da, Db, Dt}
     */
    std::vector<int> getDistanceToABT(int x, int y, std::vector<std::vector<Node>> temp_node_grid) {
      std::vector<int> result(3); // result = {Da, Db, Dt}
      result[0] = temp_node_grid[x][y].Da;
      result[1] = temp_node_grid[x][y].Db;
      result[2] = temp_node_grid[x][y].Dt;
      return result;
    }

    // helper: add edge when initialize adjacent list edges
    /*
     * curr: current node
     * neighbor: neighbor node 
     * direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void):
     */
    void addEdge(int curr, int neighbor, int direction, std::vector<std::vector<int>>& temp_edges) {
      temp_edges[curr][direction] = neighbor;
      temp_edges[neighbor][3 - direction] = curr; // if neighbor is upper of curr, 
                                                            // then curr is bottom of neighbor 
    }

    // break edge: break edge in adjacent list edges
    /*
     * curr: current node
     * direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void): 
     */
    void breakEdge(int curr, int direction) {
      // to break a edge 
      // we have to break it both for curr and neighbor
      // so before the edge is broken, i.e., the index 
      // of neighbor is set to be -1, we need to get it first
      int neighbor = edges[curr][direction]; 
      edges[curr][direction] = -1;
      // if neighbor == -1, it means it has no neighbor at the beginning(side node)
      if(neighbor != -1) {
        edges[neighbor][3 - direction] = -1;
      }
    }

    // display the node grid
    void display() {
      
//      std::cout << "this is the current grid: \n\n";
//      for(int i=0; i<grid_size; i++) {
//        for(int j=0; j<grid_size; j++) {
//          std::cout << node_grid[i][j].role << " ";
//        }
//        std::cout << "\n\n";
//      }

      std::cout << "this is the current grid(with edges): \n\n";
      for(int i=0; i<grid_size; i++) {
        for(int j=0; j<grid_size; j++) {
          std::cout << node_grid[i][j].role;
          if(edges[i*grid_size+j][2] != -1) {std::cout << "--";}
          else {std::cout << "  ";}
        }
        std::cout << "\n";
        // before it print next row of grid
        // check if current row node have some connect with next row
        for(int j=0; j<grid_size; j++) {
          if(edges[i*grid_size+j][3] != -1) {std::cout << "|  ";}
          else {std::cout << "   ";}
        }
        std::cout << "\n";
      }
    }

    // stage 1: intialize inter link with a success rate = P
    // i.e., break edges with a rate = 1-P
    void stage1() {

      std::srand(time(NULL)); // seed the random number generator
      
      for(int i=0; i<grid_size; i++) {
        for(int j=0; j<grid_size; j++) {
          // to break the edge, I only traverse the right(direction = 2)
          // and the bottom(direction = 3) edges of each node, once the 
          // edge is broken, breakEdge() will break the other edge for 
          // that neighbor 
          int curr_node = i * grid_size + j;
          for(int direction = 2; direction<4; direction++) {
            double rand_num = ((double) rand() / RAND_MAX); // generate a random number between 0 and 1
            if(rand_num > P) {
              breakEdge(curr_node, direction); 
            }
          }
        }
      } 
    }

    // stage 2 global routing: continueously traverse the grid to add shortest paths to SS
    void stage2_global() {
      int Dat = std::abs(A_index[0] - T_index[0]) + std::abs(A_index[1] - T_index[1]);
      int Dab = std::abs(A_index[0] - B_index[0]) + std::abs(A_index[1] - B_index[1]);
      int Dtb = std::abs(T_index[0] - B_index[0]) + std::abs(T_index[1] - B_index[1]);
      
      // get the shortest distances among A, B, T, if equal, randomly choose on
      int shortest_distanceABT = std::min(std::min(Dat, Dab), Dtb);
      if(shortest_distanceABT == Dat) {
        // try to form a path between A and T

      }
      else if(shortest_distanceABT == Dtb) {
        // try to form a path between T and B

      }
      else { //shortest_distanceABT == Dab
        // try to form a path between A and B

      }

    }

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
    std::vector<std::vector<int>> edges;

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



