#include "node.hpp"
#include <vector>
#include <iostream>
#include <cmath>

class Grid {
  public:

    // constructor: place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
    /*
     * Assume only one TN
     * TN_location: location of TN
     * N: size of the grid
     * return(void):  
     */
    Grid(std::vector<int> TN_location, int N): grid_size(N) {
    
      // create a temporary object to store node_grid
      std::vector<std::vector<Node>> temp_node_grid(grid_size, std::vector<Node>(grid_size));
      
      // At default, Alice is placed at left bottom corner 
      // Bob is placed at right upper corner
      std::vector<int> A_index = {N-1, 0};
      std::vector<int> B_index = {0, N-1};
      temp_node_grid[A_index[0]][A_index[1]].role = 1; 
      temp_node_grid[B_index[0]][B_index[1]].role = 2; 

      // initialize TN location 
      temp_node_grid[TN_location[0]][TN_location[1]].role = 3; 

      // initialize distance to Alice, Bob and TN for each node
      calDistanceToABT(A_index, B_index, TN_location, temp_node_grid);
     
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
      node_grid = temp_node_grid;
    }

    // calculate distance for each node to Alice, Bob, and TN and assign  
    /*
     * A_index : location of Alice in the grid
     * B_index : location of Bob in the grid
     * T_index : location of TN in the grid
     * return(void): 
     */
    void calDistanceToABT(
      std::vector<int> A_index, 
      std::vector<int> B_index, 
      std::vector<int> T_index,
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

    // display the node grid
    void display() {
      
      std::cout << "this is the current grid: \n\n";
      for(int i=0; i<grid_size; i++) {
        for(int j=0; j<grid_size; j++) {
          std::cout << node_grid[i][j].role << " ";
        }
        std::cout << "\n\n";
      } 
    }

  private: 
    // node grid
    std::vector<std::vector<Node>> node_grid;
    
    // grid size
    int grid_size;

};



