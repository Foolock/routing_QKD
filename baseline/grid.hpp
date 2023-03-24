#include "node.hpp"
#include <vector>
#include <iostream>

class Grid {
  public:
    Grid(std::vector<int> TN_location, size_t N): node_grid{grid_size, std::vector<Node>(N)}, grid_size(N) {
      
      // At default, Alice is placed at left bottom corner 
      // Bob is placed at right upper corner
      node_grid[grid_size-1][0].role = 1; 
      node_grid[0][grid_size-1].role = 2; 

      // initialize TN location 
      node_grid[TN_location[0]][TN_location[1]].role = 3; 
    }

    // display the node grid
    void display() {
      
      std::cout << "this is the current grid: \n\n";
      for(size_t i=0; i<grid_size; i++) {
        for(size_t j=0; j<grid_size; j++) {
          std::cout << node_grid[i][j].role << " ";
        }
        std::cout << "\n\n";
      } 
    }

  private: 
    // node grid
    std::vector<std::vector<Node>> node_grid;
    
    // grid size
    size_t grid_size;

};



