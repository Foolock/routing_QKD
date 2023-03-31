#ifndef GRID_HPP
#define GRID_HPP

#include "node.hpp"
#include "edge.hpp"
#include <vector>
#include <iostream>

// @brief: Grid class including operations in stage 1 to stage 5
class Grid {
  public:
    /** 
     * @brief: place Alice(1), Bob(2), TN(3) and initialize Da, Db, Dt for each node
     * Assume only one TN
     *
     * input:
     *  TN_location: location of TNs
     *  N: size of the grid
     *  P: success rate of bell state transmission in fiber channel(edge)
     * return(void): 
     *  initialize grid_size, node_grid_per_round, A_index, B_index, C_index, edges_per_round(inter edge adjacent list)
     *  and P(success rate of fiber channel transmission
     */
    Grid(std::vector<std::vector<int>>& TN_locations, int N, double P);

 private:
    // node grid
    // format: node_grid_per_round[row][col] stands for the node in row-1 row and col-1 col
    std::vector<std::vector<Node>> node_grid; // original node_grid
    std::vector<std::vector<Node>> node_grid_per_round; // node_grid copy per round

    // grid size
    int grid_size;

    // location of A, B, T
    // At default, Alice is placed at left bottom corner, Bob is placed at right upper corner
    std::vector<int> A_index = {grid_size-1, 0};
    std::vector<int> B_index = {0, grid_size-1};
    std::vector<std::vector<int>> T_indices;

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


#endif

