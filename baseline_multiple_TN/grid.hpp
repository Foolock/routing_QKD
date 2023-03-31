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
     * @brief: place Alice(1), Bob(2), TNs(3,4,...) and initialize Da, Db, Dt1, Dt2, ... for each node
     *
     * input:
     *  TN_location: location of TNs
     *  N: size of the grid
     *  P: success rate of bell state transmission in fiber channel(edge)
     * return(void): 
     *  initialize grid_size, node_grid_per_round, A_index, B_index, T_indices, edges_per_round(inter edge adjacent list)
     *  and P(success rate of fiber channel transmission
     */
    Grid(std::vector<std::vector<int>>& TN_locations, int N, double P);

    /**
     * @brief: display the node grid
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

// private:
    // node grid
    // format: node_grid_per_round[row][col] stands for the node in row-1 row and col-1 col
    std::vector<std::vector<Node>> node_grid; // original node_grid
    std::vector<std::vector<Node>> node_grid_per_round; // node_grid copy per round

    // grid size
    int grid_size;

    // location of A, B, T
    // At default, Alice and Bob is placed at diag location shown in paper 
    std::vector<int> A_index = {grid_size-2, 1};
    std::vector<int> B_index = {1, grid_size-2};
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

    // shared state buffer
    // for SS[i][j], the 1st dimension of SS stands for the Ti,
    // the 2nd dimension of SS stands for Tj
    // SS[i][j] is a vector that stores the lengths of all the paths between Ti and Tj
    std::vector<std::vector<std::vector<int>>> SS; // SS will be assigned in stage 2 when needed 
};


#endif

