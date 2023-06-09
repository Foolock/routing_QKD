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
    Grid(const std::vector<std::vector<int>>& TN_locations, int N, double P, double B, double D);

    /**
     * @brief: display the node grid
     */
    void display();

    /**
     * @brief: helper: a function to show the qubits and intra link status of a node
     */
    void displayNodeStatus();
   
    /**
     * @brief: display constructed networkflow Graph
     */
    void displayNetworkGraph();

    /**
     * @brief: reset edges_per_round[] and node_grid_per_round[] as its original copy
     */
    void reset();

    /** 
     * @brief: add inter edge when initialize adjacent list edges_per_round
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
     * @brief: break inter edge: break edge in adjacent list edges_per_round
     *
     * curr_row, curr_col: index of current node
     * direction: neighbor's direction from curr: 0 = uppper, 1 = left, 2 = right, 3 = bottom
     * return(void): 
     *  change direction in edges_per_round to -1 for broken nodes, change direction of nodes in node_grid_per_round to false
     */
    void breakEdge(int curr_row, int curr_col, int direction);

    /**
     * @brief: add intra link 
     * 
     * input:
     *  x, y: index(coordinate) of a node in the node grid
     *  q1, q2: index of 2 qubits to construct intra link.
     *  0 = uppper, 1 = left, 2 = right, 3 = bottom
     */
    void addIntraEdge(int x, int y, int q1, int q2);

    /**
     * @brief: break intra link
     * 
     * input:
     *  x, y: index(coordinate) of a node in the node grid
     */
    void breakIntraEdge(int x, int y);

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
    void stage2_global_static();

    /**
     * @brief: dynamic version of stage 2: (global routing) create intra link with a success rate = R
     *        Will construct paths between prioritized user pair first. Then try to find shortest path
     */
    void stage2_global_dynamic();

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
    void stage2_local_IA_static();
	
    /**
     * @brief: dynamic version of stage 2: (local routing: IA algorithm) create intra link with a success rate = R
     * according to the node grid from stage 1
     */
    void stage2_local_IA_dynamic();

    /**
     * @brief: apply min cost max flow in stage 2
     */
    void stage2_min_cost_max_flow();

    /**
     * @brief: construct network flow graph and get the maximum flow value
     *  1. transform SS to RK
     *  2. transform RK to SK
     *  3. based on SK(capacity), construct network flow graph
     *  4. adapt ford-fulkerson to solve a maximum flow 
     */
    int getMaxFlow(std::vector<std::vector<std::vector<int>>> SS);

    /**
     * @brief: get a set of user pair(Ti, Tj) to prioritize from the network flow graph constructed in getMaxFlow()
     */
    void getPriorityEdge(); 

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

    // create a list(vector) of users, users[i] is user index
    std::vector<std::vector<int>> users;


    // adjacent list to represent edges_per_round in the grid
    // format: edges_per_round[grid_size*grid_size][4]
    //         1st dimension: node index (integer format)
    //         2nd dimension: Edge.to = adjacent node index(integer format) in the order of upper, left, right, bottom
    // all entries are initailzed = -1 to avoid conflict with the first node(index = 0)
    std::vector<std::vector<Edge>> edges; // original initialized edges
    std::vector<std::vector<Edge>> edges_per_round; // edges copy per round

    // success rate of bell state transmission in fiber channel(inter edge)
    double P;

    // success rate of bell swap within the node(intra edge)
    double B;

    // decoherence rate of bell pair transmission in fiber channel
    double D;

    // shared state buffer, use the node's role as index
    // for SS[i][j], the 1st dimension of SS stands for the role of Ti,
    // the 2nd dimension of SS stands for the role of Tj
    // SS[i][j] is a vector that stores the lengths of all the paths between Ti and Tj
    // storing format(1st dimension): A, B, T1, T2, T3, ... 
    std::vector<std::vector<std::vector<int>>> SS_global; // SS_global will be assigned in stage 2 global  
    std::vector<std::vector<std::vector<int>>> SS_local; // SS_local will be assigned in stage 2 local 
    std::vector<std::vector<std::vector<int>>> SS_MCMF; // SS_MCMF will be assigned in stage 2 min_cost_max_flow  

    // raw key pool
    // for RK[i][j], the 1st dimension of SS stands for the Ti,
    // the 2nd dimension of SS stands for Tj
    // SS[i][j] is a vector that stores the lengths of all the paths between Ti and Tj
    // storing format(1st dimension): A, B, T1, T2, T3, ... 
    // storing format(2nd dimension): A, B, T1, T2, T3, ... 
    std::vector<std::vector<int>> RK;

    // secret key pool
    // format same as RK
    // storing format(1st dimension): A, B, T1, T2, T3, ... 
    // storing format(2nd dimension): A, B, T1, T2, T3, ... 
    std::vector<std::vector<int>> SK;

    // a network flow graph to calculate the max key value(max flow) from Alice to Bob
    // storing format(1st dimension): A, T1, T2, T3, ..., B (format finalized in getPriorityEdges()) 
    // storing format(2nd dimension): A, T1, T2, T3, ..., B (format finalized in getPriorityEdges()) 
    std::vector<std::vector<int>> _networkGraph;

    // Edges to prioritized 
    std::vector<std::vector<int>> _priorityEdges;

    // final key number from A to B
    int _key_num;
};


#endif


