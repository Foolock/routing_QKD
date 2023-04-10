#ifndef EDGE_HPP
#define EDGE_HPP

#include <climits>

// @brief: inter edge 
struct Edge {

  int to = -1;

  /*
   * inter edge is implemented in this way in case I need to add some properties later
   */

  // used in stage 2 global routing 
  // if an edge is visited, visited = true 
  bool visited = false;

  // cost of the edge, used in stage 2 MCMF
  int cost = INT_MAX;
};

#endif
