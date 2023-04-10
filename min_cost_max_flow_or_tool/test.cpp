#include <cstdint>
#include <vector>
#include <utility>

#include "ortools/graph/min_cost_flow.h"
namespace operations_research {
// MinCostFlow simple interface example.
std::pair<int, std::vector<std::vector<int>>> SimpleMinCostFlowProgram(
    const std::vector<int>& start_nodes,
    const std::vector<int>& end_nodes,
    const std::vector<int>& capacities, 
    const std::vector<int>& unit_costs,
    const std::vector<int>& supplies
    ) {
  std::pair<int, std::vector<std::vector<int>>> results; // result = {status, paths_2D}
    
  // Instantiate a SimpleMinCostFlow solver.
  SimpleMinCostFlow min_cost_flow;

  // Add each arc.
  for (int i = 0; i < start_nodes.size(); ++i) {
    int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(
        start_nodes[i], end_nodes[i], capacities[i], unit_costs[i]);
    if (arc != i) LOG(FATAL) << "Internal error";
  } 
  
  // Add node supplies.
  for (int i = 0; i < supplies.size(); ++i) {
    min_cost_flow.SetNodeSupply(i, supplies[i]);
  }

  // a 2-D vector to store edges from min cost max flow solver
  std::vector<std::vector<int>> edges_MCMF;

  // Find the min cost flow.
  int status = min_cost_flow.Solve();

  if (status == MinCostFlow::OPTIMAL) {
    LOG(INFO) << "Minimum cost flow: " << min_cost_flow.OptimalCost();
    LOG(INFO) << "";
    LOG(INFO) << " Arc   Flow / Capacity  Cost";
    for (std::size_t i = 0; i < min_cost_flow.NumArcs(); ++i) {
      int cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i);
      if(cost != 0) {
        LOG(INFO) << min_cost_flow.Tail(i) << " -> " << min_cost_flow.Head(i)
                  << "  " << min_cost_flow.Flow(i) << "  / "
                  << min_cost_flow.Capacity(i) << "       " << cost;
        edges_MCMF.push_back({min_cost_flow.Tail(i), min_cost_flow.Head(i)});
      }
    }
  } else {
//    LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
//              << status;
  }
  results = std::make_pair(status, edges_MCMF);
  return results;
}

}  // namespace operations_research

int main() {
 
  std::vector<int> start_nodes = {0, 0, 1, 1, 1, 2, 2, 3, 4};
  std::vector<int> end_nodes = {1, 2, 2, 3, 4, 3, 4, 4, 2};
  std::vector<int> capacities = {15, 8, 20, 4, 10, 15, 4, 20, 5};
  std::vector<int> unit_costs = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Define an array of supplies at each node.
  std::vector<int> supplies = {20, 0, 0, 0, -20};

  std::pair<int, std::vector<std::vector<int>>> result_MCMF = operations_research::SimpleMinCostFlowProgram(start_nodes, end_nodes,
        capacities, unit_costs, supplies); 
  return 0;
}

