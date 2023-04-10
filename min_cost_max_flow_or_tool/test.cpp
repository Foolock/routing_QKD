#include <cstdint>
#include <vector>
#include <utility>

#include "ortools/graph/min_cost_flow.h"
namespace operations_research {
// MinCostFlow simple interface example.
std::pair<int64_t, std::vector<std::vector<int64_t>>> SimpleMinCostFlowProgram(
    const std::vector<int64_t>& start_nodes,
    const std::vector<int64_t>& end_nodes,
    const std::vector<int64_t>& capacities, 
    const std::vector<int64_t>& unit_costs,
    const std::vector<int64_t>& supplies
    ) {
  std::pair<int64_t, std::vector<std::vector<int64_t>>> results; // result = {status, paths_2D}
    
  // Instantiate a SimpleMinCostFlow solver.
  SimpleMinCostFlow min_cost_flow;

  // Add each arc.
  for (int64_t i = 0; i < start_nodes.size(); ++i) {
    int64_t arc = min_cost_flow.AddArcWithCapacityAndUnitCost(
        start_nodes[i], end_nodes[i], capacities[i], unit_costs[i]);
    if (arc != i) LOG(FATAL) << "Internal error";
  } 
  
  // Add node supplies.
  for (int64_t i = 0; i < supplies.size(); ++i) {
    min_cost_flow.SetNodeSupply(i, supplies[i]);
  }

  // a 2-D vector to store edges from min cost max flow solver
  std::vector<std::vector<int64_t>> edges_MCMF;

  // Find the min cost flow.
  int64_t status = min_cost_flow.Solve();

  if (status == MinCostFlow::OPTIMAL) {
    LOG(INFO) << "Minimum cost flow: " << min_cost_flow.OptimalCost();
    LOG(INFO) << "";
    LOG(INFO) << " Arc   Flow / Capacity  Cost";
    for (std::size_t i = 0; i < min_cost_flow.NumArcs(); ++i) {
      int64_t cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i);
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
 
  std::vector<int64_t> start_nodes = {0, 0, 1, 1, 1, 2, 2, 3, 4};
  std::vector<int64_t> end_nodes = {1, 2, 2, 3, 4, 3, 4, 4, 2};
  std::vector<int64_t> capacities = {15, 8, 20, 4, 10, 15, 4, 20, 5};
  std::vector<int64_t> unit_costs = {0, 1, 0, 0, 0, 0, 0, 0, 0};

  // Define an array of supplies at each node.
  std::vector<int64_t> supplies = {20, 0, 0, 0, -20};

  std::pair<int64_t, std::vector<std::vector<int64_t>>> result_MCMF = operations_research::SimpleMinCostFlowProgram(start_nodes, end_nodes,
        capacities, unit_costs, supplies); 
  return 0;
}

