// From Bradley, Hax and Maganti, 'Applied Mathematical Programming', figure 8.1
#include <cstdint>
#include <vector>

#include "ortools/graph/min_cost_flow.h"

namespace operations_research {
// MinCostFlow simple interface example.
void SimpleMinCostFlowProgram() {
  // Instantiate a SimpleMinCostFlow solver.
  SimpleMinCostFlow min_cost_flow;

  // Define four parallel arrays: sources, destinations, capacities,
  // and unit costs between each pair. For instance, the arc from node 0
  // to node 1 has a capacity of 15.
  std::vector<int> start_nodes = {0, 0, 1, 2, 1, 3, 2, 3};
  std::vector<int> end_nodes = {1, 2, 3, 3, 0, 1, 0, 2};
  std::vector<int> capacities = {2, 1, 0, 1, 1, 1, 1, 1};
  std::vector<int> unit_costs = {-1, -1, -1, -1, 1, 1, 1, 1};

  // Define an array of supplies at each node.
  std::vector<int> supplies = {1, 0, 0, -1};

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

  // Find the min cost flow.
  int status = min_cost_flow.Solve();

  if (status == MinCostFlow::OPTIMAL) {
    LOG(INFO) << "Minimum cost flow: " << min_cost_flow.OptimalCost();
    LOG(INFO) << "";
    LOG(INFO) << " Arc   Flow / Capacity  Cost";
    for (std::size_t i = 0; i < min_cost_flow.NumArcs(); ++i) {
      int cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i);
      LOG(INFO) << min_cost_flow.Tail(i) << " -> " << min_cost_flow.Head(i)
                << "  " << min_cost_flow.Flow(i) << "  / "
                << min_cost_flow.Capacity(i) << "       " << cost;
    }
  } else { 
    LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
              << status;
  }
}

}  // namespace operations_research

int main() {
  operations_research::SimpleMinCostFlowProgram();
  return EXIT_SUCCESS;
}
