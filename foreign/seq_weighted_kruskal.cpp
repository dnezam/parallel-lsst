#include "foreign/basics.h"
#include "foreign/benchmark.h"
#include "lemon/include/lemon/list_graph.h"
#include <lemon/connectivity.h>
#include <lemon/kruskal.h>
using namespace lemon;
std::tuple<bool, double, double> seq_weighted_kruskal(lemon::ListGraph &input) {
  lemon::ListGraph output;
  ListGraph::EdgeMap<double> length(input);
  for (lemon::ListGraph::ArcIt edge{input}; edge != lemon::INVALID; ++edge) {
    length[edge] = 1;
  }
  ListGraph::EdgeMap<double> lengthNew(output);

  ListGraph::NodeMap<ListGraph::Node> TNewId(input);
  ListGraph::NodeMap<ListGraph::Node> TOldId(output);

  ListGraph::Node x0;
  for (ListGraph::NodeIt i(input); i != INVALID; ++i) {
    x0 = i;
    break;
  }

  double AvgEStretch, MaxEStretch, AvgEStretchB, MaxEStretchB, AvgPStretch,
      MaxPStretch;

  ListGraph::EdgeMap<ListGraph::Edge> oldEdgeId(output);
  ListGraph::EdgeMap<bool> inmin(input);
  UsableGraphCopy(input, output, length, lengthNew, TNewId, TOldId, oldEdgeId);
  auto start_t = std::chrono::high_resolution_clock::now();
  kruskal(input, length, inmin);

  std::list<ListGraph::Edge> toDelete;
  for (ListGraph::EdgeIt e(output); e != INVALID; ++e) {
    ListGraph::Edge edg = e;
    if (!inmin[oldEdgeId[edg]])
      toDelete.push_back(edg);
  }
  while (!toDelete.empty()) {
    output.erase(toDelete.front());
    toDelete.pop_front();
  }
  auto duration = std::chrono::high_resolution_clock::now() - start_t;
  auto duration_microseconds =
      std::chrono::duration_cast<std::chrono::microseconds>(duration);

  int n_nodes = countNodes(output);
  int n_edges = countEdges(output);

  if (n_nodes != (n_edges + 1) || !connected(output)) {
    return std::make_tuple(false, -1, -1);
  }

  CalculateStretches(input, length, output, lengthNew, TNewId, AvgEStretch,
                     MaxEStretch, AvgEStretchB, MaxEStretchB, AvgPStretch,
                     MaxPStretch);

  return std::make_tuple(true, duration_microseconds.count(), AvgEStretch);
}

lemon_main(seq_weighted_kruskal);