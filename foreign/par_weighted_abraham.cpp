#include "foreign/AbrahamBartal_parallel.h"
#include "foreign/benchmark.h"
#include "lemon/include/lemon/list_graph.h"
#include <lemon/connectivity.h>
using namespace lemon;
std::tuple<bool, double, double> par_weighted_abraham(lemon::ListGraph &input) {
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

  AbrahamBartal_par ab(input, length);
  auto start_t = std::chrono::high_resolution_clock::now();
  ab.run(x0, output, lengthNew, TNewId, TOldId);
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
lemon_main(par_weighted_abraham);