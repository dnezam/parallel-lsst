#include "gbbs/gbbs.h"
#include "gbbs/graph_io.h"
#include "lemon/include/lemon/list_graph.h"
#include "networkit/include/networkit/graph/Graph.hpp"
#include "networkit/include/networkit/io/GMLGraphWriter.hpp"
#include "networkit/include/networkit/io/GraphToolBinaryReader.hpp"
#include "networkit/include/networkit/io/GraphToolBinaryWriter.hpp"

using gbbs_graph = gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::empty>;

// READER
template <typename G> G customReader(std::string filename);

// Writer
template <typename G> void customWriter(const G &g, std::string filename);

// nx => gbbs & lemon

gbbs_graph NX2gbbs(const NetworKit::Graph &nx_g);
// Lemon Graph does not allow copy! Ultra bad API!
void NX2Lemon(const NetworKit::Graph &nx_g, lemon::ListGraph &target_g);

// gbbs & lemon => nx
NetworKit::Graph gbbs2NX(const gbbs_graph &gbbs_g);

NetworKit::Graph lemon2NX(const lemon::ListGraph &lemon_g);

// test
bool is_spanner(const NetworKit::Graph &graph, const NetworKit::Graph &tree);

// stretch calculation
double get_averge_stretch(NetworKit::Graph &graph, NetworKit::Graph &tree);

std::vector<std::string> get_files(std::string dir);