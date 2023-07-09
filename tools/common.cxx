#include "tools/common.hpp"
#include "networkit/include/networkit/components/ConnectedComponents.hpp"
#include <dirent.h>
#include <mutex>


template <> gbbs_graph customReader(std::string filename) {
  return gbbs::gbbs_io::read_unweighted_symmetric_graph(filename.c_str(), false,
                                                        false);
}

template <> NetworKit::Graph customReader(std::string filename) {
  auto reader = NetworKit::GraphToolBinaryReader{};
  auto g = reader.read(filename);
  g.indexEdges(true);
  return g;
}

template <> void customWriter(const gbbs_graph &g, std::string filename) {
  gbbs::gbbs_io::write_graph_to_file(filename.c_str(), g);
};

template <> void customWriter(const NetworKit::Graph &g, std::string filename) {
  auto writer = NetworKit::GraphToolBinaryWriter{};
  writer.write(g, filename);
};

gbbs_graph NX2gbbs(const NetworKit::Graph &nx_g) {
  std::vector<gbbs::gbbs_io::Edge<gbbs::empty>> kEdges{nx_g.numberOfEdges()};
  nx_g.parallelForEdges([&](NetworKit::node n1, NetworKit::node n2,
                            NetworKit::edgeid id) {
    kEdges[id] = gbbs::gbbs_io::Edge<gbbs::empty>{static_cast<gbbs::uintE>(n1),
                                                  static_cast<gbbs::uintE>(n2)};
  });
  return gbbs::gbbs_io::edge_list_to_symmetric_graph(kEdges);
};

void NX2Lemon(const NetworKit::Graph &nx_g, lemon::ListGraph &target_g) {
  target_g.reserveNode(nx_g.numberOfNodes());
  nx_g.forNodes([&](NetworKit::node n) { target_g.addNode(); });
  nx_g.forEdges([&](NetworKit::node n1, NetworKit::node n2) {
    target_g.addEdge(lemon::ListGraph::nodeFromId(n1),
                     lemon::ListGraph::nodeFromId(n2));
  });
};

NetworKit::Graph gbbs2NX(const gbbs_graph &gbbs_g) {
  NetworKit::Graph g{gbbs_g.num_vertices()};
  auto edges = gbbs::to_edge_array<gbbs::empty>(gbbs_g).to_seq();
  for (std::size_t i = 0; i < edges.size(); i++) {
    g.addEdge(static_cast<NetworKit::node>(std::get<0>(edges[i])),
              static_cast<NetworKit::node>(std::get<1>(edges[i])),
              NetworKit::defaultEdgeWeight, true);
  }
  g.indexEdges(true);
  return g;
};

NetworKit::Graph lemon2NX(const lemon::ListGraph &lemon_g) {
  NetworKit::Graph g{static_cast<std::size_t>(lemon::countNodes(lemon_g))};
  for (lemon::ListGraph::ArcIt edge{lemon_g}; edge != lemon::INVALID; ++edge) {
    auto from = lemon_g.id(lemon_g.source(edge));
    auto to = lemon_g.id(lemon_g.target(edge));
    g.addEdge(static_cast<NetworKit::node>(from),
              static_cast<NetworKit::node>(to), NetworKit::defaultEdgeWeight,
              true);
  }
  g.indexEdges(true);
  return g;
};

bool is_spanner(const NetworKit::Graph &graph, const NetworKit::Graph &tree) {
  if (graph.numberOfNodes() != tree.numberOfNodes()) {
    std::cout << "Lost nodes" << std::endl;
    return false;
  }
  if (tree.numberOfEdges() != tree.numberOfNodes() - 1) {
    std::cout << "Wrong # of edges" << std::endl;
    return false;
  }
  auto cc = NetworKit::ConnectedComponents{tree};
  cc.run();
  if (cc.numberOfComponents() != 1) {
    std::cout << "Not connected! Components: " << cc.numberOfComponents()
              << std::endl;
    return false;
  }
  return true;
};

void eulerian_travel(const NetworKit::node focus, int &step,
                     const std::vector<std::vector<NetworKit::node>> &children,
                     const std::vector<int> &tree_depth,
                     std::vector<int> &eulerian_depths,
                     std::vector<int> &eulerian_steps) {
  eulerian_depths.push_back(tree_depth[focus]);
  eulerian_steps[focus] = step;
  for (auto child : children[focus]) {
    step++;
    eulerian_travel(child, step, children, tree_depth, eulerian_depths,
                    eulerian_steps);
    eulerian_depths.push_back(tree_depth[focus]);
  }
  if (children[focus].size() != 0) {
    eulerian_steps[focus] = step;
    // eulerian_depths.push_back(tree_depth[focus]);
  }
  step++;
}
double get_averge_stretch(NetworKit::Graph &origin, NetworKit::Graph &tree) {
  if (origin.isWeighted() || tree.isWeighted()) {
    std::cout << "Can't do weighted" << std::endl;
    return std::numeric_limits<double>::infinity();
  }

  // Have a tree
  NetworKit::node sourse{0};

  std::vector<NetworKit::node> parent{};
  parent.resize(tree.numberOfNodes());
  std::vector<std::vector<NetworKit::node>> children{};
  children.resize(tree.numberOfNodes());
  std::vector<int> depth{};
  depth.resize(tree.numberOfNodes());

  std::queue<NetworKit::node> frontier{};
  frontier.push(sourse);
  frontier.push(std::numeric_limits<uint64_t>::max());

  int level_counter = 0;
  std::set<NetworKit::node> visited{};

  do {
    if (visited.size() == tree.numberOfNodes()) {
      break;
    }
    auto focus = frontier.front();
    frontier.pop();
    if (focus == std::numeric_limits<uint64_t>::max()) {
      level_counter++;
      frontier.push(std::numeric_limits<uint64_t>::max());
      continue;
    }
    visited.insert(focus);
    depth[focus] = level_counter;
    tree.forNeighborsOf(focus, [&](NetworKit::node n) {
      if (visited.find(n) == visited.end()) {
        parent[n] = focus;
        children[focus].push_back(n);
        frontier.push(n);
      }
    });
  } while (!frontier.empty());
  // Do an eulerian path
  std::vector<int> steps{};
  steps.resize(tree.numberOfNodes());
  std::vector<int> eulerian_depths{};
  int step = 0;
  eulerian_travel(sourse, step, children, depth, eulerian_depths, steps);

  // sparse table
  std::size_t h = floor(log2(eulerian_depths.size()));
  std::vector<std::vector<int>> ST{};
  ST.resize(h + 1);
  for (auto &st_round : ST) {
    st_round.resize(eulerian_depths.size());
  }
  for (std::size_t j = 0; j < eulerian_depths.size(); j++)
    ST[0][j] = eulerian_depths[j];

  // iterative dynamic programming approach
  for (std::size_t i = 1; i <= h; i++)
    for (std::size_t j = 0; j + (1 << i) <= eulerian_depths.size(); j++)
      ST[i][j] = std::min(ST[i - 1][j], ST[i - 1][j + (1 << (i - 1))]);
  
  double total_stretch = 0.0;
  std::mutex m{};
  origin.parallelForEdges([&](NetworKit::node n1, NetworKit::node n2) {
    int n1_i = steps[n1];
    int n2_i = steps[n2];
    int l = std::min(n1_i, n2_i);
    int r = std::max(n1_i, n2_i) + 1;
    int p = 31 - __builtin_clz(r - l);
    int lca = std::min(ST[p][l], ST[p][r - (1 << p)]);
    int dis = (depth[n1] + depth[n2] - 2 * lca);
    m.lock();
    total_stretch += dis;
    m.unlock();
  });
  return total_stretch / origin.numberOfEdges();
}

std::vector<std::string> get_files(std::string dir) {
  std::vector<std::string> files{};
  DIR *dpdf;
  struct dirent *epdf;
  dpdf = opendir(dir.c_str());
  if (dpdf != NULL) {
    while ((epdf = readdir(dpdf))) {
      if ((strcmp(epdf->d_name, ".") != 0) &&
          (strcmp(epdf->d_name, "..") != 0)) {
        files.push_back(epdf->d_name);
      }
    }
  }
  closedir(dpdf);
  return files;
}