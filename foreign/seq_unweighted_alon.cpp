#include "foreign/benchmark.h"

#include "networkit/include/networkit/graph/Graph.hpp"
#include "networkit/include/networkit/graph/GraphTools.hpp"
#include <cmath>
#include <iostream>
#include <limits>

using nodes = std::unordered_set<NetworKit::node>;
using edges = std::unordered_set<NetworKit::edgeid>;
// using node_map = std::unordered_map<NetworKit::node, NetworKit::node>;
using edge_map = std::unordered_map<NetworKit::edgeid, NetworKit::edgeid>;
std::tuple<nodes, edges> build_cluster(const NetworKit::Graph &g, double alpha,
                                       const nodes &seen) {
  //   std::cout << "BUILDING" << std::endl;
  nodes visited{};
  edges inner_edge{};
  nodes to_be_visited{};
  auto terminal = std::numeric_limits<NetworKit::node>::max();
  std::queue<NetworKit::node> q{};
  auto sourse = 0;
  do {
    sourse = NetworKit::GraphTools::randomNode(g);
  } while (seen.find(sourse) != seen.end());
  q.push(sourse);
  q.push(terminal);
  edges shell_edge{};
  while (true) {
    auto s = q.front();
    q.pop();
    if (s == terminal) {
      //   std::cout << inner_edge.size() << " " << visited.size() << std::endl;
      if (shell_edge.size() < alpha * inner_edge.size()) {
        break;
      } else {
        inner_edge.insert(shell_edge.begin(), shell_edge.end());
        shell_edge.clear();
        q.push(terminal);
      }
    } else {
      visited.insert(s);
      bool have_neighbor = false;
      g.forNeighborsOf(s, [&](NetworKit::node neighbor) {
        if (visited.find(neighbor) == visited.end() &&
            to_be_visited.find(neighbor) == to_be_visited.end() &&
            seen.find(neighbor) == seen.end()) {
          shell_edge.insert(g.edgeId(s, neighbor));
          q.push(neighbor);
          to_be_visited.insert(neighbor);
          have_neighbor = true;
        }
      });
      if (!have_neighbor && visited.size() == 1) {
        return std::make_tuple(std::move(visited), std::move(inner_edge));
      }
    }
  }
  if (visited.size() == 0) {
    throw std::exception{};
  }
  //   std::cout << inner_edge.size() << " " << visited.size() << std::endl;
  return std::make_tuple(visited, inner_edge);
}

std::tuple<std::vector<nodes>, std::vector<edges>>
build_partition(NetworKit::Graph g, double D) {
  const int C = g.numberOfEdges();
  const double alpha = 4 * std::log(C) / D;
  std::vector<nodes> vec_nodes{};
  std::vector<edges> vec_edges{};
  auto terminal = std::numeric_limits<NetworKit::node>::max();
  if (alpha >= 1) {
    // we just do a bfs
    auto sourse = NetworKit::GraphTools::randomNode(g);
    nodes visited{};
    nodes to_be_visited{};
    edges backbone{};
    std::queue<NetworKit::node> q{};
    q.push(sourse);
    q.push(terminal);
    const std::size_t total_nodes = g.numberOfNodes();
    while (visited.size() < total_nodes) {
      auto s = q.front();
      q.pop();
      if (s == terminal) {
        q.push(terminal);
      } else {
        visited.insert(s);
        g.forNeighborsOf(s, [&](NetworKit::node neighbor) {
          if (visited.find(neighbor) == visited.end() &&
              to_be_visited.find(neighbor) == to_be_visited.end()) {
            // not visited
            backbone.insert(g.edgeId(s, neighbor));
            q.push(neighbor);
            to_be_visited.insert(neighbor);
          }
        });
        g.removeNode(s);
      }
    }
    vec_nodes.push_back(std::move(visited));
    vec_edges.push_back(std::move(backbone));
  } else {
    // we build smaller pieces
    nodes seen{};
    while (seen.size() != g.numberOfNodes()) {
      auto cluster = build_cluster(g, alpha, seen);
      auto &ns = std::get<0>(cluster);
      seen.insert(ns.begin(), ns.end());
      auto &es = std::get<1>(cluster);
      //   std::cout << "==================" << std::endl;
      //   std::cout << ns.size() << " " << es.size() << std::endl;
      //   std::cout << "==================" << std::endl;
      //   if (ns.size() != es.size() + 1) {
      //     std::cout << "BAD" << std::endl;
      //     throw std::exception{};
      //   }
      vec_nodes.push_back(std::move(ns));
      vec_edges.push_back(std::move(es));
    }
  }
  return std::make_tuple(std::move(vec_nodes), std::move(vec_edges));
}

std::tuple<NetworKit::Graph, edge_map>
contract(const NetworKit::Graph &g, const std::vector<nodes> &vec_nodes) {
  NetworKit::Graph new_graph{vec_nodes.size(), false, false, true};
  edge_map mapping{};
  std::unordered_map<NetworKit::node, std::size_t> identity{};
  for (std::size_t i = 0; i < vec_nodes.size(); i++) {
    const auto &cluster = vec_nodes[i];
    for (auto node : cluster) {
      identity[node] = i;
    }
  }
  for (std::size_t i = 0; i < vec_nodes.size(); i++) {
    const auto &home_cluster = vec_nodes[i];
    for (auto home_node : home_cluster) {
      g.forNeighborsOf(home_node, [&](NetworKit::node to_node) {
        std::size_t to_cluster = identity[to_node];
        if (to_cluster != i && !new_graph.hasEdge(i, to_cluster)) {
          new_graph.addEdge(i, to_cluster);
          mapping[new_graph.edgeId(i, to_cluster)] =
              g.edgeId(home_node, to_node);
        }
      });
    }
  }
  return std::make_tuple(std::move(new_graph), std::move(mapping));
}

edges low_stretch_tree(const NetworKit::Graph &g) {
  // assume g is edge indexed
  const int C = g.numberOfEdges();
  const double D = std::exp(std::pow(std::log(C), 0.51)) * 1.5;
  auto clusters = build_partition(g, D);
  auto &clus_nodes = std::get<0>(clusters);
  auto &clus_edges = std::get<1>(clusters);
  auto all_edges = edges{};
  for (auto &clu_edge : clus_edges) {
    all_edges.insert(clu_edge.cbegin(), clu_edge.cend());
  }
  if (clus_nodes.size() == 1) {
    return all_edges;
  } else {
    auto contracted_info = contract(g, clus_nodes);
    auto &contracted_g = std::get<0>(contracted_info);
    auto &mapping = std::get<1>(contracted_info);
    auto super_edges = low_stretch_tree(contracted_g);
    for (auto edge : super_edges) {
      all_edges.insert(mapping[edge]);
    }
    return all_edges;
  }
}

NetworKit::Graph unweighted_alon(NetworKit::Graph &g) {
  g.indexEdges(true);
  auto kept_edges = low_stretch_tree(g);
  NetworKit::Graph new_graph{g.numberOfNodes()};
  g.forEdges(
      [&](NetworKit::node n1, NetworKit::node n2, NetworKit::edgeid e_id) {
        if (kept_edges.find(e_id) != kept_edges.end()) {
          // is kept
          new_graph.addEdge(n1, n2);
        }
      });
  return new_graph;
}

nx_main(unweighted_alon);