#include <cmath>
#include <iostream>
#include <lemon/connectivity.h>
#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <limits>
#include <list>
#include <map>
#include <stdlib.h>

#include "basics.h"

using namespace lemon;
using namespace std;

void SelectCluster(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                   const ListGraph::EdgeMap<int> &whichE, int lenE, int Jj,
                   double x, ListGraph::EdgeMap<bool> &Edges,
                   ListGraph::NodeMap<bool> &Nodes) {
  ListGraph::Node x0;
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    x0 = i;
    break;
  }
  if (countNodes(G) == 1) {
    Nodes[x0] = true;
    return;
  }

  ListGraph::NodeMap<int> dist(G);

  Bfs<ListGraph> bfs(G);
  bfs.distMap(dist);
  bfs.init();
  bfs.addSource(x0);
  bfs.start();

  ListGraph::NodeMap<long> id(G);
  NodeStructI *order = new NodeStructI[countNodes(G)];
  ListGraph::Node *inv = new ListGraph::Node[countNodes(G)];

  long j = 0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    id[i] = j;
    order[j].id = j;
    order[j].dist = dist[i];
    ListGraph::Node a = i;
    inv[j] = a;
    Nodes[i] = false;
    ++j;
  }

  int *inEij = new int[Jj];
  int *inEijLp1 = new int[Jj];
  bool *empty = new bool[Jj];

  j = 0;
  while (j < Jj) {
    empty[j] = true;
    inEij[j] = 0;
    ++j;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    Edges[e] = false;
    if (whichE[e] < Jj)
      empty[whichE[e]] = false;
  }

  qsort(order, countNodes(G), sizeof(NodeStructI), compareNSI);

  j = 0;
  while (j < countNodes(G) && order[j].dist == 0)
    ++j;
  long until = j - 1;

  while (until < countNodes(G)) {
    long newOnes = 1;
    ++until;
    int l = order[until].dist;
    while (until < countNodes(G) - 1 && order[until + 1].dist == l) {
      ++until;
      ++newOnes;
    }
    if (until == countNodes(G) - 1)
      break;

    j = 0;
    while (j < Jj) {
      inEijLp1[j] = 0;
      ++j;
    }

    j = newOnes;
    while (j > 0) {
      ListGraph::Node &node = inv[order[until - j + 1].id];
      for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
        if (G.u(e) == node) {
          if (dist[G.v(e)] < l ||
              (dist[G.v(e)] == l && id[G.u(e)] < id[G.v(e)])) {
            if (whichE[e] < Jj)
              ++inEijLp1[whichE[e]];
          }
        } else {
          if (dist[G.u(e)] < l ||
              (dist[G.u(e)] == l && id[G.v(e)] < id[G.u(e)])) {
            if (whichE[e] < Jj)
              ++inEijLp1[whichE[e]];
          }
        }
      }
      --j;
    }

    bool accepted = true;

    j = 0;
    while (j < Jj) {
      if (empty[j]) {
        ++j;
        continue;
      }
      if ((double)inEijLp1[j] > (double)inEij[j] / x)
        accepted = false;
      ++j;
    }

    if (accepted && until - newOnes > 1) {
      until -= newOnes;
      break;
    }

    j = 0;
    while (j < Jj) {
      inEij[j] += inEijLp1[j];
      ++j;
    }

    ++j;
  }

  j = 0;
  while (j <= until) {
    Nodes[inv[order[j].id]] = true;
    ListGraph::Node a = inv[order[j].id];
    ++j;
  }

  j = 0;
  while (j <= until) {
    ListGraph::Node a = inv[order[j].id];
    if (a != x0)
      Edges[bfs.predArc(a)] = true;
    ++j;
  }

  delete[] inEij;
  delete[] inEijLp1;
  delete[] empty;
  delete[] order;
  delete[] inv;
}

void AK_Decompose(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                  ListGraph &Tree, ListGraph::EdgeMap<double> &lengthTNew,
                  ListGraph::NodeMap<ListGraph::Node> &TNewId,
                  ListGraph::NodeMap<ListGraph::Node> &TOldId) {
  if (countNodes(G) <= 2) {
    ListGraph::NodeMap<ListGraph::Node> &refnewId = TNewId;
    ListGraph::NodeMap<ListGraph::Node> &refoldId = TOldId;
    UsableGraphCopy(G, Tree, length, lengthTNew, refnewId, refoldId);
    if (countEdges(Tree) > 1) {
      ListGraph::Edge minE;
      std::list<ListGraph::Edge> toDelete;
      double min = INF;
      for (ListGraph::EdgeIt e(Tree); e != INVALID; ++e) {
        if (lengthTNew[e] < min) {
          min = lengthTNew[e];
          minE = e;
        }
      }
      for (ListGraph::EdgeIt e(Tree); e != INVALID; ++e) {
        if (e != minE)
          toDelete.push_back(e);
      }
      while (!toDelete.empty()) {
        Tree.erase(toDelete.front());
        toDelete.pop_front();
      }
    }
    return;
  }

  double x, ro, mu, y;
  long N, M;

  N = countNodes(G);
  M = countEdges(G);
  x = exp(sqrt(log((double)N) * log(log((double)N))));
  ro = 3 * log((double)N) / log(x);
  mu = 9 * ro * log((double)N);
  y = x * mu;

  double min = INF, max = -1;
  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    if (length[e] < min)
      min = length[e];
    if (length[e] > max)
      max = length[e];
  }

  // create E lists
  double LogVal = log(min) / log(y);
  int below = (int)floor(LogVal);
  LogVal = log(max) / log(y);
  int above = (int)ceil(LogVal);

  int lenE = above - below + 1;

  ListGraph::EdgeMap<int> whichE(G);

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    LogVal = log(length[e]) / log(y);
    whichE[e] = (int)floor(LogVal) - below;
  }

  // Create Tree
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = Tree.addNode();
    TNewId[i] = a;
    TOldId[a] = i;
  }

  // Copy Graph
  ListGraph Gj;
  ListGraph::NodeMap<ListGraph::Node> newId(G);
  ListGraph::NodeMap<ListGraph::Node> &refnewId = newId;
  ListGraph::NodeMap<ListGraph::Node> oldIdGj(Gj);
  ListGraph::NodeMap<ListGraph::Node> &refoldId = oldIdGj;
  ListGraph::EdgeMap<ListGraph::Edge> oldEIdGj(Gj);
  ListGraph::EdgeMap<double> lengthGj(Gj);
  UsableGraphCopy(G, Gj, length, lengthGj, newId, oldIdGj, oldEIdGj);

  // Execute decomposition
  int J = 1;
  while (countEdges(Gj) > 0) {
    ListGraph::NodeMap<long> comp(Gj);
    int c = connectedComponents(Gj, comp);
    ListGraph *Gnew = new ListGraph[c];
    ListGraph::EdgeMap<double> **lengthGnew =
        new ListGraph::EdgeMap<double> *[c];
    ListGraph::NodeMap<ListGraph::Node> **oldId =
        new ListGraph::NodeMap<ListGraph::Node> *[c];
    ListGraph::EdgeMap<ListGraph::Edge> **oldEdgeId =
        new ListGraph::EdgeMap<ListGraph::Edge> *[c];
    long j = 0;
    while (j < c) {
      lengthGnew[j] = new ListGraph::EdgeMap<double>(Gnew[j]);
      oldId[j] = new ListGraph::NodeMap<ListGraph::Node>(Gnew[j]);
      oldEdgeId[j] = new ListGraph::EdgeMap<ListGraph::Edge>(Gnew[j]);
      ++j;
    }

    ListGraph::NodeMap<ListGraph::Node> newId(G);

    GraphSplit(Gj, comp, c, length, lengthGnew, Gnew, newId, oldId, oldEdgeId);

    j = 0;
    while (j < c) {
      ListGraph::NodeMap<bool> nodes(Gnew[j]);
      ListGraph::EdgeMap<bool> edges(Gnew[j]);

      ListGraph::EdgeMap<int> whichEGj(Gnew[j]);
      for (ListGraph::EdgeIt e(Gnew[j]); e != INVALID; ++e) {
        whichEGj[e] = whichE[oldEIdGj[(*oldEdgeId[j])[e]]];
      }

      SelectCluster(Gnew[j], *(lengthGnew[j]), whichEGj, lenE, J, x, edges,
                    nodes);

      for (ListGraph::EdgeIt e(Gnew[j]); e != INVALID; ++e) {
        if (edges[e]) {
          ListGraph::Edge GEdge = oldEIdGj[(*oldEdgeId[j])[e]];
          double ind = length[GEdge];
          ListGraph::Edge newE =
              Tree.addEdge(TNewId[G.u(GEdge)], TNewId[G.v(GEdge)]);
          lengthTNew[newE] = length[GEdge];
        }
      }

      ListGraph::Node x0;

      for (ListGraph::NodeIt i(Gnew[j]); i != INVALID; ++i) {
        if (nodes[i]) {
          x0 = i;
          break;
        }
      }

      for (ListGraph::NodeIt i(Gnew[j]); i != INVALID; ++i) {
        ListGraph::Node a = i;
        if (nodes[a] && a != x0) {
          Gj.contract((*oldId[j])[x0], (*oldId[j])[a], true);
        }
      }

      ++j;
    }

    j = 0;
    while (j < c) {
      delete lengthGnew[j];
      delete oldId[j];
      delete oldEdgeId[j];
      ++j;
    }

    delete[] lengthGnew;
    delete[] oldId;
    delete[] oldEdgeId;
    delete[] Gnew;
    long h = countNodes(Gj);
    long bb = countEdges(Tree);
    if (J < lenE)
      ++J;
  }
}
