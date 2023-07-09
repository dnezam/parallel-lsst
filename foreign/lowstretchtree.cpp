#include <iostream>
#include <lemon/connectivity.h>
#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <list>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "EmelElkin.h"
#include "basics.h"
#include "decompose.h"

using namespace lemon;
using namespace std;

void EmelElkinBase::GetConeDistance(const ListGraph &G,
                                    const ListGraph::EdgeMap<double> &length,
                                    const ListGraph::NodeMap<bool> &set,
                                    const ListGraph::Node x0,
                                    ListGraph::NodeMap<double> &dist) {
  ListDigraph G2;
  ListGraph::NodeMap<ListDigraph::Node> corr(G);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListDigraph::Node a = G2.addNode();
    corr[i] = a;
  }

  ListGraph::NodeMap<double> baseDist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(baseDist);
  dijkstra.init();
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (set[i]) {
      dijkstra.addSource(i);
    }
  }
  dijkstra.start();

  ListDigraph::ArcMap<double> coneLength(G2);
  ListDigraph::NodeMap<double> coneDist(G2);

  for (ListGraph::EdgeIt i(G); i != INVALID; ++i) {
    ListDigraph::Arc uv = G2.addArc(corr[G.u(i)], corr[G.v(i)]);
    ListDigraph::Arc vu = G2.addArc(corr[G.v(i)], corr[G.u(i)]);

    if (dijkstra.reached(G.u(i)) &&
        baseDist[G.u(i)] + length[i] == baseDist[G.v(i)])
      coneLength[uv] = 0;
    else
      coneLength[uv] = length[i];

    if (dijkstra.reached(G.v(i)) &&
        baseDist[G.v(i)] + length[i] == baseDist[G.u(i)])
      coneLength[vu] = 0;
    else
      coneLength[vu] = length[i];
  }

  Dijkstra<ListDigraph, ListDigraph::ArcMap<double>> dijkstra2(G2, coneLength);
  dijkstra2.distMap(coneDist);
  dijkstra2.init();
  dijkstra2.addSource(corr[x0]);
  dijkstra2.start();

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (!dijkstra2.reached(corr[i]))
      dist[i] = INF;
    else
      dist[i] = coneDist[corr[i]];
  }
}

double EmelElkinBase::BallCut(const ListGraph &G, const ListGraph::Node x0,
                              const ListGraph::EdgeMap<double> &length,
                              const double ro, const double delta,
                              ListGraph::NodeMap<bool> &Ball,
                              ListGraph::NodeMap<bool> &Shell) {
  if (countNodes(G) == 1) {
    Ball[x0] = true;
    Shell[x0] = false;
    return 0;
  }

  double r = ro * delta;
  double coeff =
      (log((double)(actM + 1)) / log((double)2)) / ((1 - 2 * delta) * ro);

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  ListGraph::NodeMap<long> id(G);
  NodeStruct *order = new NodeStruct[countNodes(G)];
  ListGraph::Node *inv = new ListGraph::Node[countNodes(G)];

  long j = 0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    id[i] = j;
    order[j].id = j;
    order[j].dist = dist[i];
    ListGraph::Node a = i;
    inv[j] = a;
    Ball[i] = false;
    Shell[i] = false;
    ++j;
  }

  qsort(order, countNodes(G), sizeof(NodeStruct), compareNS);

  j = 0;
  while (j < countNodes(G) && order[j].dist <= r + ERR)
    ++j;

  long until = j - 1;

  // Set initial values for vol and bound
  long vol = 0;
  double bound = 0;

  j = 0;
  while (j < until + 1) {
    ListGraph::Node node = inv[order[j].id];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.u(e) == node) {
        if (dist[G.v(e)] > r || id[G.u(e)] < id[G.v(e)])
          ++vol;
        if (dist[G.v(e)] > r)
          bound += 1 / length[e];
      } else {
        if (dist[G.u(e)] > r || id[G.v(e)] < id[G.u(e)])
          ++vol;
        if (dist[G.u(e)] > r)
          bound += 1 / length[e];
      }
    }
    ++j;
  }

  // Main loop
  while (bound > (vol + 1) * coeff) {
    long newOnes = 1;
    ++until;
    r = order[until].dist;
    while (until < countNodes(G) - 1 && order[until + 1].dist <= r + ERR) {
      ++until;
      ++newOnes;
    }
    if (until == countNodes(G) - 1)
      break;

    // refresh vol
    j = newOnes;
    while (j > 0) {
      ListGraph::Node &node = inv[order[until - j + 1].id];
      for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
        if (G.u(e) == node) {
          if (dist[G.v(e)] > r)
            ++vol;
        } else {
          if (dist[G.u(e)] > r)
            ++vol;
        }
      }
      --j;
    }

    // refresh bound
    bound = 0;
    j = 0;
    while (j != until + 1) {
      ListGraph::Node &node = inv[order[j].id];
      for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
        if (G.u(e) == node) {
          if (dist[G.v(e)] > r)
            bound += 1 / length[e];
        } else {
          if (dist[G.u(e)] > r)
            bound += 1 / length[e];
        }
      }
      ++j;
    }
  }

  j = 0;
  while (j <= until) {
    Ball[inv[order[j].id]] = true;
    ++j;
  }

  j = until + 1;
  while (j < countNodes(G)) {
    ListGraph::Node &node = inv[order[j].id];
    double jlk = dist[node];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      double lgt = length[e];
      if (G.u(e) == node) {
        if (dist[G.v(e)] <= r && dist[node] == dist[G.v(e)] + lgt)
          Shell[node] = true;
      } else {
        if (dist[G.u(e)] <= r && dist[node] == dist[G.u(e)] + lgt)
          Shell[node] = true;
      }
    }
    ++j;
  }

  delete[] order;
  delete[] inv;

  return r;
}

double EmelElkinBase::ConeCut(const ListGraph &G, const ListGraph::Node x0,
                              const ListGraph::EdgeMap<double> &length,
                              const double lambda1, const double lambda2,
                              const ListGraph::NodeMap<bool> &S,
                              ListGraph::NodeMap<bool> &Cone) {
  if (countNodes(G) == 1) {
    Cone[x0] = true;
    return 0;
  }

  double r = lambda1;

  ListGraph::NodeMap<double> dist(G);
  ListGraph::NodeMap<double> &distRef = dist;

  GetConeDistance(G, length, S, x0, distRef);

  ListGraph::NodeMap<long> id(G);
  NodeStruct *order = new NodeStruct[countNodes(G)];
  ListGraph::Node *inv = new ListGraph::Node[countNodes(G)];

  long j = 0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    id[i] = j;
    order[j].id = j;
    order[j].dist = dist[i];
    ListGraph::Node a = i;
    inv[j] = a;
    Cone[i] = false;
    ++j;
  }

  qsort(order, countNodes(G), sizeof(NodeStruct), compareNS);

  j = 0;
  while (j < countNodes(G) && order[j].dist <= r + ERR)
    ++j;

  long until = j - 1;

  // Set parameters (vol, volE) and initial boundary
  long vol = 0, volE = 0;
  double bound = 0;

  j = 0;
  while (j < until + 1) {
    ListGraph::Node &node = inv[order[j].id];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.u(e) == node) {
        if (dist[G.v(e)] > r || id[G.u(e)] < id[G.v(e)])
          ++vol;
        if (id[G.u(e)] < id[G.v(e)] && dist[G.v(e)] <= r)
          ++volE;
        if (dist[G.v(e)] > r)
          bound += 1 / length[e];
      } else {
        if (dist[G.v(e)] > r || id[G.v(e)] < id[G.u(e)])
          ++vol;
        if (id[G.v(e)] < id[G.u(e)] && dist[G.u(e)] <= r)
          ++volE;
        if (dist[G.u(e)] > r)
          bound += 1 / length[e];
      }
    }
    ++j;
  }

  double mu;
  if (volE == 0)
    mu = (vol + 1) * (log((double)actM + 1) / log((double)2));
  else
    mu = (vol) * (log((double)actM / (double)volE) / log((double)2));
  double coeff = mu / (lambda2 - lambda1);

  // main cycle
  while (bound > coeff) {
    long newOnes = 1;
    ++until;
    if (until == countNodes(G) ||
        order[until].dist >
            0.9 * INF) { // All nodes, or all nodes of the component
      --until;
      break;
    }
    r = order[until].dist;
    while (until < countNodes(G) - 1 && order[until + 1].dist <= r + ERR) {
      ++until;
      ++newOnes;
    }
    if (until == countNodes(G) - 1)
      break;

    // refresh bound
    bound = 0;
    j = 0;
    while (j < until + 1) {
      ListGraph::Node &node = inv[order[j].id];
      for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
        if (G.u(e) == node) {
          if (dist[G.v(e)] > r)
            bound += 1 / length[e];
        } else {
          if (dist[G.u(e)] > r)
            bound += 1 / length[e];
        }
      }
      ++j;
    }
  }

  // specify cone
  j = 0;
  while (j <= until) {
    Cone[inv[order[j].id]] = true;
    ++j;
  }

  delete[] order;

  return r;
}

Decomp EmelElkin::ConeDecomp(const ListGraph &G,
                             const ListGraph::EdgeMap<double> &length,
                             const ListGraph::NodeMap<bool> &S,
                             const double delta,
                             ListGraph::NodeMap<long> &whComp) {
  Decomp Dec;

  long k = 0;
  long inS = 0;

  ListGraph G2;
  ListGraph &G2ref = G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::NodeMap<ListGraph::Node> &refoldId = oldId;
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, oldId);

  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> chosenAsXi(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    S2[corr[i]] = S[i];
    if (S[i])
      ++inS;
  }
  for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
    whCompG2[i] = -1;
    chosenAsXi[i] = -1;
  }
  ListGraph::NodeMap<bool> &S2Ref = S2;

  ListGraph::NodeMap<bool> inCone(G2);
  ListGraph::NodeMap<bool> &inConeRef = inCone;

  list<ListGraph::Node> Xs;

  // Loop
  while (inS > 0) {
    ++k;

    int u = countNodes(G2);
    ListGraph::Node temp;

    for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
      if (S2[i] && whCompG2[i] == -1) {
        chosenAsXi[i] = k;
        temp = i;
        ListGraph::Node xi = i;
        Xs.push_back(oldId[xi]);
        break;
      }
    }

    ListGraph::Node &newx0 = temp;

    ConeCut(G2ref, newx0, lengthG2, 0, delta, S2Ref, inConeRef);

    list<ListGraph::Node> toDelete;

    for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
      if (inCone[i]) {
        whComp[oldId[i]] = k;
        if (S2[i])
          --inS;
        toDelete.push_back(i);
      }
    }
    while (!toDelete.empty()) {
      ListGraph::Node del = toDelete.front();
      toDelete.pop_front();
      G2.erase(del);
    }
  }

  Dec.x = new ListGraph::Node[k + 1];
  Dec.y = new ListGraph::Node[k + 1];
  Dec.k = k;
  long j = 1;
  while (j <= Dec.k) {
    Dec.x[j] = Xs.front();
    Xs.pop_front();
    ++j;
  }

  return Dec;
}

Decomp EmelElkin::StarDecomp(const ListGraph &G, ListGraph::Node x0,
                             const ListGraph::EdgeMap<double> &length,
                             const double delta, const double eps,
                             ListGraph::NodeMap<long> &whComp) {
  ListGraph::NodeMap<bool> inBall(G);
  ListGraph::NodeMap<bool> Shell(G);
  ListGraph::NodeMap<bool> S(G);
  double ro = GetRadius(G, length, x0);
  double r = BallCut(G, x0, length, ro, delta, inBall, Shell);

  ListGraph G2;
  ListGraph &G2ref = G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, oldId);

  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  int c = 0;
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      G2.erase(corr[i]);
    else
      S2[corr[i]] = Shell[i];
  }

  Decomp Dec1 = ConeDecomp(G2, lengthG2, S2, eps * ro / 2, whCompG2);

  Decomp Dec;
  Dec.k = Dec1.k;
  Dec.x = new ListGraph::Node[Dec.k + 1];
  Dec.y = new ListGraph::Node[Dec.k + 1];

  Dec.x[0] = x0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      whComp[i] = 0;
    else {
      long c = whCompG2[corr[i]];
      whComp[i] = whCompG2[corr[i]];
    }
  }
  long j = 1;
  while (j <= Dec1.k) {
    Dec.x[j] = oldId[Dec1.x[j]];
    Dec.y[j] = INVALID;
    ++j;
  }

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  j = 1;
  while (j <= Dec1.k) {
    ListGraph::Node &node = Dec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (inBall[G.u(e)] && approx(dist[G.v(e)], dist[G.u(e)] + length[e])) {
          Dec.y[j] = G.u(e);
          break;
        }
      } else {
        if (inBall[G.v(e)] && approx(dist[G.u(e)], dist[G.v(e)] + length[e])) {
          Dec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  j = 1;
  while (j <= Dec1.k) {
    if (Dec.y[j] != INVALID) {
      ++j;
      continue;
    }
    ListGraph::Node &node = Dec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (inBall[G.u(e)]) {
          Dec.y[j] = G.u(e);
          break;
        }
      } else {
        if (inBall[G.v(e)]) {
          Dec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  return Dec;
}

void EmelElkinBase::LowStrecthTree(
    ListGraph &G, ListGraph::Node x0, const ListGraph::EdgeMap<double> &length,
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

  actN = countNodes(G);
  actM = countEdges(G);

  ListGraph::NodeMap<long> whComp(G);

  double ro = GetRadius(G, length, x0);

  ListGraph G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> OldCorr(G2);
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, OldCorr);

  ListGraph::NodeMap<Nlist> from(G2);

  Tree.clear();
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = Tree.addNode();
    TNewId[i] = a;
    from[corr[i]].push_front(i);
  }
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    TOldId[TNewId[i]] = i;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    if (length[e] < beta * ro / N) {
      from[corr[G.u(e)]].splice(from[corr[G.u(e)]].end(), from[corr[G.v(e)]]);
      G2.contract(corr[G.u(e)], corr[G.v(e)], true);
      corr[G.v(e)] = corr[G.u(e)];
    }
  }

  ListGraph::NodeMap<long> whCompG2(G2);
  Decomp sd =
      StarDecomp(G2, corr[x0], lengthG2, (double)1 / (double)3, beta, whCompG2);

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    whComp[i] = whCompG2[corr[i]];
  }

  ListGraph *Gnew = new ListGraph[sd.k + 1];
  ListGraph::EdgeMap<double> **lengthGnew =
      new ListGraph::EdgeMap<double> *[sd.k + 1];
  ListGraph::NodeMap<ListGraph::Node> **oldId =
      new ListGraph::NodeMap<ListGraph::Node> *[sd.k + 1];
  long j = 0;
  while (j < sd.k + 1) {
    lengthGnew[j] = new ListGraph::EdgeMap<double>(Gnew[j]);
    oldId[j] = new ListGraph::NodeMap<ListGraph::Node>(Gnew[j]);
    ++j;
  }

  ListGraph::NodeMap<ListGraph::Node> newId(G);

  GraphSplit(G, whComp, sd.k, length, lengthGnew, Gnew, newId, oldId);

  j = 0;
  while (j < sd.k + 1) {
    ListGraph GjTree;
    ListGraph::EdgeMap<double> lengthNew(GjTree);
    ListGraph::NodeMap<ListGraph::Node> treeId(Gnew[j]);
    ListGraph::NodeMap<ListGraph::Node> origId(GjTree);

    ListGraph::EdgeMap<double> &lGnewj = *lengthGnew[j];

    ListGraph::Node a = newId[sd.x[j]];

    LowStrecthTree(Gnew[j], newId[OldCorr[sd.x[j]]], lGnewj, GjTree, lengthNew,
                   treeId, origId);

    for (ListGraph::EdgeIt e(GjTree); e != INVALID; ++e) {
      ListGraph::Node a = (*oldId[j])[origId[GjTree.u(e)]];
      ListGraph::Node b = (*oldId[j])[origId[GjTree.v(e)]];
      ListGraph::Edge e2 = Tree.addEdge(TNewId[a], TNewId[b]);
      lengthTNew[e2] = lengthNew[e];
    }

    ++j;
  }

  // Find xi-yi edges
  j = 1;
  while (j < sd.k + 1) {
    ListGraph::Node x, y;
    double min = INF;
    int count = 0;
    for (std::list<ListGraph::Node>::const_iterator it = from[sd.x[j]].begin(),
                                                    end = from[sd.x[j]].end();
         it != end; ++it) {
      ListGraph::Node st = *it;
      for (ListGraph::IncEdgeIt e(G, st); e != INVALID; ++e) {
        if (G.u(e) == st && Contains(from[sd.y[j]], G.v(e))) {
          if (length[e] < min) {
            min = length[e];
            x = st;
            y = G.v(e);
          }
          ++count;
        }
        if (G.v(e) == st && Contains(from[sd.y[j]], G.u(e))) {
          if (length[e] < min) {
            min = length[e];
            x = st;
            y = G.u(e);
          }
          ++count;
        }
      }
    }
    int h = countNodes(G);
    ListGraph::Edge e = Tree.addEdge(TNewId[x], TNewId[y]);
    lengthTNew[e] = min;

    ++j;
  }

  j = 0;
  while (j < sd.k + 1) {
    delete lengthGnew[j];
    delete oldId[j];
    ++j;
  }

  delete[] Gnew;
}

Decomp EmelElkinPlus::StarDecomp(const ListGraph &G, ListGraph::Node x0,
                                 const ListGraph::EdgeMap<double> &length,
                                 const double delta, const double eps,
                                 ListGraph::NodeMap<long> &whComp) {
  ListGraph::NodeMap<bool> inBall(G);
  ListGraph::NodeMap<bool> Shell(G);
  ListGraph::NodeMap<bool> S(G);
  double ro = GetRadius(G, length, x0);
  double r = BallCut(G, x0, length, ro, delta, inBall, Shell);

  ListGraph G2;
  ListGraph &G2ref = G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, oldId);

  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  int c = 0;
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      G2.erase(corr[i]);
    else
      S2[corr[i]] = Shell[i];
  }

  ListGraph::Node useless = INVALID;

  Decomp Dec1 =
      ImpConeDecomp(G2, lengthG2, S2, eps * ro / 2,
                    (long)round(2 * log(log((double)M))), whCompG2, useless);

  Decomp Dec;
  Dec.k = Dec1.k;
  Dec.x = new ListGraph::Node[Dec.k + 1];
  Dec.y = new ListGraph::Node[Dec.k + 1];

  Dec.x[0] = x0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      whComp[i] = 0;
    else {
      long c = whCompG2[corr[i]];
      whComp[i] = whCompG2[corr[i]];
    }
  }
  long j = 1;
  while (j <= Dec1.k) {
    Dec.x[j] = oldId[Dec1.x[j]];
    Dec.y[j] = INVALID;
    ++j;
  }

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  j = 1;
  while (j <= Dec1.k) {
    ListGraph::Node &node = Dec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (inBall[G.u(e)] && approx(dist[G.v(e)], dist[G.u(e)] + length[e])) {
          Dec.y[j] = G.u(e);
          break;
        }
      } else {
        if (inBall[G.v(e)] && approx(dist[G.u(e)], dist[G.v(e)] + length[e])) {
          Dec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  j = 1;
  while (j <= Dec1.k) {
    if (Dec.y[j] != INVALID) {
      ++j;
      continue;
    }
    ListGraph::Node &node = Dec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (inBall[G.u(e)]) {
          Dec.y[j] = G.u(e);
          break;
        }
      } else {
        if (inBall[G.v(e)]) {
          Dec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  return Dec;
}

Decomp EmelElkinPlus::ImpConeDecomp(const ListGraph &G,
                                    const ListGraph::EdgeMap<double> &length,
                                    const ListGraph::NodeMap<bool> &S,
                                    const double delta, const long t,
                                    ListGraph::NodeMap<long> &whComp,
                                    ListGraph::Node &x1, bool givenX1) {

  Decomp Dec;

  long k = 0;
  long inS = 0;

  ListGraph G2;
  ListGraph &G2ref = G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::NodeMap<ListGraph::Node> &refoldId = oldId;
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, oldId);

  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> chosenAsXi(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    S2[corr[i]] = S[i];
    if (S[i])
      ++inS;
  }
  for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
    whCompG2[i] = -1;
    chosenAsXi[i] = -1;
  }
  ListGraph::NodeMap<bool> &S2Ref = S2;

  ListGraph::NodeMap<bool> inCone(G2);
  ListGraph::NodeMap<bool> &inConeRef = inCone;

  list<ListGraph::Node> Xs;

  // Loop
  while (inS > 0) {
    ++k;

    ListGraph::Node temp;

    if (k == 1 && givenX1) {
      chosenAsXi[corr[x1]] = k;
      temp = corr[x1];
      Xs.push_back(x1);
    } else {
      for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
        if (S2[i] && whCompG2[i] == -1) {
          chosenAsXi[i] = k;
          temp = i;
          ListGraph::Node xi = i;
          Xs.push_back(oldId[xi]);
          break;
        }
      }
    }

    ListGraph::Node &newx0 = temp;

    long p = t - 1;
    if (p < 1)
      p = 1;

    while (p > 0) {
      int tp = t - p - 1;
      if (tp < 0)
        tp = 0;

      ConeCut(G2ref, newx0, lengthG2, ((double)tp / (double)t) * delta,
              ((double)(tp + 1) / (double)t) * delta, S2Ref, inConeRef);
      long volE = 0;
      for (ListGraph::EdgeIt e(G2); e != INVALID; ++e) {
        if (inCone[G2.u(e)] && inCone[G2.v(e)])
          ++volE;
      }
      if ((double)volE <=
          (double)actM / (pow(2, pow(log((double)M), (double)p / (double)t))))
        break;
      else
        --p;
    }

    list<ListGraph::Node> toDelete;

    for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
      if (inCone[i]) {
        whComp[oldId[i]] = k;
        if (S2[i])
          --inS;
        toDelete.push_back(i);
      }
    }
    while (!toDelete.empty()) {
      ListGraph::Node del = toDelete.front();
      toDelete.pop_front();
      G2.erase(del);
    }
  }

  Dec.x = new ListGraph::Node[k + 1];
  Dec.y = new ListGraph::Node[k + 1];
  Dec.k = k;
  long j = 1;
  while (j <= Dec.k) {
    Dec.x[j] = Xs.front();
    Xs.pop_front();
    ++j;
  }

  return Dec;
}
