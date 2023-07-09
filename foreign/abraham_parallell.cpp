#include <iostream>
#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>
#include <list>


// #include <parlay/io.h>
#include "parlay/io.h"
#include "parlay/primitives.h"
#include "parlay/random.h"

#include "AbrahamBartal_parallel.h"

using namespace lemon;
using namespace std;

void AbrahamBartal_par::FirstHierarStarPart(
    ListGraph &G, ListGraph::Node x0, const ListGraph::EdgeMap<double> &length,
    list<ListGraph::Node> Q, ListGraph &Tree,
    ListGraph::EdgeMap<double> &lengthTNew,
    ListGraph::NodeMap<ListGraph::Node> &TNewId,
    ListGraph::NodeMap<ListGraph::Node> &TOldId) {
  if (countNodes(G) <= 16 * c) {
    ListGraph::NodeMap<ListGraph::Node> &refnewId = TNewId;
    ListGraph::NodeMap<ListGraph::Node> &refoldId = TOldId;

    for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
      ListGraph::Node a = Tree.addNode();
      TNewId[i] = a;
      TOldId[a] = i;
    }

    ListGraph::NodeMap<double> dist(G);

    Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
    dijkstra.distMap(dist);
    dijkstra.init();
    dijkstra.addSource(x0);
    dijkstra.start();

    for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
      ListGraph::Node a = i;
      if (i == x0)
        continue;
      ListGraph::Edge e = dijkstra.predArc(i);
      ListGraph::Edge NewE = Tree.addEdge(TNewId[G.u(e)], TNewId[G.v(e)]);
      lengthTNew[NewE] = length[e];
    }

    return;
  }

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

  list<ListGraph::Node> Q2;
  ListGraph::NodeMap<bool> PutIn(G2);
  for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
    PutIn[i] = false;
  }

  for (std::list<ListGraph::Node>::iterator it = Q.begin(); it != Q.end();
       it++) {
    ListGraph::Node a = *it;
    ListGraph::Node b = corr[a];
    if (!PutIn[b]) {
      Q2.push_back(corr[a]);
      PutIn[corr[a]] = false;
    }
  }

  ListGraph::NodeMap<long> whCompG2(G2);

  QDecomp sd = StarPart(G2, corr[x0], lengthG2, Q2, whCompG2);

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

  ListGraph *GjTree = new ListGraph[sd.k + 1];
  ListGraph::EdgeMap<double> **lengthNew =
      new ListGraph::EdgeMap<double> *[sd.k + 1];
  ListGraph::NodeMap<ListGraph::Node> **origId =
      new ListGraph::NodeMap<ListGraph::Node> *[sd.k + 1];
  j = 0;
  while (j < sd.k + 1) {
    lengthNew[j] = new ListGraph::EdgeMap<double>(GjTree[j]);
    origId[j] = new ListGraph::NodeMap<ListGraph::Node>(GjTree[j]);
    ++j;
  }

  parlay::parallel_for(0, sd.k + 1, [&] (size_t j) {
    ListGraph::NodeMap<ListGraph::Node> treeId(Gnew[j]);

    ListGraph::EdgeMap<double> &lGnewj = *lengthGnew[j];

    ListGraph::Node a = newId[sd.x[j]];

    list<ListGraph::Node> Qj;

    for (std::list<ListGraph::Node>::iterator it = sd.Q[j].begin();
         it != sd.Q[j].end(); it++) {
      ListGraph::Node a = *it;

      for (std::list<ListGraph::Node>::iterator it2 = from[a].begin(),
                                                end = from[a].end();
           it2 != end; ++it2) {
        ListGraph::Node b = *it2;
        Qj.push_back(newId[b]);
      }
    }

    FirstHierarStarPart(Gnew[j], newId[OldCorr[sd.x[j]]], lGnewj, Qj, GjTree[j],
                        *lengthNew[j], treeId, *origId[j]);
    }, 1);

  for (int j = 0; j < sd.k + 1; ++j) {
    for (ListGraph::EdgeIt e(GjTree[j]); e != INVALID; ++e) {
      ListGraph::Node a = (*oldId[j])[(*origId[j])[GjTree[j].u(e)]];
      ListGraph::Node b = (*oldId[j])[(*origId[j])[GjTree[j].v(e)]];
      ListGraph::Edge e2 =
          Tree.addEdge(TNewId[a], TNewId[b]); // This is what fails

      lengthTNew[e2] = (*lengthNew[j])[e];
    }
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
    if (lengthGnew[j]) {
      delete lengthGnew[j];
      lengthGnew[j] = NULL;
    }
    if (oldId[j]) {
      delete oldId[j];
      oldId[j] = NULL;
    }
    ++j;
  }

  if (Gnew) {
    delete[] Gnew;
    Gnew = NULL;
  }
}

QDecomp AbrahamBartal_par::StarPart(const ListGraph &G, ListGraph::Node x0,
                                    const ListGraph::EdgeMap<double> &length,
                                    const list<ListGraph::Node> Q,
                                    ListGraph::NodeMap<long> &whComp) {

  double r0 = 1 / (8 * (double)c);
  double eps = 1 / (32 * log(log((double)countNodes(G))));

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  ListGraph G2;
  ListGraph::NodeMap<ListGraph::Node> newId(G);
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, newId, oldId);

  ListGraph::NodeMap<double> distV(G2);
  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> whCompG2(G2);

  double rad = 0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = i;
    if (dist[a] > rad)
      rad = dist[a];
  }

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = i;
    if (dist[a] <= r0 * rad + ERR) {
      G2.erase(newId[a]);
      whComp[a] = 0;
    } else {
      whComp[a] = -1;
      S2[newId[a]] = false;
      distV[newId[a]] = dist[a];
    }
  }

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = i;
    if (whComp[a] == -1) {
      for (ListGraph::IncEdgeIt e(G, a); e != INVALID; ++e) {
        if ((G.u(e) == a && whComp[G.v(e)] == 0) ||
            (G.v(e) == a && whComp[G.u(e)] == 0)) {
          S2[newId[a]] = true;
          break;
        }
      }
    }
  }

  int cc1 = whComp[Q.front()];

  EmelElkinCorrPlus ee(G2, lengthG2);

  ListGraph::Node bridge1 = INVALID;

  Decomp Dec1;
  ListGraph::Node gg = Q.front();
  if (whComp[Q.front()] == 0) {
    Dec1 = ee.ImpConeDecomp(G2, lengthG2, S2, rad / log(log((double)N)),
                            (long)round(log(log((double)N))), distV, whCompG2,
                            bridge1);
  } else {
    ListGraph::Node act = Q.front();
    while (whComp[dijkstra.predNode(act)] == -1)
      act = dijkstra.predNode(act);
    bridge1 = newId[act];
    Dec1 = ee.ImpConeDecomp(G2, lengthG2, S2, rad / log(log((double)N)),
                            (long)round(log(log((double)N))), distV, whCompG2,
                            bridge1, true);
  }

  // Create queues

  QDecomp qdec;
  qdec.k = Dec1.k;
  qdec.Q = new list<ListGraph::Node>[Dec1.k + 1];
  qdec.x = new ListGraph::Node[Dec1.k + 1];
  qdec.y = new ListGraph::Node[Dec1.k + 1];

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = i;
    if (whComp[a] == -1) {
      whComp[a] = whCompG2[newId[a]];
    }
  }

  qdec.x[0] = x0;

  long j = 1;
  while (j <= Dec1.k) {
    qdec.x[j] = oldId[Dec1.x[j]];
    qdec.y[j] = INVALID;
    ++j;
  }

  j = 1;
  while (j <= qdec.k) {
    ListGraph::Node &node = qdec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (whComp[G.u(e)] == 0 &&
            approx(dist[G.v(e)], dist[G.u(e)] + length[e])) {
          qdec.y[j] = G.u(e);
          break;
        }
      } else {
        if (whComp[G.v(e)] == 0 &&
            approx(dist[G.u(e)], dist[G.v(e)] + length[e])) {
          qdec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  j = 1;
  while (j <= qdec.k) {
    if (qdec.y[j] != INVALID) {
      ++j;
      continue;
    }
    ListGraph::Node &node = qdec.x[j];
    for (ListGraph::IncEdgeIt e(G, node); e != INVALID; ++e) {
      if (G.v(e) == node) {
        if (whComp[G.u(e)] == 0) {
          qdec.y[j] = G.u(e);
          break;
        }
      } else {
        if (whComp[G.v(e)] == 0) {
          qdec.y[j] = G.v(e);
          break;
        }
      }
    }
    ++j;
  }

  list<ListGraph::Node> Qball, Qfat, Qreg;

  long *cnt = new long[qdec.k + 1];

  j = 0;
  while (j < qdec.k + 1) {
    cnt[j] = 0;
    ++j;
  }

  j = 1;
  for (std::list<ListGraph::Node>::const_iterator it = Q.begin(), end = Q.end();
       it != end; ++it) {
    ListGraph::Node a = *it;
    if (whComp[a] == 0)
      Qball.push_back(a);
    else {
      ListGraph::Node frt = qdec.y[whComp[a]];
      if (a != qdec.x[whComp[a]])
        qdec.Q[whComp[a]].push_back(a);
      if (!Contains(Qreg, qdec.y[whComp[a]]))
        Qreg.push_back(qdec.y[whComp[a]]);

      ++cnt[whComp[a]];
      if ((double)cnt[whComp[a]] > sqrt((double)j) &&
          !Contains(Qfat, qdec.y[whComp[a]]))
        Qfat.push_back(qdec.y[whComp[a]]);
      ++j;
    }
  }

  if (whComp[Q.front()] == 0)
    qdec.Q[0].push_front(Q.front());
  else
    qdec.Q[0].push_front(qdec.y[1]);

  j = 1;
  while (!Qball.empty() || !Qfat.empty() || !Qreg.empty()) {
    if (!Qball.empty()) {
      ListGraph::Node a = Qball.front();
      Qball.pop_front();
      if (!Contains(qdec.Q[0], a))
        qdec.Q[0].push_back(a);
    }
    if (!Qreg.empty()) {
      ListGraph::Node a = Qreg.front();
      Qreg.pop_front();
      if (!Contains(qdec.Q[0], a))
        qdec.Q[0].push_back(a);
    }
    if (!Qfat.empty()) {
      ListGraph::Node a = Qfat.front();
      Qfat.pop_front();
      if (!Contains(qdec.Q[0], a))
        qdec.Q[0].push_back(a);
    }
  }

  if (cnt) {
    delete[] cnt;
    cnt = NULL;
  }
  return qdec;
}
