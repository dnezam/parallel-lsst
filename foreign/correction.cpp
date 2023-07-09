#include <iostream>
#include <lemon/connectivity.h>
#include <lemon/core.h>
#include <lemon/dijkstra.h>
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

Decomp EmelElkinCorr::StarDecomp(const ListGraph &G, ListGraph::Node x0,
                                 const ListGraph::EdgeMap<double> &length,
                                 const double delta, const double eps,
                                 ListGraph::NodeMap<long> &whComp) {
  ListGraph::NodeMap<bool> inBall(G);
  ListGraph::NodeMap<bool> Shell(G);
  ListGraph::NodeMap<bool> S(G);

  double ro = GetRadius(G, length, x0);
  double r = BallCut(G, x0, length, ro, delta, inBall, Shell);

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  ListGraph G2;
  ListGraph::NodeMap<ListGraph::Node> corr(G);
  ListGraph::NodeMap<ListGraph::Node> &refCorr = corr;
  ListGraph::NodeMap<ListGraph::Node> oldId(G2);
  ListGraph::EdgeMap<double> lengthG2(G2);
  UsableGraphCopy(G, G2, length, lengthG2, refCorr, oldId);

  ListGraph::NodeMap<double> distV(G2);
  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  int c = 0;
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      G2.erase(corr[i]);
    else {
      S2[corr[i]] = Shell[i];
      distV[corr[i]] = dist[i];
    }
  }

  Decomp Dec1 = ConeDecomp(G2, lengthG2, S2, eps * ro, distV, whCompG2);

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

Decomp EmelElkinCorr::ConeDecomp(const ListGraph &G,
                                 const ListGraph::EdgeMap<double> &length,
                                 const ListGraph::NodeMap<bool> &S,
                                 const double delta,
                                 const ListGraph::NodeMap<double> &distV,
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

  ListGraph::NodeMap<double> distV2(G2);
  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> chosenAsXi(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    S2[corr[i]] = S[i];
    distV2[corr[i]] = distV[i];
    if (S[i])
      ++inS;
    whComp[i] = -1;
  }
  for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
    whCompG2[i] = -1;
    chosenAsXi[i] = -1;
  }

  ListGraph::NodeMap<bool> &S2Ref = S2;

  ListGraph::NodeMap<bool> inCone(G2);
  ListGraph::NodeMap<bool> &inConeRef = inCone;

  ListGraph::NodeMap<int> com(G2);

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

    ConeCut(G2ref, newx0, lengthG2, 0, delta, S2Ref, distV2, inConeRef);

    list<ListGraph::Node> toDelete;

    for (ListGraph::NodeIt i(G2); i != INVALID; ++i) {
      if (inCone[i]) {
        whComp[oldId[i]] = k;
        whCompG2[i] = k;
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

double EmelElkinCorr::ConeCut(const ListGraph &G, const ListGraph::Node x0,
                              const ListGraph::EdgeMap<double> &length,
                              const double lambda1, const double lambda2,
                              const ListGraph::NodeMap<bool> &S,
                              const ListGraph::NodeMap<double> &distV,
                              ListGraph::NodeMap<bool> &Cone) {
  if (countNodes(G) == 1) {
    Cone[x0] = true;
    return 0;
  }

  double r = lambda1;

  ListGraph::NodeMap<double> realdist(G);
  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(realdist);
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
    if (dijkstra.reached(i))
      dist[i] = distV[x0] + realdist[i] - distV[i];
    else
      dist[i] = INF;
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

//
//
//

Decomp EmelElkinCorrPlus::StarDecomp(const ListGraph &G, ListGraph::Node x0,
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

  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  ListGraph::NodeMap<double> distV(G2);
  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  int c = 0;
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (inBall[i])
      G2.erase(corr[i]);
    else {
      S2[corr[i]] = Shell[i];
      distV[corr[i]] = dist[i];
    }
  }

  ListGraph::Node useless = INVALID;

  Decomp Dec1 = ImpConeDecomp(G2, lengthG2, S2, eps * ro,
                              (long)round(2 * log(log((double)M))), distV,
                              whCompG2, useless);

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

Decomp EmelElkinCorrPlus::ImpConeDecomp(
    const ListGraph &G, const ListGraph::EdgeMap<double> &length,
    const ListGraph::NodeMap<bool> &S, const double delta, const long t,
    const ListGraph::NodeMap<double> &distV, ListGraph::NodeMap<long> &whComp,
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

  ListGraph::NodeMap<double> distV2(G2);
  ListGraph::NodeMap<bool> S2(G2);
  ListGraph::NodeMap<long> chosenAsXi(G2);
  ListGraph::NodeMap<long> whCompG2(G2);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    S2[corr[i]] = S[i];
    if (S[i])
      ++inS;
    distV2[corr[i]] = distV[i];
    whComp[i] = -1;
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
              ((double)(tp + 1) / (double)t) * delta, S2Ref, distV2, inConeRef);
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
        whCompG2[i] = k;
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

double EmelElkinCorrPlus::ConeCut(const ListGraph &G, const ListGraph::Node x0,
                                  const ListGraph::EdgeMap<double> &length,
                                  const double lambda1, const double lambda2,
                                  const ListGraph::NodeMap<bool> &S,
                                  const ListGraph::NodeMap<double> &distV,
                                  ListGraph::NodeMap<bool> &Cone) {
  if (countNodes(G) == 1) {
    Cone[x0] = true;
    return 0;
  }

  double r = lambda1;

  ListGraph::NodeMap<double> realdist(G);
  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(realdist);
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
    if (dijkstra.reached(i))
      dist[i] = distV[x0] + realdist[i] - distV[i];
    else
      dist[i] = INF;
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
