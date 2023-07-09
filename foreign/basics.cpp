#include <lemon/connectivity.h>
#include <lemon/dijkstra.h>
#include <lemon/random.h>
#include <time.h>

#include "basics.h"

int compareNS(const void *a, const void *b) {
  double k = ((NodeStruct *)a)->dist - ((NodeStruct *)b)->dist;
  if (k < 0)
    return -1;
  if (k > 0)
    return 1;
  return 0;
};

int compareNSI(const void *a, const void *b) {
  int k = ((NodeStructI *)a)->dist - ((NodeStructI *)b)->dist;
  if (k < 0)
    return -1;
  if (k > 0)
    return 1;
  return 0;
};

int compareNP(const void *a, const void *b) {
  double k = ((NodePair *)a)->dist - ((NodePair *)b)->dist;
  if (k < 0)
    return -1;
  if (k > 0)
    return 1;
  return 0;
};

bool Contains(list<ListGraph::Node> l, ListGraph::Node a) {
  for (list<ListGraph::Node>::const_iterator it = l.begin(), end = l.end();
       it != end; ++it) {
    if (*it == a)
      return true;
  }
  return false;
}

bool approx(double a, double b) {
  if (a - b >= (-1) * ERR && a - b <= ERR)
    return true;
  else
    return false;
}

double round(double x) { return x >= 0.0 ? floor(x + 0.5) : ceil(x - 0.5); }

double GetRadius(const ListGraph &G, const ListGraph::EdgeMap<double> &length,
                 const ListGraph::Node x0) {
  ListGraph::NodeMap<double> dist(G);

  Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
  dijkstra.distMap(dist);
  dijkstra.init();
  dijkstra.addSource(x0);
  dijkstra.start();

  double max = 0;

  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    if (dist[i] > max)
      max = dist[i];
  }
  return max;
}

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId) {
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = G2.addNode();
    newId[i] = a;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    ListGraph::Edge e2 = G2.addEdge(newId[G.u(e)], newId[G.v(e)]);
    lengthG2[e2] = length[e];
  }
}

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId,
                     ListGraph::NodeMap<ListGraph::Node> &oldId) {
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = G2.addNode();
    newId[i] = a;
    oldId[a] = i;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    ListGraph::Edge e2 = G2.addEdge(newId[G.u(e)], newId[G.v(e)]);
    lengthG2[e2] = length[e];
  }
}

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId,
                     ListGraph::NodeMap<ListGraph::Node> &oldId,
                     ListGraph::EdgeMap<ListGraph::Edge> &oldEdgeId) {
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = G2.addNode();
    newId[i] = a;
    oldId[a] = i;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    ListGraph::Edge e2 = G2.addEdge(newId[G.u(e)], newId[G.v(e)]);
    lengthG2[e2] = length[e];
    oldEdgeId[e2] = e;
  }
}

void GraphSplit(const ListGraph &G, const ListGraph::NodeMap<long> &whComp,
                const long k, const ListGraph::EdgeMap<double> &length,
                ListGraph::EdgeMap<double> **lengthGnew, ListGraph *Gnew,
                ListGraph::NodeMap<ListGraph::Node> &newId,
                ListGraph::NodeMap<ListGraph::Node> **oldId) {
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = Gnew[whComp[i]].addNode();
    newId[i] = a;
    (*(oldId[whComp[i]]))[a] = i;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    if (whComp[G.u(e)] == whComp[G.v(e)]) {
      ListGraph::Edge e2 =
          Gnew[whComp[G.u(e)]].addEdge(newId[G.u(e)], newId[G.v(e)]);
      (*(lengthGnew[whComp[G.u(e)]]))[e2] = length[e];
    }
  }
}

void GraphSplit(const ListGraph &G, const ListGraph::NodeMap<long> &whComp,
                const long k, const ListGraph::EdgeMap<double> &length,
                ListGraph::EdgeMap<double> **lengthGnew, ListGraph *Gnew,
                ListGraph::NodeMap<ListGraph::Node> &newId,
                ListGraph::NodeMap<ListGraph::Node> **oldId,
                ListGraph::EdgeMap<ListGraph::Edge> **oldEdgeId) {
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    ListGraph::Node a = Gnew[whComp[i]].addNode();
    newId[i] = a;
    (*(oldId[whComp[i]]))[a] = i;
  }

  for (ListGraph::EdgeIt e(G); e != INVALID; ++e) {
    if (whComp[G.u(e)] == whComp[G.v(e)]) {
      ListGraph::Edge e2 =
          Gnew[whComp[G.u(e)]].addEdge(newId[G.u(e)], newId[G.v(e)]);
      (*(lengthGnew[whComp[G.u(e)]]))[e2] = length[e];
      (*(oldEdgeId[whComp[G.u(e)]]))[e2] = e;
    }
  }
}

void Create01Graph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const double dist) {
  G.clear();
  int j = 0;
  while (j < n) {
    G.addNode();
    ++j;
  }

  int parts = 3;
  if (n > 200)
    parts = 5;

  ListGraph::NodeMap<double> x(G), y(G);
  ListGraph::NodeMap<int> px(G), py(G);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    x[i] = (double)rand() / (double)RAND_MAX;
    y[i] = (double)rand() / (double)RAND_MAX;
    px[i] = (int)floor((double)parts * x[i]);
    py[i] = (int)floor((double)parts * y[i]);
  }

  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      double d =
          sqrt((x[u] - x[v]) * (x[u] - x[v]) + (y[u] - y[v]) * (y[u] - y[v]));
      if (u < v && d <= dist) {
        ListGraph::Edge e = G.addEdge(u, v);
        length[e] = d;
        if (approx(length[e], 0))
          length[e] = 0.0005;
      }
    }
  }

  ListGraph::NodeMap<int> com(G);
  int cNum = connectedComponents(G, com);
  if (cNum == 1)
    return;

  long prs = 0;
  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      if (u < v && com[u] != com[v] && abs(px[u] - px[v]) <= 1 &&
          abs(py[u] - py[v]) <= 1) {
        ++prs;
      }
    }
  }

  NodePair *NP = new NodePair[prs];
  long c = 0;
  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      if (u < v && com[u] != com[v] && abs(px[u] - px[v]) <= 1 &&
          abs(py[u] - py[v]) <= 1) {
        NP[c].u = a;
        NP[c].v = b;
        NP[c].dist =
            sqrt((x[u] - x[v]) * (x[u] - x[v]) + (y[u] - y[v]) * (y[u] - y[v]));
        ++c;
      }
    }
  }

  qsort(NP, prs, sizeof(NodePair), compareNP);

  j = 0;

  while (j < prs) {
    if (com[NP[j].u] != com[NP[j].v]) {
      ListGraph::Edge e = G.addEdge(NP[j].u, NP[j].v);
      length[e] = NP[j].dist;
      if (approx(length[e], 0))
        length[e] = 0.0005;

      cNum = connectedComponents(G, com);
      if (cNum == 1)
        break;
    }
    ++j;
  }

  delete[] NP;

  cNum = connectedComponents(G, com);
  if (cNum == 1)
    return;

  // In the completely unprobable case of not being connected,
  // because components were completely separated:

  while (cNum > 1) {
    double min = 2;
    cNum = connectedComponents(G, com);
    ListGraph::Node n1, n2;
    for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
      ListGraph::Node a = u;
      for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
        ListGraph::Node b = v;
        if (u < v && com[u] != com[v] &&
            sqrt((x[u] - x[v]) * (x[u] - x[v]) +
                 (y[u] - y[v]) * (y[u] - y[v])) < min) {
          n1 = a;
          n2 = b;
          min = sqrt((x[u] - x[v]) * (x[u] - x[v]) +
                     (y[u] - y[v]) * (y[u] - y[v]));
        }
      }
    }
    ListGraph::Edge e = G.addEdge(n1, n2);
    length[e] = min;
    if (approx(length[e], 0))
      length[e] = 0.0005;

    --cNum;
  }
}

void Create01Graph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const int k) {
  G.clear();
  int j = 0;
  while (j < n) {
    G.addNode();
    ++j;
  }

  int parts = 3;
  if (n > 200)
    parts = 5;

  ListGraph::NodeMap<double> x(G), y(G);
  ListGraph::NodeMap<int> px(G), py(G);
  for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
    x[i] = (double)rand() / (double)RAND_MAX;
    y[i] = (double)rand() / (double)RAND_MAX;
    px[i] = (int)floor((double)parts * x[i]);
    py[i] = (int)floor((double)parts * y[i]);
  }

  ListGraph::NodeMap<int> com(G);
  int cNum = connectedComponents(G, com);

  long prs = 0;
  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      if (u < v && abs(px[u] - px[v]) <= 1 && abs(py[u] - py[v]) <= 1) {
        ++prs;
      }
    }
  }

  NodePair *NP = new NodePair[prs];
  long c = 0;
  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      if (u < v && abs(px[u] - px[v]) <= 1 && abs(py[u] - py[v]) <= 1) {
        NP[c].u = a;
        NP[c].v = b;
        NP[c].dist =
            sqrt((x[u] - x[v]) * (x[u] - x[v]) + (y[u] - y[v]) * (y[u] - y[v]));
        ++c;
      }
    }
  }

  qsort(NP, prs, sizeof(NodePair), compareNP);

  j = 0;
  long Edg = 0;

  while (j < prs) {
    ListGraph::Edge e = G.addEdge(NP[j].u, NP[j].v);
    length[e] = NP[j].dist;
    if (approx(length[e], 0))
      length[e] = 0.0005;
    ++Edg;
    if (2 * Edg >= n * (long)k)
      break;

    ++j;
  }

  cNum = connectedComponents(G, com);
  if (cNum == 1) {
    delete[] NP;
    return;
  }

  while (j < prs) {
    if (com[NP[j].u] != com[NP[j].v]) {
      ListGraph::Edge e = G.addEdge(NP[j].u, NP[j].v);
      length[e] = NP[j].dist;
      if (approx(length[e], 0))
        length[e] = 0.0005;

      cNum = connectedComponents(G, com);
      if (cNum == 1)
        break;
    }
    ++j;
  }

  delete[] NP;

  cNum = connectedComponents(G, com);
  if (cNum == 1)
    return;

  // In the completely unprobable case of not being connected,
  // because components were completely separated:

  while (cNum > 1) {
    double min = 2;
    cNum = connectedComponents(G, com);
    ListGraph::Node n1, n2;
    for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
      ListGraph::Node a = u;
      for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
        ListGraph::Node b = v;
        if (u < v && com[u] != com[v] &&
            sqrt((x[u] - x[v]) * (x[u] - x[v]) +
                 (y[u] - y[v]) * (y[u] - y[v])) < min) {
          n1 = a;
          n2 = b;
          min = sqrt((x[u] - x[v]) * (x[u] - x[v]) +
                     (y[u] - y[v]) * (y[u] - y[v]));
        }
      }
    }
    ListGraph::Edge e = G.addEdge(n1, n2);
    length[e] = min;
    if (approx(length[e], 0))
      length[e] = 0.0005;

    --cNum;
  }
}

void CreateStrangeGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                        const int ring, const int side) {
  ListGraph::Node a = G.addNode();
  ListGraph::Node b = G.addNode();
  ListGraph::Edge e = G.addEdge(a, b);
  length[e] = 3;
  int i = 0;
  while (i < side) {
    ListGraph::Node t1 = G.addNode();
    ListGraph::Edge e1 = G.addEdge(a, t1);
    length[e1] = 1;
    ListGraph::Node t2 = G.addNode();
    ListGraph::Edge e2 = G.addEdge(a, t2);
    length[e2] = 1;
    ++i;
  }
  ListGraph::Node prev = a;
  ListGraph::Node next;
  i = 0;
  while (i < ring) {
    next = G.addNode();
    e = G.addEdge(prev, next);
    length[e] = 2;
    prev = next;
    ++i;
  }
  e = G.addEdge(prev, b);
  length[e] = 2;
}

void CreateERGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const double p) {
  G.clear();

  int j = 0;
  while (j < n) {
    G.addNode();
    ++j;
  }

  // Generate Hamiltonian cycle

  ListGraph::Node first = INVALID;
  ListGraph::Node last = INVALID;

  ListGraph::NodeMap<ListGraph::Node> prev(G);
  ListGraph::NodeMap<ListGraph::Node> next(G);

  for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
    if (first == INVALID) {
      first = v;
      last = v;
    } else {
      ListGraph::Edge e = G.addEdge(last, v);
      length[e] = rnd();
      if (approx(length[e], 0))
        length[e] = 0.0005;

      last = v;
      next[last] = v;
      prev[v] = last;
    }
  }

  ListGraph::Edge e = G.addEdge(first, last);
  length[e] = rnd();
  if (approx(length[e], 0))
    length[e] = 0.0005;

  prev[first] = first;
  next[last] = first;

  // Generate the rest of the edges

  // p has to be modified, since we already have some edges

  double p2 = (((double)n - 1) * p - 2) / ((double)n - 3);
  if (p2 < 0)
    p2 = 0;

  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      if (u < v && rnd() < p2) {
        ListGraph::Edge e = G.addEdge(u, v);
        length[e] = rnd();
        if (approx(length[e], 0))
          length[e] = 0.0005;
      }
    }
  }
}

void CreateERGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const int k) {
  CreateERGraph(G, length, n, (double)k / (double)(n - 1));
}