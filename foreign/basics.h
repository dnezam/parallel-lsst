#ifndef BASICSH_H
#define BASICSH_H

#include <iostream>
#include <lemon/list_graph.h>
#include <limits>
#include <list>

using namespace lemon;
using namespace std;

const double INF = std::numeric_limits<double>::max();
const double ERR = 0.0000000000001;

typedef list<ListGraph::Node> Nlist;

struct NodeStruct {
  long id;
  double dist;
};

struct NodeStructI {
  long id;
  int dist;
};

struct NodePair {
  ListGraph::Node u;
  ListGraph::Node v;
  double dist;
};

bool approx(double a, double b);

// double round(double x);

bool Contains(list<ListGraph::Node> l, ListGraph::Node a);

int compareNS(const void *a, const void *b);

int compareNSI(const void *a, const void *b);

int compareNP(const void *a, const void *b);

double GetRadius(const ListGraph &G, const ListGraph::EdgeMap<double> &length,
                 const ListGraph::Node x0);

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId);

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId,
                     ListGraph::NodeMap<ListGraph::Node> &oldId);

void UsableGraphCopy(const ListGraph &G, ListGraph &G2,
                     const ListGraph::EdgeMap<double> &length,
                     ListGraph::EdgeMap<double> &lengthG2,
                     ListGraph::NodeMap<ListGraph::Node> &newId,
                     ListGraph::NodeMap<ListGraph::Node> &oldId,
                     ListGraph::EdgeMap<ListGraph::Edge> &oldEdgeId);

void GraphSplit(const ListGraph &G, const ListGraph::NodeMap<long> &whComp,
                const long k, const ListGraph::EdgeMap<double> &length,
                ListGraph::EdgeMap<double> **lengthGnew, ListGraph *Gnew,
                ListGraph::NodeMap<ListGraph::Node> &newId,
                ListGraph::NodeMap<ListGraph::Node> **oldId);

void GraphSplit(const ListGraph &G, const ListGraph::NodeMap<long> &whComp,
                const long k, const ListGraph::EdgeMap<double> &length,
                ListGraph::EdgeMap<double> **lengthGnew, ListGraph *Gnew,
                ListGraph::NodeMap<ListGraph::Node> &newId,
                ListGraph::NodeMap<ListGraph::Node> **oldId,
                ListGraph::EdgeMap<ListGraph::Edge> **oldEdgeId);

void Create01Graph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const double dist);

void Create01Graph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const int k);

void CreateStrangeGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                        const int ring, const int side);

void CreateERGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const double p);

void CreateERGraph(ListGraph &G, ListGraph::EdgeMap<double> &length,
                   const long n, const int k);

void CalculateStretches(ListGraph &G, ListGraph::EdgeMap<double> &length,
                        ListGraph &Tree, ListGraph::EdgeMap<double> &Tlength,
                        ListGraph::NodeMap<ListGraph::Node> &TId,
                        double &AvgEStretch, double &MaxEStretch,
                        double &AvgEStretchB, double &MaxEStretchB,
                        double &AvgPStretch, double &MaxPStretch);

#endif // BASICSH_H