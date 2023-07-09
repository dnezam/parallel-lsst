#ifndef ABRAHAMBARTALH_H
#define ABRAHAMBARTALH_H

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <list>

#include "EmelElkin.h"
#include "basics.h"

using namespace std;
struct QDecomp {
  long k;
  ListGraph::Node *x;
  ListGraph::Node *y;
  list<ListGraph::Node> *Q;
};

class AbrahamBartal_seq {
  double c;
  double beta;
  long N;

  ListGraph &Gr;
  ListGraph::EdgeMap<double> &length;

public:
  AbrahamBartal_seq(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : Gr(G), length(l) {
    c = 2;
    N = countNodes(G);
    beta =
        1 / (2 * (log((double)N + (double)32) / log(((double)4 / (double)3))));
  }

  void run(ListGraph::Node &x0, ListGraph &Tree,
           ListGraph::EdgeMap<double> &lengthTNew,
           ListGraph::NodeMap<ListGraph::Node> &TNewId,
           ListGraph::NodeMap<ListGraph::Node> &TOldId) {
    list<ListGraph::Node> Q;
    for (ListGraph::NodeIt i(Gr); i != INVALID; ++i) {
      ListGraph::Node a = i;
      if (a == x0)
        continue;
      Q.push_back(a);
    }
    HierarStarPart(Gr, x0, length, Q, Tree, lengthTNew, TNewId, TOldId);
  }

  void HierarStarPart(ListGraph &G, ListGraph::Node x0,
                      const ListGraph::EdgeMap<double> &length,
                      list<ListGraph::Node> Q, ListGraph &Tree,
                      ListGraph::EdgeMap<double> &lengthTNew,
                      ListGraph::NodeMap<ListGraph::Node> &TNewId,
                      ListGraph::NodeMap<ListGraph::Node> &TOldId);

  Decomp DiffConeDecomp(const ListGraph &G,
                        const ListGraph::EdgeMap<double> &length,
                        const ListGraph::NodeMap<bool> &S, const double delta,
                        ListGraph::NodeMap<long> &whComp);

  QDecomp StarPart(const ListGraph &G, ListGraph::Node x0,
                   const ListGraph::EdgeMap<double> &length,
                   const list<ListGraph::Node> Q,
                   ListGraph::NodeMap<long> &whComp);
};

#endif // ABRAHAMBARTALH_H
