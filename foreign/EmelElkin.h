#ifndef EMELELKINH_H
#define EMELELKINH_H

#include <cmath>
#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <list>
#include <map>
using namespace lemon;
struct Decomp {
  long k;
  ListGraph::Node *x;
  ListGraph::Node *y;
};

class EmelElkinBase {
protected:
  double beta;
  long N, M;
  long actN, actM;

  ListGraph &Gr;
  ListGraph::EdgeMap<double> &length;

public:
  EmelElkinBase(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : Gr(G), length(l) {
    N = countNodes(G);
    M = countEdges(G);
    beta =
        1 / (2 * (log((double)N + (double)32) / log(((double)4 / (double)3))));
  }

  void run(ListGraph::Node &x0, ListGraph &Tree,
           ListGraph::EdgeMap<double> &lengthTNew,
           ListGraph::NodeMap<ListGraph::Node> &TNewId,
           ListGraph::NodeMap<ListGraph::Node> &TOldId) {
    LowStrecthTree(Gr, x0, length, Tree, lengthTNew, TNewId, TOldId);
  }

  void GetConeDistance(const ListGraph &G,
                       const ListGraph::EdgeMap<double> &length,
                       const ListGraph::NodeMap<bool> &set,
                       const ListGraph::Node x0,
                       ListGraph::NodeMap<double> &dist);

  double BallCut(const ListGraph &G, const ListGraph::Node x0,
                 const ListGraph::EdgeMap<double> &length, const double ro,
                 const double delta, ListGraph::NodeMap<bool> &Ball,
                 ListGraph::NodeMap<bool> &Shell);

  virtual double ConeCut(const ListGraph &G, const ListGraph::Node x0,
                         const ListGraph::EdgeMap<double> &length,
                         const double lambda1, const double lambda2,
                         const ListGraph::NodeMap<bool> &S,
                         ListGraph::NodeMap<bool> &Cone);

  virtual Decomp StarDecomp(const ListGraph &G, ListGraph::Node x0,
                            const ListGraph::EdgeMap<double> &length,
                            const double delta, const double eps,
                            ListGraph::NodeMap<long> &whComp) = 0;

  void LowStrecthTree(ListGraph &G, ListGraph::Node x0,
                      const ListGraph::EdgeMap<double> &length, ListGraph &Tree,
                      ListGraph::EdgeMap<double> &lengthTNew,
                      ListGraph::NodeMap<ListGraph::Node> &TNewId,
                      ListGraph::NodeMap<ListGraph::Node> &TOldId);
};

class EmelElkin : public EmelElkinBase {

public:
  EmelElkin(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : EmelElkinBase(G, l) {}

  Decomp ConeDecomp(const ListGraph &G,
                    const ListGraph::EdgeMap<double> &length,
                    const ListGraph::NodeMap<bool> &S, const double delta,
                    ListGraph::NodeMap<long> &whComp);

  Decomp StarDecomp(const ListGraph &G, ListGraph::Node x0,
                    const ListGraph::EdgeMap<double> &length,
                    const double delta, const double eps,
                    ListGraph::NodeMap<long> &whComp);
};

class EmelElkinPlus : public EmelElkinBase {

public:
  EmelElkinPlus(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : EmelElkinBase(G, l) {}

  Decomp ImpConeDecomp(const ListGraph &G,
                       const ListGraph::EdgeMap<double> &length,
                       const ListGraph::NodeMap<bool> &S, const double delta,
                       const long t, ListGraph::NodeMap<long> &whComp,
                       ListGraph::Node &x1, bool givenX1 = false);

  Decomp StarDecomp(const ListGraph &G, ListGraph::Node x0,
                    const ListGraph::EdgeMap<double> &length,
                    const double delta, const double eps,
                    ListGraph::NodeMap<long> &whComp);
};

class EmelElkinCorr : public EmelElkinBase {
public:
  EmelElkinCorr(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : EmelElkinBase(G, l) {}

  Decomp ConeDecomp(const ListGraph &G,
                    const ListGraph::EdgeMap<double> &length,
                    const ListGraph::NodeMap<bool> &S, const double delta,
                    const ListGraph::NodeMap<double> &distV,
                    ListGraph::NodeMap<long> &whComp);

  Decomp StarDecomp(const ListGraph &G, ListGraph::Node x0,
                    const ListGraph::EdgeMap<double> &length,
                    const double delta, const double eps,
                    ListGraph::NodeMap<long> &whComp);

  double ConeCut(const ListGraph &G, const ListGraph::Node x0,
                 const ListGraph::EdgeMap<double> &length, const double lambda1,
                 const double lambda2, const ListGraph::NodeMap<bool> &S,
                 const ListGraph::NodeMap<double> &distV,
                 ListGraph::NodeMap<bool> &Cone);
};

class EmelElkinCorrPlus : public EmelElkinBase {
public:
  EmelElkinCorrPlus(ListGraph &G, ListGraph::EdgeMap<double> &l)
      : EmelElkinBase(G, l) {}

  Decomp ImpConeDecomp(const ListGraph &G,
                       const ListGraph::EdgeMap<double> &length,
                       const ListGraph::NodeMap<bool> &S, const double delta,
                       const long t, const ListGraph::NodeMap<double> &distV,
                       ListGraph::NodeMap<long> &whComp, ListGraph::Node &x1,
                       bool givenX1 = false);

  Decomp StarDecomp(const ListGraph &G, ListGraph::Node x0,
                    const ListGraph::EdgeMap<double> &length,
                    const double delta, const double eps,
                    ListGraph::NodeMap<long> &whComp);

  double ConeCut(const ListGraph &G, const ListGraph::Node x0,
                 const ListGraph::EdgeMap<double> &length, const double lambda1,
                 const double lambda2, const ListGraph::NodeMap<bool> &S,
                 const ListGraph::NodeMap<double> &distV,
                 ListGraph::NodeMap<bool> &Cone);
};

#endif // EMELELKINH_H