#ifndef ALANKARPH_H
#define ALANKARPH_H

#include <lemon/core.h>
#include <lemon/list_graph.h>

#include "basics.h"

void SelectCluster(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                   const ListGraph::EdgeMap<int> &whichE, int lenE, int Jj,
                   double x, ListGraph::EdgeMap<bool> &Edges,
                   ListGraph::NodeMap<bool> &Nodes);

void AK_Decompose(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                  ListGraph &Tree, ListGraph::EdgeMap<double> &lengthTNew,
                  ListGraph::NodeMap<ListGraph::Node> &TNewId,
                  ListGraph::NodeMap<ListGraph::Node> &TOldId);

class AlanKarp {
  ListGraph &Gr;
  ListGraph::EdgeMap<double> &length;

public:
  AlanKarp(ListGraph &G, ListGraph::EdgeMap<double> &l) : Gr(G), length(l){};

  void run(ListGraph::Node &x0, ListGraph &Tree,
           ListGraph::EdgeMap<double> &lengthTNew,
           ListGraph::NodeMap<ListGraph::Node> &TNewId,
           ListGraph::NodeMap<ListGraph::Node> &TOldId) {
    AK_Decompose(Gr, length, Tree, lengthTNew, TNewId, TOldId);
  }
};

#endif // ALANKARPH_H