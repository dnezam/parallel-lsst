#ifndef DECOMPOSEH_H
#define DECOMPOSEH_H

#include <lemon/list_graph.h>

#include "basics.h"

using namespace lemon;

void SelectCluster(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                   const ListGraph::EdgeMap<int> &whichE, int lenE, int Jj,
                   double x, ListGraph::EdgeMap<bool> &Edges,
                   ListGraph::NodeMap<bool> &Nodes);

void AK_Decompose(ListGraph &G, const ListGraph::EdgeMap<double> &length,
                  ListGraph &Tree, ListGraph::EdgeMap<double> &lengthTNew,
                  ListGraph::NodeMap<ListGraph::Node> &TNewId,
                  ListGraph::NodeMap<ListGraph::Node> &TOldId);

void test4();

#endif // DECOMPOSEH_H