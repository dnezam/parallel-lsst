#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>

#include "basics.h"

void CalculateStretches(ListGraph &G, ListGraph::EdgeMap<double> &length,
                        ListGraph &Tree, ListGraph::EdgeMap<double> &Tlength,
                        ListGraph::NodeMap<ListGraph::Node> &TId,
                        double &AvgEStretch, double &MaxEStretch,
                        double &AvgEStretchB, double &MaxEStretchB,
                        double &AvgPStretch, double &MaxPStretch) {
  AvgEStretch = 0;
  MaxEStretch = 0;
  AvgEStretchB = 0;
  MaxEStretchB = 0;
  AvgPStretch = 0;
  MaxPStretch = 0;

  for (ListGraph::NodeIt u(G); u != INVALID; ++u) {
    ListGraph::Node a = u;

    ListGraph::NodeMap<double> dist(G);

    Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(G, length);
    dijkstra.distMap(dist);
    dijkstra.init();
    dijkstra.addSource(a);
    dijkstra.start();

    ListGraph::NodeMap<double> Tdist(Tree);

    Dijkstra<ListGraph, ListGraph::EdgeMap<double>> Tdijkstra(Tree, Tlength);
    Tdijkstra.distMap(Tdist);
    Tdijkstra.init();
    Tdijkstra.addSource(TId[a]);
    Tdijkstra.start();

    for (ListGraph::IncEdgeIt e(G, a); e != INVALID; ++e) {
      if (G.u(e) == a) {
        double stretch = Tdist[TId[G.v(e)]] / length[e];
        AvgEStretch += stretch;
        if (stretch > MaxEStretch)
          MaxEStretch = stretch;

        stretch = Tdist[TId[G.v(e)]] / dist[G.v(e)];
        AvgEStretchB += stretch;
        if (stretch > MaxEStretchB)
          MaxEStretchB = stretch;
      }
    }

    for (ListGraph::NodeIt v(G); v != INVALID; ++v) {
      ListGraph::Node b = v;
      if (u < v) {
        double stretch = Tdist[TId[b]] / dist[b];
        AvgPStretch += stretch;
        if (stretch > MaxPStretch)
          MaxPStretch = stretch;
      }
    }
  }

  AvgEStretch /= (double)countEdges(G);
  AvgEStretchB /= (double)countEdges(G);
  AvgPStretch /= ((double)countNodes(G) * ((double)countNodes(G) - 1) / 2);
}