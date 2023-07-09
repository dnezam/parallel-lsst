#include <chrono>
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

#include "AbrahamBartal_seq.h"
#include "EmelElkin.h"
#include "basics.h"
#include "decompose.h"

#include "gbbs/gbbs.h"

using namespace lemon;
using namespace std;
using namespace std::chrono;

void testAlgorithms(int Algorithm, int GraphModel) {
  if (Algorithm < 1 || Algorithm > 4 || GraphModel < 1 || GraphModel > 4)
    return;

  int nvals[11] = {100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400, 3200};
  int kvals[4] = {6, 8, 10, 12};
  double pvals[4] = {0.02, 0.05, 0.1, 0.2};
  double dvals[4] = {0.02, 0.05, 0.1, 0.2};

  FILE *out;
  if ((out = fopen("./outf.txt", "wt")) == NULL) {
    cout << "Unable 2 create output file";
    return;
  }

  int runs = 5;

  int i = 0;
  while (i < 4) {
    auto tmp = high_resolution_clock::now();
    auto fullDuration = duration_cast<microseconds>(tmp - tmp);
    // cout << "-----------DURATION CHECK-----------" << std::endl;
    // cout << fullDuration.count() << endl;
    // cout << "----------------------" << std::endl;

    int j = 0;
    while (j < 11) {
      int nval = nvals[j];
      double p;
      int kval;
      if (GraphModel == 1)
        p = pvals[i];
      if (GraphModel == 2 || GraphModel == 4)
        kval = kvals[i];
      if (GraphModel == 3)
        p = dvals[i];

      double sums[7];
      int k = 0;
      while (k < 7) {
        sums[k] = 0;
        ++k;
      }
      k = 0;
      while (k != runs) {
        ListGraph G;
        ListGraph::EdgeMap<double> length(G);
        if (GraphModel == 1)
          CreateERGraph(G, length, nval, p);
        if (GraphModel == 2)
          CreateERGraph(G, length, nval, kval);
        if (GraphModel == 3)
          Create01Graph(G, length, nval, p);
        if (GraphModel == 4)
          Create01Graph(G, length, nval, kval);

        gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::empty> gbbs_graph{};
        // for (auto e: length) {
        // 	cout << "HI" << std::endl;
        // }
        cout << "Nodes: " << countNodes(G) << std::endl;
        cout << "Edges: " << countEdges(G) << std::endl;
        cout << "Perc: "
             << 2 * (double)countEdges(G) /
                    (double)((countNodes(G) - 1) * countNodes(G))
             << std::endl
             << std::endl;

        sums[0] += countEdges(G);

        ListGraph Tree;
        ListGraph::EdgeMap<double> lengthNew(Tree);

        ListGraph::NodeMap<ListGraph::Node> TNewId(G);
        ListGraph::NodeMap<ListGraph::Node> TOldId(Tree);

        ListGraph::Node x0;
        for (ListGraph::NodeIt i(G); i != INVALID; ++i) {
          x0 = i;
          break;
        }

        double AvgEStretch, MaxEStretch, AvgEStretchB, MaxEStretchB,
            AvgPStretch, MaxPStretch;

        if (Algorithm == 1) {
          EmelElkin ee(G, length);
          ee.run(x0, Tree, lengthNew, TNewId, TOldId);
        }
        if (Algorithm == 2) {
          AK_Decompose(G, length, Tree, lengthNew, TNewId, TOldId);
        }
        if (Algorithm == 3) {
          cout << "Parallel algorithm: with i: " << i << ", j: " << j
               << std::endl;

          auto start2 = high_resolution_clock::now();
          AbrahamBartal_seq ab(G, length);
          ab.run(x0, Tree, lengthNew, TNewId, TOldId);
          auto stop2 = high_resolution_clock::now();

          auto duration = duration_cast<microseconds>(stop2 - start2);
          fullDuration = duration_cast<microseconds>(fullDuration + duration);
          // cout << "-TIME-" << std::endl;
          // cout << duration.count() << endl;
          // cout << "-TIME-" << std::endl;
        }
        if (Algorithm == 4) {
          ListGraph::EdgeMap<ListGraph::Edge> oldEdgeId(Tree);
          ListGraph::EdgeMap<bool> inmin(G);
          UsableGraphCopy(G, Tree, length, lengthNew, TNewId, TOldId,
                          oldEdgeId);
          kruskal(G, length, inmin);

          std::list<ListGraph::Edge> toDelete;
          for (ListGraph::EdgeIt e(Tree); e != INVALID; ++e) {
            ListGraph::Edge edg = e;
            if (!inmin[oldEdgeId[edg]])
              toDelete.push_back(edg);
          }
          while (!toDelete.empty()) {
            Tree.erase(toDelete.front());
            toDelete.pop_front();
          }
        }

        cout << "Nodes: " << countNodes(Tree) << std::endl;
        cout << "Edges: " << countEdges(Tree) << std::endl;
        if (!connected(Tree))
          cout << "not ";
        cout << "connected. " << std::endl << std::endl;

        CalculateStretches(G, length, Tree, lengthNew, TNewId, AvgEStretch,
                           MaxEStretch, AvgEStretchB, MaxEStretchB, AvgPStretch,
                           MaxPStretch);

        cout << "Avg Edge Stretch: " << AvgEStretch << std::endl;
        cout << "Max Edge Stretch: " << MaxEStretch << std::endl;
        cout << "Avg Edge StretchB: " << AvgEStretchB << std::endl;
        cout << "Max Edge StretchB: " << MaxEStretchB << std::endl;
        cout << "Avg Path Stretch: " << AvgPStretch << std::endl;
        cout << "Max Path Stretch: " << MaxPStretch << std::endl;
        sums[1] += AvgEStretch;
        sums[4] += MaxEStretch;
        sums[2] += AvgEStretchB;
        sums[5] += MaxEStretchB;
        sums[3] += AvgPStretch;
        sums[6] += MaxPStretch;

        Tree.clear();

        ++k;
      }
      fprintf(out, "%d\t%d\t\t", Algorithm, GraphModel);
      if (GraphModel == 2 || GraphModel == 4)
        fprintf(out, "%d\t%d\t\t", nval, kval);
      else
        fprintf(out, "%d\t%lf\t\t", nval, p);
      k = 0;
      while (k < 7) {
        fprintf(out, "%lf\t", sums[k] / runs);
        if (k == 0 || k == 3)
          fprintf(out, "\t");
        ++k;
      }
      fprintf(out, "\n");

      ++j;
    }
    fprintf(out, "\n\n");

    // Calculating time
    cout << "-----------TIME-----------" << std::endl;
    cout << fullDuration.count() << endl;
    cout << "-----------TIME-----------" << std::endl;

    ++i;
  }

  fclose(out);
}

int main() {
  srand(0);
  // testStarDecomp(4, 4);
  // 3 -> the parallel graph implementation
  // 2 -> sparse Erdos Renyi random graph as input
  testAlgorithms(3, 2);
  getchar();
  return 0;
}