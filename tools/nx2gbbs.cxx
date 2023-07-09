#include "tools/common.hpp"

int main(int argc, char **argv) {
  gbbs::commandLine P(argc, argv, " <in_nx_file> <out_gbbs_file>");
  char *iFile = P.getArgument(1);
  char *oFile = P.getArgument(0);
  NetworKit::Graph g = customReader<NetworKit::Graph>(iFile);
  gbbs_graph g_target = NX2gbbs(g);
  customWriter(g_target, oFile);
  return 0;
}