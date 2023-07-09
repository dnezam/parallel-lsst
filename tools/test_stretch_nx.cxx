#include "tools/common.hpp"
#include <filesystem>

int main(int argc, char **argv) {
  gbbs::commandLine P(argc, argv, "<nx_in_file> <nx_out_dir>");
  char *graph_file = P.getArgument(1);
  char *tree_dir = P.getArgument(0);
  NetworKit::Graph nx_g = customReader<NetworKit::Graph>(graph_file);
  int total_output = 0;
  int is_spanner_output = 0;
  std::vector<double> vec_stretch{};
  auto files = get_files(tree_dir);
  for (const auto &tree_file : files) {
    total_output += 1;
    std::string tree_path{tree_dir + tree_file};
    NetworKit::Graph nx_tree = customReader<NetworKit::Graph>(tree_path);
    if (is_spanner(nx_g, nx_tree)) {
      is_spanner_output += 1;
      double stretch = get_averge_stretch(nx_g, nx_tree);
      vec_stretch.push_back(stretch);
    }
  }
  std::cout << "#Total#"
            << " " << total_output << std::endl;
  std::cout << "#Is_Spanner#"
            << " " << is_spanner_output << std::endl;
  std::cout << "#Stretch#";
  for (double stretch : vec_stretch) {
    std::cout << " " << stretch;
  }
  std::cout << std::endl;
}