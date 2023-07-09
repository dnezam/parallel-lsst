#include "lemon/include/lemon/list_graph.h"
#include "tools/common.hpp"

#define lemon_main(APP)                                                        \
  int main(int argc, char **argv) {                                            \
    gbbs::commandLine P(argc, argv, "<nx_input_graph>");                       \
    char *iFile = P.getArgument(0);                                            \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    auto nx_g = customReader<NetworKit::Graph>(iFile);                         \
    lemon::ListGraph lemon_g;                                                  \
    NX2Lemon(nx_g, lemon_g);                                                   \
    std::vector<double> vec_times{};                                           \
    std::vector<double> vec_stretch{};                                         \
    int total_output = 0;                                                      \
    int spanner_count = 0;                                                     \
    for (size_t r = 0; r < rounds; r++) {                                      \
      total_output += 1;                                                       \
      auto stat = APP(lemon_g);                                                \
      auto is_spanner = std::get<0>(stat);                                     \
      auto time = std::get<1>(stat);                                           \
      auto stretch = std::get<2>(stat);                                        \
      vec_times.push_back(time);                                               \
      if (is_spanner) {                                                        \
        spanner_count += 1;                                                    \
        vec_stretch.push_back(stretch);                                        \
      }                                                                        \
    }                                                                          \
    std::cout << "#Total#"                                                     \
              << " " << total_output << std::endl;                             \
    std::cout << "#Is_Spanner#"                                                \
              << " " << spanner_count << std::endl;                            \
    std::cout << "#Stretch#";                                                  \
    for (double stretch : vec_stretch) {                                       \
      std::cout << " " << stretch;                                             \
    }                                                                          \
    std::cout << std::endl;                                                    \
    std::cout << "#Time#";                                                     \
    for (auto time : vec_times) {                                              \
      std::cout << " " << time;                                                \
    }                                                                          \
    std::cout << std::endl;                                                    \
  }

#define nx_main(APP)                                                           \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, "<nx_input_graph> <nx_output_dir>");       \
    char *iFile = P.getArgument(1);                                            \
    char *oDir = P.getArgument(0);                                             \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    std::string output_base{std::string{oDir} + #APP + std::string{"_"}};      \
    auto nx_g = customReader<NetworKit::Graph>(iFile);                         \
    std::vector<double> vec_times{};                                           \
    int total_output = 0;                                                      \
    for (size_t r = 0; r < rounds; r++) {                                      \
      total_output += 1;                                                       \
      auto G_copy = nx_g;                                                      \
      auto start_t = std::chrono::high_resolution_clock::now();                \
      auto T_out = APP(G_copy);                                                \
      auto duration = std::chrono::high_resolution_clock::now() - start_t;     \
      auto duration_microseconds =                                             \
          std::chrono::duration_cast<std::chrono::microseconds>(duration);     \
      vec_times.push_back(duration_microseconds.count());                      \
      std::string output_file{output_base + std::to_string(r)};                \
      customWriter(T_out, output_file);                                        \
    }                                                                          \
    std::cout << "#Total#"                                                     \
              << " " << total_output << std::endl;                             \
    std::cout << "#Time#";                                                     \
    for (auto time : vec_times) {                                              \
      std::cout << " " << time;                                                \
    }                                                                          \
    std::cout << std::endl;                                                    \
  }
