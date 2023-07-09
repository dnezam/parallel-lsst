// Utilities for creating main files that read graphs.
#pragma once

#include "assert.h"
#include "graph_io.h"

#include <chrono>
#include <filesystem>
#include <time.h>

#define run_app(G, APP, mutates, rounds)                                       \
  double total_time = 0.0;                                                     \
  for (size_t r = 0; r < rounds; r++) {                                        \
    if (mutates) {                                                             \
      auto G_copy = G;                                                         \
      total_time += APP(G_copy, P);                                            \
    } else {                                                                   \
      total_time += APP(G, P);                                                 \
    }                                                                          \
  }                                                                            \
  auto time_per_iter = total_time / rounds;                                    \
  std::cout << "# time per iter: " << time_per_iter << "\n";

/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_main(APP, mutates)                                        \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_once_main(APP, mutates)                                   \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * either symmetric or asymmetric graph inputs */
#define generate_main(APP, mutates)                                            \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool binary = P.getOptionValue("-b");                                      \
    bool mmap = P.getOptionValue("-m");                                        \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only asymmetric graph inputs */
#define generate_asymmetric_main(APP, mutates)                                 \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(   \
          iFile, mmap);                                                        \
      run_app(G, APP, mutates, rounds)                                         \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,    \
                                                               binary);        \
      run_app(G, APP, mutates, rounds)                                         \
    }                                                                          \
  }

#define generate_custom_main(APP)                                              \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, "<gbbs_input_graph> <gbbs_output_dir>");   \
    char *iFile = P.getArgument(1);                                            \
    char *oDir = P.getArgument(0);                                             \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    std::string output_base{std::string{oDir} + #APP + std::string{"_"}};      \
    auto gbbs_g =                                                              \
        gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, false, false);   \
    std::vector<double> vec_times{};                                           \
    int total_output = 0;                                                      \
    for (size_t r = 0; r < rounds; r++) {                                      \
      total_output += 1;                                                       \
      auto G_copy = gbbs_g;                                                    \
      auto start_t = std::chrono::high_resolution_clock::now();                \
      auto G_seq = APP(G_copy, P);                                             \
      auto duration = std::chrono::high_resolution_clock::now() - start_t;     \
      auto duration_microseconds =                                             \
          std::chrono::duration_cast<std::chrono::microseconds>(duration);     \
      vec_times.push_back(duration_microseconds.count());                      \
      auto G_out = gbbs::gbbs_io::edge_list_to_symmetric_graph(G_seq);         \
      std::string output_file{output_base + std::to_string(r)};                \
      gbbs::gbbs_io::write_graph_to_file(output_file.c_str(), G_out);          \
    }                                                                          \
    std::cout << "#Total#"                                                     \
              << " " << total_output << std::endl;                             \
    std::cout << "#Time#";                                                     \
    for (auto time : vec_times) {                                              \
      std::cout << " " << time;                                                \
    }                                                                          \
    std::cout << std::endl;                                                    \
  }

/* Run all tests in test_inputs */
#define generate_testing_main(APP)                                             \
  int main(int argc, char *argv[]) {                                           \
    namespace fs = std::filesystem;                                            \
    size_t rounds = 50;                                                        \
                                                                               \
    for (const fs::directory_entry &dir_entry :                                \
         fs::recursive_directory_iterator("test_inputs")) {                    \
      if (!(fs::is_regular_file(dir_entry))) {                                 \
        continue;                                                              \
      }                                                                        \
                                                                               \
      std::string input{dir_entry.path().string()};                            \
      std::cout << "Testing: " << input << std::endl;                          \
                                                                               \
      bool is_span = true;                                                     \
      auto Gs = gbbs::gbbs_io::read_unweighted_symmetric_graph_custom(input);  \
      auto gbbs_g{std::move(std::get<1>(Gs))};                                 \
                                                                               \
      for (size_t r = 0; r < rounds; r++) {                                    \
        std::cout << "Round: " << r << "/" << rounds - 1 << std::endl;         \
                                                                               \
        auto G_copy = gbbs_g;                                                  \
        auto G_seq = APP(G_copy);                                              \
        auto G_out = gbbs::gbbs_io::edge_list_to_symmetric_graph(G_seq);       \
        auto G_out_nx = gbbs::gbbs_io::gbbs_2_nk(G_out);                       \
        is_span =                                                              \
            (is_span && gbbs::gbbs_io::is_spanner(std::get<0>(Gs), G_out_nx)); \
        if (!is_span) {                                                        \
          std::cout << "NOT A SPANNER - Stopping." << std::endl;               \
          return -1;                                                           \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    std::cout << "Testing finished successfully." << std::endl;                \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_main(APP, mutates)                                                    \
  int main(int argc, char *argv[]) {                                                             \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                                           \
    char *iFile = P.getArgument(0);                                                              \
    bool symmetric = P.getOptionValue("-s");                                                     \
    bool compressed = P.getOptionValue("-c");                                                    \
    bool mmap = P.getOptionValue("-m");                                                          \
    bool binary = P.getOptionValue("-b");                                                        \
    if (!symmetric) {                                                                            \
      std::cout                                                                                  \
          << "# The application expects the input graph to be symmetric (-s "                    \
             "flag)."                                                                            \
          << std::endl;                                                                          \
      std::cout << "# Please run on a symmetric input." << std::endl;                            \
    }                                                                                            \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                                          \
    if (compressed) {                                                                            \
      std::cout << "Compressed graphs were disabled in benchmark.h"                              \
                << std::endl;                                                                    \
      std::cout << "# Please run on a uncompressed input." << std::endl;                         \
                                                                                                 \
      /*auto G =                                                                                 \
       * gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(*/ /*      iFile, mmap); */ \
      /*  run_app(G, APP, mutates, rounds)  */                                                   \
    } else {                                                                                     \
      auto G =                                                                                   \
          gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, binary);                   \
      run_app(G, APP, mutates, rounds)                                                           \
    }                                                                                            \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_once_main(APP, mutates)                             \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    if (!symmetric) {                                                          \
      std::cout                                                                \
          << "# The application expects the input graph to be symmetric (-s "  \
             "flag)."                                                          \
          << std::endl;                                                        \
      std::cout << "# Please run on a symmetric input." << std::endl;          \
    }                                                                          \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
          iFile, mmap);                                                        \
      run_app(G, APP, mutates, 1)                                              \
    } else {                                                                   \
      auto G =                                                                 \
          gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, binary); \
      run_app(G, APP, mutates, 1)                                              \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_weighted_main(APP, mutates)                                   \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>(   \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::intE>(  \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(     \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<gbbs::intE>(    \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_float_main(APP, mutates)                                      \
  int(int argc, char *argv[]) {                                                \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<float>(iFile,  \
                                                                       mmap);  \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<float>(iFile, \
                                                                        mmap); \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(          \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<float>(         \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * only symmetric graph inputs */
#define generate_symmetric_weighted_main(APP, mutates)                         \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric););        \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>(     \
          iFile, mmap);                                                        \
      run_app(G, APP, mutates, rounds)                                         \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(       \
          iFile, mmap, binary);                                                \
      run_app(G, APP, mutates, rounds)                                         \
    }                                                                          \
  }

/* Macro to generate binary for floating-point weighted graph applications that
 * can ingest only symmetric graph inputs */
#define generate_symmetric_float_weighted_main(APP)                            \
  int main(int argc, char *argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char *iFile = P.getArgument(0);                                            \
    debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric););        \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      ABORT("Graph compression not yet implemented for float weights");        \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(            \
          iFile, mmap, binary);                                                \
      run_app(G, APP, mutates, rounds)                                         \
    }                                                                          \
  }
