#pragma once

#include "gbbs/gbbs.h"
// Dear all, please pay attention to this include path
// ALWAYS start from the project root(where the WORKSPACE file is defined)
#include "benchmarks/BFS/NonDeterministicBFS/BFS.h"

#include "benchmarks/LowAverageStretchTree_ParallelStarDecomp/BFS_radius_grow.h"

#include "gbbs/contract.h"
#include <iostream>
#include <omp.h>
#include <cmath>

namespace gbbs {
    
    namespace LowAvgStretch_utils {

        constexpr size_t small_cluster_size = 2048;
        using edge = std::tuple<uintE, uintE>;
        using edge_id = std::tuple<uintE, uintE, uintE>;
        using edge_entry = std::tuple<edge, gbbs::empty>;


        template <
                template <class W> class vertex, class W,
                // This means W(weight) is not defined, thus unweighted graph
                typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline parlay::sequence<edge> SplitGraph_ParaStarDecomp(symmetric_graph<vertex, W>& G, double p){
            // Later during the algo I need to have a random permutation, for now I will choose a random number
            // and then start my sequence from that random number
            uint random_starter = 0;
            int n = G.n;
            int initial_n = n;

            using unweighted_edge = std::tuple<uintE, uintE, gbbs::empty>;
            using index = size_t;
            using indexed_edge = std::tuple<uintE, uintE, index>;
            
            auto edge_list = to_edge_array<gbbs::empty>(G);
            //std::cout << "Initial edges: source_dest\n";
            sequence<unweighted_edge> original_edges = G.edges();
            auto initial_edges = sequence<indexed_edge>(G.m);
            auto current_edges = sequence<indexed_edge>(G.m);
            auto solution_ids = sequence<index>();

            size_t original_edges_size = original_edges.size();
            parallel_for(0, original_edges_size, [&](int i) {
                indexed_edge ie = std::make_tuple(std::get<0>(original_edges[i]),
                                                  std::get<1>(original_edges[i]),
                                                  i);
                initial_edges[i] = ie;
                current_edges[i] = ie;
            });
            
            size_t nr_centers = 1;

            int R = (int) ceil((p/(2*log2(n))));
            //std::cout << "\nR is: " << p << "\n";
            // Check loop condition
            sequence<index> sol_ids;
            for(int i = 1; i <= 10;i++){
                n = G.n;
                double t_delta = 12*std::pow(n,i/(2*log2(n)) - 1)*G.n*log2(n);
                
                //int nr_centers = 8;// G.n < t_delta ? G.n : (int)t_delta;
                nr_centers = nr_centers >= G.n ? 1:nr_centers;
                int center_index_start = rand() % n; 
                //std::cout << "there exist: " << nr_centers << " centers and n: " << n << "\n";
                // For each center draw random distance
                std::vector<std::pair<int,int>> delta_s(nr_centers); // center,radius
                
                parallel_for(0, nr_centers, [&](size_t j) {

                    delta_s[j].second = (int)(rand() % (R+1)); // max_radius
                    delta_s[j].first = (center_index_start + j) % n; // center_id
            });


                //std::cout << "\n";
                double radius_general = (2*log2(n)-i+1)*R;
                // For each center compute the ball
                std::vector<std::pair<sequence<uintE>,std::vector<int>>> balls(nr_centers); // Parents, dist

                parallel_for(0, nr_centers, [&](size_t j) {

                    int offset = j;
                    // TODO should be only until radius size
                    //std::cout << "BFS\n";
                    auto res = BFS_grow_rad(G,delta_s[j].first, delta_s[j].second); // parents , distance
                    // res.second holds distance to center
                    balls[j] = res;
                });

                std::vector<sequence<std::tuple<int,uintE,uintE>>> distance_to_center(n, sequence<std::tuple<int,uintE,uintE>>(nr_centers)); // distance,parent,center


                parallel_for(0, nr_centers, [&](size_t j) {
                    parallel_for(0, n, [&](size_t k) {
                    distance_to_center[k][j] = std::make_tuple((balls[j].second)[k],(balls[j].first)[k],(center_index_start + j) % n);
                   
                });
                });
                
                
                std::vector<std::tuple<int,int,uintE>> min_dist(n); // distance, parent , center

                parallel_for(0, n, [&](uintE j) {
                    min_dist[j] = parlay::reduce(distance_to_center[j],parlay::minm<std::tuple<int,uintE,uintE>>());
                    min_dist[j] = std::get<0>(min_dist[j]) == INT_MAX ? std::make_tuple(INT_MAX,INT_MAX,j) : min_dist[j];
                });

                // Assign each vertex that was reached the value that minimizes dist(center,vertex) + delta_s[center]
                sequence<uintE> clusters(n);

                parallel_for(0, n, [&](int j) {
                    clusters[j] = std::get<2>(min_dist[j]);
                    std::get<2>(min_dist[j]) = j;
                });

                auto no_duplicates_clusters = parlay::remove_duplicates_ordered(clusters,std::less<uintE>{});
                size_t num_clusters = no_duplicates_clusters.size();
                //std::cout << "clusters:" << num_clusters << "\n";
                // Contract the graph into its clusters
                size_t n_clusters = contract::RelabelIds(clusters);

                
                //std::cout << "\n";
                auto res_contract = contract::contract(G,clusters,num_clusters);
                //std::cout << "N = " << num_clusters << "\n";
                G = std::move(std::get<0>(res_contract));
            sequence<std::tuple<int,uintE,uintE>> sol;

            size_t until = min_dist.size();


            
            auto filtered_edges = parlay::filter(min_dist, [&](const std::tuple<int,int,uintE>& e) {
                    return std::get<0>(e) > 0 && std::get<0>(e) < INT_MAX;
                });

            sequence<index> new_edges = parlay::map(filtered_edges, [&](const std::tuple<int,int,uintE>& e) {
                    const auto& m = *std::lower_bound(current_edges.cbegin(), current_edges.cend(),
                                                    std::make_tuple(std::get<2>(e), std::get<1>(e), 0));
                    return std::get<2>(m);
                });

            
            sol_ids = parlay::append(sol_ids, new_edges);
            
        
            // Prepare current_edges for next iteration
                current_edges = parlay::map(current_edges, [&](const indexed_edge& ie) {
                    uintE u = clusters[std::get<0>(ie)];
                    uintE v = clusters[std::get<1>(ie)];
                    if (u == v) {
                        return std::make_tuple(UINT_E_MAX, UINT_E_MAX, SIZE_MAX);
                    } else {
                        return std::make_tuple(u, v, std::get<2>(ie));
                    }
                });

                current_edges = parlay::filter(current_edges, [&](const indexed_edge& ie) {
                    return std::get<2>(ie) != SIZE_MAX;
                });

                current_edges = parlay::sort(current_edges);


            // For completely reduced graphs we can terminate early
            if(1 >= G.n){
                break;
            }
            }
            
            sequence<edge> ret_val(sol_ids.size());

            parallel_for(0, sol_ids.size(), [&](size_t i) {
                    ret_val[i] = std::make_tuple(std::get<0>(initial_edges[sol_ids[i]]),std::get<1>(initial_edges[sol_ids[i]]));

                });



            return ret_val;
        }

    template <
                template <class W> class vertex, class W,
                // This means W(weight) is not defined, thus unweighted graph
                typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline sequence<edge> Partition_ParaStarDecomp(symmetric_graph<vertex, W>& G, double p){
            int k = 1; // nr_edge_classes
            int n = G.n;
            auto res = SplitGraph_ParaStarDecomp(G,p);
            double threshold = G.num_edges() * 272 * k * pow(log(n),3) / p;
            // If there exists component that has more than threshold edges on its border, redo SplitGraph
            // Since an unweighted graph has only a single component, this will always hold
            std::cout << "size result: " << res.size() << "\n";
            return res;
        }
    }

    namespace LowAvgStretch {
        
        using edge = std::tuple<uintE, uintE>;
        template <
                template <class W> class vertex, class W,
                // This means W(weight) is not defined, thus unweighted graph
                typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline sequence<edge> ParaStarDecomp(symmetric_graph<vertex, W>& G) {
            // Unweighted Algo here
            
            // 1) Normalize edgeset -> not necessary, unweighted graph
            // 2) Bunch of definitions
            double c1 = 272; // given in the paper
            int n = G.n;
            double y = std::pow(2,sqrt(6*log2(n)*log2(log2(n))));
            double theta = ceil(3*log2(n)/log2(y));
            double z = 4*y*theta*std::pow(log2(n),3) * c1;
            //std::cout << z << std::endl;
            // 3) Create subclasses (k) of edges depending on weight
            // for unweighted graphs all edges fall into first subclass
            /*
            4) While(G not single vertex){
                (C1,...) = Partition(V^j,....)
                T = T u BFS(Ci)
                Contract edges within Components
                remove selfloops, keep parallel edges
            }
            return T
            */

            // Testing to contract a graph into two seperate clusters

            auto res_func = LowAvgStretch_utils::Partition_ParaStarDecomp(G,z/4);
            //std::cout << "gets out of part\n";
            //sequence<std::tuple<uintE, uintE, W>> res_func_hand(res_func.size());

            //edge_array<W> res = edge_array<W>(std::move(res_func_hand),n);
            
            return res_func;
        }

        template <
                template <class W> class vertex, class W,
                // This means W(weight) is not undefined, thus weighted graph
                typename std::enable_if<!std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline sequence<edge> ParaStarDecomp(symmetric_graph<vertex, W>& GA) {
            // Weighted Algo here

            return sequence<edge>();
        }

    }
}