#pragma once

#include "gbbs/gbbs.h"

#include "benchmarks/BFS/NonDeterministicBFS/BFS.h"

#include <cmath>
#include <random>
#include <tuple>

#define SSSP_CUTOFF 0

namespace gbbs
{
    namespace LowAvgStretch
    {
        template <class T>
        inline bool contains(sequence<T> &seq, T elem)
        {
            for (int i = 0; i < seq.size(); i++)
            {
                if (elem == seq[i])
                {
                    return true;
                }
            }
            return false;
        }
        inline bool printIfFalse(const bool assertion, const char *message)
        {
            if (!assertion)
            {
                std::cout << message << std::endl;
            }
            return assertion;
        }

        inline bool isReachable(sequence<uintE> &parents, uintE dest)
        {
            return (parents[dest] != UINT_E_MAX);
        }

        template <class Graph>
        inline bool isConnected(Graph &G, uintE src)
        {
            for (int i = 0; i < G.n; i++)
            {
                sequence<uint> pars = BFS(G, src);
                for (int j = 0; j < G.n; j++)
                {
                    if (j == i)
                    {
                        continue;
                    }
                    bool reach = isReachable(pars, j);
                    // debug(assert(reach););
                    if (!reach)
                    {
                        return false;
                    }
                }
            }
        }
        template <class Graph>
        inline bool isDisjoint(Graph &G, sequence<uintE> &srcs)
        {
            for (int i = 0; i < srcs.size(); i++)
            {
                sequence<uintE> pars = BFS(G, srcs[i]);
                for (int j = 0; j < srcs.size(); j++)
                {
                    if (j == i)
                    {
                        continue;
                    }
                    bool reach = isReachable(pars, srcs[j]);
                    // debug(assert(!reach););
                    if (reach)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        template <class Graph>
        inline void isBallCut(Graph &G, sequence<uintE> &ball, sequence<uintE> &ballShell)
        {
            for (int i = 0; i < ball.size(); i++)
            {
                auto pars = BFS(G, ball[i]);
                for (int j = 0; j < ballShell.size(); j++)
                {
                    // debug(assert(!isReachable(pars, ballShell[j])););
                }
            }
        }

        template <class Graph>
        inline void correctCluster(Graph &G, sequence<uintE> &srcs, sequence<uintE> &cluster)
        {
            char buff[256];
            // // debug(std::cout << "correctCluster with srcs.size(): " << srcs.size() << std::endl;);
            for (int i = 0; i < srcs.size(); i++)
            {
                sequence<uintE> pars = BFS(G, srcs[i]);
                for (uintE j = 0; j < G.n; j++)
                {
                    sprintf(buff, "\ni: %d, j: %d, cluster[srcs[i]]: %d, cluster[j]: %d, src: %d\n", i, j, cluster[srcs[i]], cluster[j], srcs[i]);
                    if (isReachable(pars, j))
                    {
                        // debug(assert(printIfFalse(cluster[srcs[i]] == cluster[j], buff)););
                    }
                    else
                    {
                        // debug(assert(printIfFalse(cluster[srcs[i]] != cluster[j], buff)););
                    }
                }
            }
        }

        template <class T>
        inline bool notContained(sequence<T> &seq, T elem)
        {
            for (int i = 0; i < seq.size(); i++)
            {
                bool eq = seq[i] == elem;
                // debug(assert(!eq););
                if (eq)
                {
                    return false;
                }
            }
            return true;
        }

        template <class T>
        inline bool inRange(sequence<T> &seq, sequence<T> offset, double lower, double upper)
        {
            for (int i = 0; i < seq.size(); i++)
            {
                bool in = (floor(offset[i] * lower) <= seq[i] && seq[i] <= floor(offset[i] * upper));
                // debug(assert(in););
                if (!in)
                {
                    return false;
                }
            }
            return true;
        }

        inline void correctClusterIndex(sequence<size_t> &index, sequence<uintE> &srcs)
        {
            char buff[256];
            for (int i = 0; i < index.size(); i++)
            {
                if (index[i] != SIZE_MAX)
                {
                    sprintf(buff, "\ni: %d, index[i]: %d, srcs[index[i]]: %d\n", i, index[i], srcs[index[i]]);
                    // debug(assert(printIfFalse(srcs[index[i]] == i, buff)););
                }
            }
        }

        inline void noDoubleEdges(sequence<std::tuple<uintE, uintE, gbbs::empty>> &edges)
        {
            uintE u1, u2;
            uintE v1, v2;
            gbbs::empty w1, w2;
            for (int i = 0; i < edges.size(); i++)
            {
                std::tie(u1, v1, w1) = edges[i];
                for (int j = i + 1; j < edges.size(); j++)
                {
                    std::tie(u2, v2, w2) = edges[j];
                    // debug(assert(!(((u1 == u2) && (v1 == v2)) || ((u1 == v2) && (v1 == u2)))););
                }
            }
        }

        typedef std::pair<uintE, size_t> PPD;
        PPD defaultPPD = std::make_pair(UINT_E_MAX, SIZE_MAX);

        template <class W>
        struct BFS_R
        {
            // uintE *Parents;
            // size_t *dists;
            PPD *Parents;
            size_t dist;
            BFS_R(PPD *_Parents, size_t _dist) : Parents(_Parents)
            {
                dist = _dist;
            }
            inline bool update(uintE s, uintE d, W w)
            {
                // std::cout << "updating " << d << " with parent " << s << ", and dist " << dist << std::endl;
                if (Parents[d].first == UINT_E_MAX && Parents[d].second == SIZE_MAX)
                {
                    Parents[d] = std::make_pair(s, dist);
                    return 1;
                }
                else
                {
                    return 0;
                }
            }
            inline bool updateAtomic(uintE s, uintE d, W w)
            {
                // std::cout << "updating atomic " << d << " with parent " << s << ", and dist " << dist << std::endl;
                // std::cout << "Address: " << &Parents[d] << std::endl;
                return (gbbs::atomic_compare_and_swap((&Parents[d]), {UINT_E_MAX, SIZE_MAX}, std::make_pair(s, Parents[s].second)));
            }
            inline bool cond(uintE d) { return (Parents[d].first == UINT_E_MAX && Parents[d].second == SIZE_MAX); }
        };
        template <class W>
        struct BFS_dist
        {
            size_t *dists;
            BFS_dist(size_t *_Parents) : dists(_Parents)
            {
            }
            inline bool update(uintE s, uintE d, W w)
            {
                // std::cout << "updating " << d << " with parent " << s << ", and dist " << dist << std::endl;
                if (dists[d] == SIZE_MAX)
                {
                    dists[d] = dists[s] + 1;
                    return 1;
                }
                else
                {
                    return 0;
                }
            }
            inline bool updateAtomic(uintE s, uintE d, W w)
            {
                return (gbbs::atomic_compare_and_swap(&dists[d], SIZE_MAX, dists[s] + 1));
            }
            inline bool cond(uintE d) { return (dists[d] == SIZE_MAX); }
        };

        template <class W>
        struct BFS_Parents
        {
            uintE *Parents;
            BFS_Parents(uintE *_Parents) : Parents(_Parents) {}
            inline bool update(uintE s, uintE d, W w)
            {
                if (Parents[d] == UINT_E_MAX)
                {
                    Parents[d] = s;
                    return 1;
                }
                else
                {
                    return 0;
                }
            }
            inline bool updateAtomic(uintE s, uintE d, W w)
            {
                return (gbbs::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, s));
            }
            inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
        };

        template <class Graph>
        inline void BFS_radius(Graph &G,
                               sequence<uintE> &srcs,
                               sequence<size_t> &radius,
                               sequence<uintE> &parents,
                               sequence<uintE> &cluster,
                               sequence<size_t> &cluster_index,
                               sequence<size_t> &dists,
                               size_t &max_rad,
                               sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> &edges) // sequence<PPD> &parents)
        {
            // debug(std::cout << "entering BFS_radius\n";);
            // debug(assert(!srcs.empty()););
            // debug(assert(cluster.size() == G.n););
            // debug(assert(cluster_index.size() == G.n););
            using W = typename Graph::weight_type;
            // debug(auto srcsSize_tmp = srcs.size(););
            // debug(auto edgesSize_tmp = edges.size(););
            // debug(correctCluster(G, srcs, cluster););
            // debug(correctClusterIndex(cluster_index, srcs););
            // debug(noDoubleEdges(edges););
            // parents = sequence<PPD>::uninitialized(G.n);
            // parallel_for(0, parents.size(), [&](size_t i)
            //              { parents[i] = std::make_pair(UINT_E_MAX, SIZE_MAX); });
            // parents[src] = std::make_pair(src, 0);

            // std::cout << "srcs.size(): " << srcs.size() << std::endl;
            parents = sequence<uintE>::from_function(G.n, [&](size_t i)
                                                     { return UINT_E_MAX; });
            dists = sequence<size_t>::from_function(G.n, [&](size_t i)
                                                    { return SIZE_MAX; });
            radius = sequence<size_t>::from_function(srcs.size(), [&](size_t i)
                                                     { return 0; });

            parallel_for(0, srcs.size(), [&](size_t i)
                         {dists[srcs[i]]=0; parents[srcs[i]]=srcs[i]; });

            vertexSubset Frontier(G.n);
            add_to_vsubset(Frontier, srcs.begin(), srcs.size());

            max_rad = 1;
            while (!Frontier.isEmpty())
            {
                
                edgeMap(G, Frontier, BFS_dist<W>(dists.begin()), -1, sparse_blocked | dense_parallel);
                Frontier = edgeMap(G, Frontier, BFS_Parents<W>(parents.begin()), -1,
                                   sparse_blocked | dense_parallel);

                for (uintE i = 0; i < G.n; i++) {
                    if (dists[i] == max_rad) {
                        radius[cluster_index[i]] = max_rad;
                    }
                }

                max_rad++;
            }
            max_rad--;


            // debug(assert(notContained(radius, SIZE_MAX)););
            // TODO: remove if radius <= Cutoff

            // auto new_edges = sequence<std::tuple<uintE, uintE, typename Graph::weight_type>>::from_function(G.n, [&](uintE i)
            //                                                                                                 {
            //     if (radius[cluster_index[cluster[i]]] > SSSP_CUTOFF) {
            //         return std::make_tuple(UINT_E_MAX, UINT_E_MAX, gbbs::empty());
            //     } else {
            //         return std::make_tuple(i, parents[i], gbbs::empty());
            //     } });

            // new_edges = parlay::remove_if(new_edges, [&](auto e)
            //                               { return ((std::get<0>(e) != std::get<1>(e)) && (std::get<0>(e) != UINT_E_MAX) && (std::get<1>(e) != UINT_E_MAX)); });

            // // debug(std::cout << "Appending " << new_edges.size() << " new edges\n";);
            // edges.append(parlay::make_slice(new_edges));

            // // debug(std::cout << "srcs.size(): " << srcs.size() << std::endl;);
            srcs = parlay::filter(srcs, [&](size_t i)
                                  { size_t tmp = radius[cluster_index[i]]; if (tmp <= SSSP_CUTOFF) {
                                    radius[cluster_index[i]] = SIZE_MAX;
                                    return false;
                                    }
                                    return true; });
            // // debug(std::cout << "srcs.size(): " << srcs.size() << std::endl;);

            // // debug(std::cout << "radius.size(): " << radius.size() << std::endl;);
            auto filtered_radius = parlay::filter(radius, [&](size_t i)
                                                  { return i != SIZE_MAX; });
            // // debug(std::cout << "radius.size(): " << filtered_radius.size() << std::endl;);

            // debug(assert(notContained(filtered_radius, SIZE_MAX)););

            parallel_for(0, cluster_index.size(), [&](size_t i)
                         { cluster_index[i] = -1; });
            parallel_for(0, srcs.size(), [&](size_t i)
                         { cluster_index[srcs[i]] = i; });

            radius = filtered_radius;
            // debug(assert(notContained(radius, SIZE_MAX)););

            // debug(assert(srcs.size() == radius.size()););
            // debug(assert(srcs.empty() == radius.empty()););
            // debug(
            // if (srcs.size() < srcsSize_tmp) {
            //     assert(edges.size() >= edgesSize_tmp);
            // } else {
            //     assert(edges.size() == edgesSize_tmp);
            // });

            // debug(correctCluster(G, srcs, cluster););
            // debug(correctClusterIndex(cluster_index, srcs););
            // debug(noDoubleEdges(edges););

            // debug(std::cout << "exiting BFS_radius\n";);

            return;
        }

        template <class Graph>
        inline void BallGrowing(Graph &G,
                                sequence<uintE> &srcs,
                                sequence<uintE> &ball,
                                sequence<uintE> &ballShell,
                                sequence<uintE> &parents,
                                sequence<uintE> &cluster,
                                sequence<size_t> &cluster_index,
                                sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> &edges,
                                sequence<size_t> &radius,
                                size_t &max_rad,
                                parlay::random_generator &gen,
                                sequence<size_t> &dists,
                                sequence<size_t> &r0s)
        {
            using W = typename Graph::weight_type;
            // debug(std::cout << "entering BallGrowing\n";);
            // debug(assert(!srcs.empty()););
            // debug(assert(cluster.size() == G.n););
            // debug(assert(cluster_index.size() == G.n););
            // debug(correctCluster(G, srcs, cluster););
            // debug(correctClusterIndex(cluster_index, srcs););
            // debug(noDoubleEdges(edges););

            // sequence<size_t> dists;
            BFS_radius(G, srcs, radius, parents, cluster, cluster_index, dists, max_rad, edges);
            // std::cout << "returned from BFS_radius, srcs.size(): " << srcs.size() << " radius.size(): " << radius.size() << std::endl;
            if (srcs.empty())
            {
                // debug(noDoubleEdges(edges););
                // std::cout << "We are finished\n";
                return;
            }

            // debug(assert(parents.size() == G.n););
            // debug(assert(dists.size() == G.n););
            // debug(assert(srcs.size() == radius.size()););
            // debug(noDoubleEdges(edges););

            sequence<uintE> vertices = sequence<uintE>::from_function(G.n, [&](uintE i)
                                                                      { return i; });
            // std::cout << "genereted vertex array\n";

            // std::random_device rd;
            // std::mt19937 gen(rd());
            std::uniform_real_distribution<double> distr(0.5, 2 / 3.0);

            // sequence<size_t> r0s = sequence<size_t>::from_function(srcs.size(), [&](size_t i)
            //                                                        { auto r = gen[i]; return (int)floor(radius[i] * distr(r)); });
            r0s = sequence<size_t>::from_function(srcs.size(), [&](size_t i)
                                                  { auto r = gen[i]; return (int)floor(radius[i] * distr(r)); });

            // debug(assert(r0s.size() == radius.size()););

            // debug(assert(inRange(r0s, radius, 0.5, 2.0 / 3.0)););

            auto inBall = [&](const uintE &v)
            { return dists[v] <= r0s[cluster_index[cluster[v]]]; };
            auto inBallShell = [&](const uintE &v)
            { return (r0s[cluster_index[cluster[v]]] < dists[v] && dists[v] <= r0s[cluster_index[cluster[v]]] + 1); };
            auto cut = [&](const uintE &src, const uintE &ngh, const W &wgh)
            {
                if ((inBall(src) && inBallShell(ngh)) || (inBall(ngh) && inBallShell(src)))
                {
                    return 1;
                }
                else
                {
                    return 0;
                }
            };

            ball = parlay::filter(vertices, inBall);
            ballShell = parlay::filter(vertices, inBallShell);

            // Isolate Ball from rest of graph to make life easier later on. Might be inefficient
            // std::cout << "FilterEdges\n";

            // debug(for (uintE i = 0; i < G.n; i++) {
            //     if(inBall(i)) {
            //         assert(contains(ball, i));
            //     }
            //     if (inBallShell(i)) {
            //         assert(contains(ballShell, i));
            //     }
            // });

            // debug(correctCluster(G, srcs, cluster););
            // debug(correctClusterIndex(cluster_index, srcs););

            // // debug(std::cout << "cluster[0]: " << cluster[0] << ", cluster[4]: " << cluster[4] << std::endl;);
            // // debug(std::cout << "dists[0]: " << dists[0] << ", dists[4]: " << dists[4] << std::endl;);
            // // debug(std::cout << "r0 for cluster[0]: " << r0s[cluster_index[cluster[0]]] << std::endl;);
            // // debug(std::cout << "is inBallShell(4): " << inBallShell(4) << std::endl;);
            // // debug(std::cout << "contains(ballShell, 4): " << contains(ballShell, (uintE) 4) << std::endl;);
            // // debug(if (inBallShell(4)) {std::cout << "true\n\n";});

            // // debug(
            //     auto pars = BFS(G, 0);
            //     if (isReachable(pars, 4)) {
            //         std::cout << "Is Reachable\n";
            //     } else {
            //         std::cout << "Not reachable\n";
            //     });

            filterEdges(G, cut);

            // // debug(std::cout << "cluster[0]: " << cluster[0] << ", cluster[4]: " << cluster[4] << std::endl;);

            // // debug(
            //     pars = BFS(G, 0);
            //     if (isReachable(pars, 4)) {
            //         std::cout << "Is Reachable\n";
            //     } else {
            //         std::cout << "Not reachable\n";
            //     });

            // debug(isBallCut(G, ball, ballShell););

            // debug(for (int i = 0; i < cluster.size(); i++) {
            //     if (parents[i] != UINT_E_MAX)
            //     {
            //         assert(cluster[i] == cluster[parents[i]]);
            //     }
            // });
            // debug(assert(isDisjoint(G, srcs)););
            // debug(noDoubleEdges(edges););

            // debug(std::cout << "exiting BallGrowing\n";);

            return;
        }

    } // namespace LowAvgStrech
} // namespace gbbs