// ugly include, just to make it work
#include "networkit/include/networkit/graph/Graph.hpp"
#include "networkit/include/networkit/io/GraphToolBinaryWriter.hpp"
#include <iostream>



#include <iostream>
#include <vector>
#include <cmath>


// Create a tree such that if BFS starts at the root, it is bad
std::vector<std::pair<int,int>> Create_Bad_Tree(int n, double p, int root, int &heighest_idx){
    std::vector<std::pair<int,int>> tree;
    int height = (int)ceil(n*p / 4);
    int left_branch_idx = root;
    int right_branch_idx = root;
    for(int i = heighest_idx; i < heighest_idx + 2*height; i += 2){
        // Add new edges to tree
        tree.push_back(std::make_pair(left_branch_idx,i));
        tree.push_back(std::make_pair(right_branch_idx,i+1));
        // update parent
        left_branch_idx = i;
        right_branch_idx = i+1;
    }

    heighest_idx += 2 * height;

    // Add edges still in the tree
    for(int i = heighest_idx; i < heighest_idx + 2*height; i+= 2){
        tree.push_back(std::make_pair(left_branch_idx,i));
        tree.push_back(std::make_pair(right_branch_idx,i+1));
    }

    // add edges reaching to the other side
    for(int i = heighest_idx; i < heighest_idx + 2*height; i+= 2){
        for(int j = heighest_idx; j < heighest_idx + 2*height; j+= 2){
        tree.push_back(std::make_pair(i,j+1));
    }
    }
    heighest_idx += 2 * height;

    return tree;
}

// Creates a graph that has bad average stretch for random BFS with probability 1-p
std::vector<std::pair<int,int>> Create_Bad_Graph(int n, double p){
    std::vector<std::pair<int,int>> res;

    int list_size = (int) ceil(1/p);

    std::cout << "nr trees: 1/p = " << list_size << "\n";
    std::cout << "height of tree: np/4 = " << (int)ceil(n*p / 4) << "\n";

    // create list of size 1/p
    for(int i = 0; i < list_size-1; i++){
        res.push_back(std::make_pair(i,i+1));
    }
    // for each vertex in the list create tree of size np
    int idx = list_size;
    for(int i = 0; i < list_size-1; i++){
        auto tree = Create_Bad_Tree(n,p, i, idx);
        res.insert(res.end(),tree.begin(),tree.end());
    }



    return res;
}


// Disclaimer, don't run with p = 1

int main() {
    int n = 1000000;
    double p = 0.001;
    auto res = Create_Bad_Graph(n,p);
    NetworKit::Graph g{n};
    for(int i = 0; i < res.size(); i++){
        g.addEdge(res[i].first, res[i].second);
    }
    auto writer = NetworKit::GraphToolBinaryWriter{};
    writer.write(g, "BadGraph_"+std::to_string(g.numberOfNodes())+"_"+std::to_string(g.numberOfEdges())+"_"+"BIN");


    return 0;
}