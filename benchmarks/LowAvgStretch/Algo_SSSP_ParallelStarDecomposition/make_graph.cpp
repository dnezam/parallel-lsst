#include "gbbs/gbbs.h"
#include "gbbs/graph_io.h"

int main()
{
    const std::vector<gbbs::gbbs_io::Edge<gbbs::empty>> kEdges{
        {0, 3}, {0, 5}, {0, 6}, {0, 7}, //0,0
        {1, 2}, {1, 3}, //0,2
        {2, 3}, {2, 8}, //0,2
        {3, 5}, //0,1
        {4, 9}, {4, 10}, //10,1
        {5, 6},  //0,1
        {6, 7}, //0,1
        {7, 8}, //0,1
        {8, 9}, //0,2;10,3
        //9: 0,3;10,2
        {10, 11}, {10, 12}, //10,0
        {11, 13}, {11, 14}, //10,1
        {12, 15}, //10,1
        //13: 10,2
        //14: 10,2
        //15: 10,2
    };
    
    gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::empty> graph{gbbs::gbbs_io::edge_list_to_symmetric_graph(kEdges)};

    gbbs::gbbs_io::write_graph_to_file("./inputs/rMatGraph_Test_10_20", graph);
}