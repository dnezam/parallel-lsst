import networkit as nk
import snap


def snap_to_nx(input_fn):
    g_snap = snap.LoadEdgeList(snap.TUNGraph, input_fn)
    g_nx = nk.Graph()
    max_id = 0
    for n in g_snap.Nodes():
        n = n.GetId()
        if n > max_id:
            max_id = n
    g_nx.addNodes(max_id)
    for e in g_snap.Edges():
        try:
            g_nx.addEdge(e.GetSrcNId(), e.GetDstNId())
        except RuntimeError:
            print('skr')
    g_nx = nk.graphtools.toUndirected(g_nx)
    cc = nk.components.ConnectedComponents(g_nx)
    g_connect = cc.extractLargestConnectedComponent(g_nx, True)
    # nk.graphtools.getCompactedGraph(g_connect, [0]*g_connect.numberOfNodes())
    nk.writeGraph(g_connect, 'tools/cite-_{}_{}_BIN'.format(g_connect.numberOfNodes(), g_connect.numberOfEdges()), nk.Format.GraphToolBinary)

snap_to_nx('tools/cit-Patents.txt')