import networkit as nk
import networkx as nx
import pathlib

base_dir = pathlib.Path('../test_inputs')

def _get_connected(g):
    if type(g) == nx.Graph:
        g = nk.nxadapter.nx2nk(g)
    cc = nk.components.ConnectedComponents(g)
    return cc.extractLargestConnectedComponent(g, True)

def gen_and_store_graph(algo_name, number_vertices, parameter):
    algo_base = base_dir.joinpath(algo_name)

    g = None
    if algo_name == 'ErdosRenyi':
        g = _get_connected(nk.generators.ErdosRenyiGenerator(number_vertices, parameter).generate())
    elif algo_name == 'BarabasiAlbert':
        g = _get_connected(nk.generators.BarabasiAlbertGenerator(parameter, number_vertices, n0=1).generate())
    else:
        raise ValueError("Invalid algo_name")

    n_connected = g.numberOfNodes()
    m_connected = g.numberOfEdges()

    algo_base.mkdir(parents=True, exist_ok=True)

    path = str(algo_base.joinpath('{}_{}_{}_{}_{}_BIN'.format(algo_name, number_vertices, parameter, n_connected, m_connected)))
    nk.writeGraph(g, path, nk.Format.GraphToolBinary)


ls = [8, 10, 12, 14, 15]  # n = 2^l

for l in ls:
    n = pow(2, l)

    gen_and_store_graph('ErdosRenyi', n, 0.01)
    gen_and_store_graph('ErdosRenyi', n, 10/n)
    gen_and_store_graph('BarabasiAlbert', n, 4)


    

