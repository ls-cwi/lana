import itertools
import random

n_graphs = 5
n_nodes = 5000
p_edge_graph = 0.0008
p_edge_compatibility = 0.002
data_folder = "data/gen/"
output_folder = "output/"


def generate_random_graph(num_nodes, p_edge, name):
    print "Generating graph '%s'." % name
    nodes = range(0,num_nodes)
    edges = filter(lambda x: x[0] < x[1] and random.random() < p_edge, itertools.product(nodes, nodes))
    return (nodes, edges, name)

def generate_compatibility_graph(g1, g2, p_edge):
    print "Generating compatibility graph '%s_%s'." % (g1[2], g2[2])
    g1_labels = [g1[2] + str(x) for x in g1[0]]
    g2_labels = [g2[2] + str(x) for x in g2[0]]
    return dict([(x, filter(lambda y: random.random() < p_edge, g2_labels)) for x in g1_labels])

def save_graph(g, file_name):
    print "Saving graph '%s'." % file_name

    file = open(file_name, "w")

    file.write( "graph [\n"
                "   comment \"%s graph\"\n"
                "   directed 0\n"
                "\n" % g[2])

    for node in g[0]:
        file.write( "   node [\n"
                    "      id %d\n"
                    "       label \"%s%d\"\n"
                    "   ]\n" % (node, g[2], node))

    for edge in g[1]:
        file.write( "   edge [\n"
                    "       source %d\n"
                    "       target %d\n"
                    "       label \"%s%d_%s%d\"\n"
                    "   ]\n" % (edge[0], edge[1], g[2], edge[0], g[2], edge[1]))

    file.write("] \n")

    file.close()

def save_compatibility_graph(g, file_name):
    print "Saving compatibility graph '%s'." % file_name
    file = open(file_name, "w")

    for (v1, vs2) in g.items():
        if vs2:
            file.write("%s\t%s\n" % (v1, ' '.join([x for x in vs2])))


def generate_run_script(graphs, data_folder, output_folder, file_name):
    print "Generating run script (%s)." % file_name
    file = open(file_name, 'w')

    for (g1, g2) in itertools.combinations(graphs, 2):
        file.write("build/lana -if1 0 -if2 0 -ifm 0 -g1 %s%s.gml -g2 %s%s.gml -gm %s%s_%s.seqSim -freq %s%s_%s-freq.csv $1; \n"
                % (data_folder, g1[2], data_folder, g2[2], data_folder, g1[2], g2[2], output_folder, g1[2], g2[2]))





def generate_data(n_graphs, n_nodes, p_edge_graph, p_edge_compatibility, data_folder, output_folder):
    if n_graphs > 26:
        print "Can't make more than 26 graphs!"

    graphs = [generate_random_graph(n_nodes, p_edge_graph, str(chr(ord('a')+x))) for x in range(0, n_graphs)]

    for g in graphs:
        save_graph(g, "%s%s.gml" % (data_folder, g[2]))

    for (g1, g2) in itertools.combinations(graphs, 2):
        gc = generate_compatibility_graph(g1, g2, p_edge_compatibility)
        save_compatibility_graph(gc, "%s%s_%s.seqSim" % (data_folder, g1[2], g2[2]))

    generate_run_script(graphs, data_folder, output_folder, "run_all.sh")


generate_data(n_graphs, n_nodes, p_edge_graph, p_edge_compatibility, data_folder, output_folder)
