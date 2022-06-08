from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
import pandas as pd
import sys
from copy import deepcopy
from collections import Counter
from argparse import ArgumentParser


def ab2qpg(G, fname):
    known_nodes = set()
    s = []
    for i, (node, data) in enumerate(G.items()):


        if data[0] is not None and data[1] is None:
            if data[5] is None and data[6] is None: #samples
                s.append("label\t%s\t%s" % (node, node))
            s.append("edge\te_%s\t%s\t%s" % (i, data[0], node))

        if data[0] is not None and data[1] is not None: #admixture
            tpl = node, data[0], data[1], data[2], (1. - data[2])
            s.append("admix\t%s\t%s\t%s\t%s\t%s" % tpl)

        known_nodes.add(node)
        if data[0] is not None:
            known_nodes.add(data[0])
        if data[1] is not None:
            known_nodes.add(data[1])

    for node in G:
        known_nodes.remove(node)

    assert len(known_nodes) == 1
    node = known_nodes.pop()
    s.append("root\t%s" % node)

    s.sort(reverse=True)

    with open(fname, 'w') as f:
        f.write("\n".join(s))
    print('written qp graph to file ', fname)

def main(args):

    parser = ArgumentParser(usage='pipeline for plotting posterior distribution summaries.', version='1.0.0')

    parser.add_argument('--posterior_distribution_file', required=True, type=str, help='The file containing posterior distributions from the "AdmixtureBayes posterior" command. It needs the two columns "pops" and topology.')
    parser.add_argument('--no_topologies_to_plot', default=10, type=int, help='The number of the most posterior topologies to transform to qpgraphs')
    parser.add_argument('--consensus_threshold', default=[0.25, 0.5, 0.75, 0.9, 0.95, 0.99], type=float, nargs='+',
                        help='The posterior thresholds for which to draw different consensus trees.')
    parser.add_argument('--sep', default=',', type=str, help='the separator used in the input file')
    parser.add_argument('--outfile_prefix',default='', type=str, help='beginning of all files where the qp graphs are saved')

    options=parser.parse_args(args)
    df = pd.read_csv(options.posterior_distribution_file, sep=options.sep, usecols=['string_tree','topology'])
    stree_list = df['string_tree'].tolist()

    nodes=stree_list[0].split('=')[:-1]

    topologies=df['topology'].tolist()

    counter=Counter(topologies)

    rd=counter.most_common(options.no_topologies_to_plot)

    for n,(string_topology, common_ness) in enumerate(rd):
        index=topologies.index(string_topology)
        stree=stree_list[index]
        Rtree=identifier_to_tree_clean(stree.split('=')[-1], leaves=generate_predefined_list_string(deepcopy(nodes)))
        ab2qpg(Rtree, options.outfile_prefix+'qp'+str(n+1)+'.graph')
