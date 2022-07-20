from node_structure import Node
from collections import Counter
import pandas as pd
from argparse import ArgumentParser
from tree_statistics import generate_predefined_list_string,topological_identifier_to_tree_clean, identifier_to_tree
from copy import deepcopy
from Rtree_operations import node_is_admixture, rename_key, get_admixture_proportion_from_key, get_all_admixture_origins
import sys

def get_numeric(string_tree):
    branch_lengths_string, admixture_proportion_string= string_tree.split(';')[1:]
    branch_lengths=list(map(float,branch_lengths_string.split('-')))
    if admixture_proportion_string:
        admixture_proportions=list(map(float, admixture_proportion_string.split('-')))
    else:
        admixture_proportions=[]
    return branch_lengths, admixture_proportions

def branch_and_proportion_quantiles(list_of_string_trees):
    ''' extracts the branch lengths and admixture proportions and returns them as tuples of four. the list should not be empty'''
    branches=[]
    admixtures=[]
    for string_tree in list_of_string_trees:
        b,a=get_numeric(string_tree)
        branches.append(b)
        admixtures.append(a)
    bdat=pd.DataFrame.from_records(branches)
    adat=pd.DataFrame.from_records(admixtures)
    bresults=[]
    aresults=[]
    for n,(lower,mean, upper) in enumerate(zip(bdat.quantile(0.025),bdat.mean(),bdat.quantile(0.975))):
        bresults.append(('c'+str(n+1), lower,mean,upper))
    if len(admixtures[0])>0:
        for n,(lower,mean, upper) in enumerate(zip(adat.quantile(0.025),adat.mean(),adat.quantile(0.975))):
            aresults.append(('ax'+str(n+1), lower,mean,upper))
    return bresults, aresults


from tree_to_data import unzip
    
def read_one_line(filename):
    if filename.endswith('.gz'):
        filename=unzip(filename)
    with open(filename, 'r') as f:
        return f.readline().rstrip().split()

def main(args):
    parser = ArgumentParser(usage='pipeline for plotting posterior distribution summaries.')

    parser.add_argument('--posterior', required=True, type=str, help='The file containing posterior distributions from the "AdmixtureBayes posterior" command. It needs the two columns "pops" and topology.')
    parser.add_argument('--plot', choices=['consensus_trees', 'top_minimal_topologies', 'top_trees','estimates'], required=True,
                        help='The type of plot to make. Choose between: 1) consensus_trees. '
                             'It plots an admixture graph based on all nodes that have a higher (marginal) posterior probability of X. '
                             'Different X\'s can be supplied with the command --consensus_threshold \n'
                             '2) top_minimal_topologies. It plots the X highest posterior combinations of node types '
                             'and creates the corresponding minimal topologies.  X can be supplied through the command --top_minimal_topologies_to_plot'
                             '3) top_trees. It plots the X highest posterior topologies. X can be supplied by the command --top_trees_to_plot.'
                             '4) estimates. It creates a table  with continuous parameters estimated from the posterior sample'
                             'It also plots the concerned topologies with labels. It does this for either the X highest posterior topologies '
                             'or the topologies specified by --estimate_topologies.'
                             'If no --estimate_topologies is set, X can be set by top_trees_to_estimate (default=3). ')
    parser.add_argument('--outgroup', default='outgroup', help='name of the outgroup to plot')
    parser.add_argument('--prefix', default='', type=str, help='string to prepend before each file created by this routine. '
                                                               'That means that any rankings written to a file by setting --write_rankings or --write_estimates_to_file will have this prefix prepended.')
    parser.add_argument('--consensus_thresholds', default=[0.25, 0.5, 0.75, 0.9, 0.95, 0.99], type=float, nargs='+',
                        help='The posterior thresholds for which to draw different consensus trees.')
    parser.add_argument('--top_minimal_topologies_to_plot', type=int, default=3,
                        help='The number of node trees (or minimal topologies) to plot')
    parser.add_argument('--top_trees_to_plot', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--top_trees_to_estimate', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--estimate_topologies', default=[], type=str, nargs='+',
                        help='The topologies whose conitnuous parameters should be estimated with "--plot estimates"')
    parser.add_argument('--write_estimates_to_file', default=[], type=str, nargs='+',
                        help='The file in which to put the tables when plotting estimates. ')
    parser.add_argument('--write_rankings', type=str, default='', help='if a file is supplied here, the natural rankings for each of the plots is written here.')
    parser.add_argument('--rankings_to_write_to_file', type=int, default=1000,
                        help='the number of rankings(nodes, min topology or topology depending on --plot) to write to the ranking file.')
    parser.add_argument('--dont_annotate_node_posterior', default=False, action='store_true',
                        help='This will not color the nodes according to their posterior probability.')
    parser.add_argument('--popup', default=False, action='store_true')
    parser.add_argument('--sep', default=',', type=str, help='the separator used in the input file')

    options= parser.parse_args(args)

    def combine_nodes(node_structure, new_node, seen_sets):
        candidate=new_node.name
        seen=[]
        for lists_of_fixed_size in seen_sets[::-1]:
            for attached_branch in lists_of_fixed_size:
                if( attached_branch.issubset(candidate) and
                   ((not attached_branch.issubset(seen)) or (not node_structure[attached_branch].has_parent()))):
                    seen.extend(list(attached_branch))
                    new_node.add_child(node_structure[attached_branch])
                    node_structure[attached_branch].add_parent(new_node)
        return node_structure

    def node_combinations_to_node_structure(node_combinations):
        length_sorted={}
        for node_combination in node_combinations:
            leaves=frozenset(node_combination.split('.'))
            k=len(leaves)
            if k in length_sorted:
                length_sorted[k].append(leaves)
            else:
                length_sorted[k]=[leaves]
        length_sorted_list=[length_sorted.get(k,[]) for k in range(1,max(length_sorted.keys())+1)]
        #length_sorted_list is of the form [[[A],[B],[C]],[[A,B],[B,C]],...,[[A,B,C]]]
        node_structure={}
        for leaf_node in length_sorted_list[0]:
            node_structure[leaf_node]=Node(leaf_node)
        added_sets=[length_sorted_list[0]]
        for lists_of_fixed_size in length_sorted_list[1:]:
            for branch_set in lists_of_fixed_size:
                new_node=Node(branch_set)
                combine_nodes(node_structure, new_node, added_sets)
                node_structure[branch_set]=new_node
            added_sets.append(lists_of_fixed_size)
        return node_structure

    if options.plot=='consensus_trees' or options.plot=='top_minimal_topologies':
        df = pd.read_csv(options.posterior, sep=options.sep, usecols=['pops'])
        nodes_list = df['pops'].tolist()
        seen_combinations = {}
        for nodes in nodes_list:
            for node in nodes.split('-'):
                seen_combinations[node] = seen_combinations.get(node, 0) + 1
        N = len(nodes_list)
        if options.plot=='consensus_trees':
            node_combinations = []
            for threshold in options.consensus_thresholds:
                total_threshold = int(N * threshold)
                final_node_combinations = [k for k, v in list(seen_combinations.items()) if v > total_threshold]
                node_combinations.append(final_node_combinations)
            if not options.dont_annotate_node_posterior:
                node_count_dic={frozenset(k.split('.')):float(v)/N for k,v in list(seen_combinations.items())}
            else:
                node_count_dic=None
            for i, final_node_combinations in enumerate(node_combinations):
                final_node_structure = node_combinations_to_node_structure(final_node_combinations)
                from tree_plotting import plot_node_structure_as_directed_graph
                plot_node_structure_as_directed_graph(final_node_structure, drawing_name=options.prefix+'consensus_'+str(int(100*options.consensus_thresholds[i]))+'.png', node_dic=node_count_dic,  popup=options.popup)
            if options.write_rankings:
                with open(options.prefix+options.write_rankings, 'w') as f:
                    c = Counter(seen_combinations)
                    to_write = c.most_common(options.rankings_to_write_to_file)
                    for node, frequency in to_write:
                        f.write(node+','+str(float(frequency)/N)+'\n')
        elif options.plot=='top_minimal_topologies':
            c=Counter(nodes_list)
            to_plots=c.most_common(options.top_minimal_topologies_to_plot)
            if options.write_rankings:
                with open(options.prefix+options.write_rankings, 'w') as f:
                    for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                        f.write(tree + ',' + str(float(frequency) / N) + '\n')
            if not options.dont_annotate_node_posterior:
                c=Counter(seen_combinations)
                node_count_dic={frozenset(key.split('.')):float(count)/N for key,count in c.most_common(1000)}
            else:
                node_count_dic=None
            from tree_plotting import plot_node_structure_as_directed_graph
            for i, (to_plot,count) in enumerate(to_plots):
                node_structure = node_combinations_to_node_structure(to_plot.split('-'))
                plot_node_structure_as_directed_graph(node_structure, drawing_name=options.prefix+'minimal_topology_' +str(i+1)+'.png',
                                                          node_dic=node_count_dic,  popup=options.popup)
    elif options.plot=='top_trees':
        df = pd.read_csv(options.posterior, sep=options.sep, usecols=['pops','topology'])
        trees_list = df['topology'].tolist()
        no_leaves=len(trees_list[0].split('-')[0].split('.'))
        N=len(trees_list)
        c = Counter(trees_list)
        to_plots = c.most_common(options.top_trees_to_plot)

        if True:
            nodes=df['pops'].tolist()[0].split('-')
            leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
            if len(leaves)==no_leaves:
                pass #everything is good
            elif len(leaves)==no_leaves-1:
                leaves.append(options.outgroup)
            else:
                assert False, 'The number of leaves could not be obtained'
            leaves=sorted(leaves)

        if options.write_rankings:
            with open(options.prefix+options.write_rankings, 'w') as f:
                for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                    f.write(tree + ',' + str(float(frequency) / N) + '\n')

        from tree_plotting import plot_as_directed_graph
        for i, (to_plot, count) in enumerate(to_plots):
            tree=topological_identifier_to_tree_clean(to_plot, leaves=generate_predefined_list_string(deepcopy(leaves)))
            plot_as_directed_graph(tree,drawing_name=options.prefix+'topology_' + str(i + 1) + '.png', popup=options.popup)
    elif options.plot=='estimates':
        try:
            df = pd.read_csv(options.posterior, sep=options.sep, usecols=['string_tree', 'topology', 'pops'])
        except ValueError as e:
            raise Exception('Unexpected columns in the posterior_distribution file. Did you turn on the --faster flag in AdmixtureBayes posterior?')

        topologies_list = df['topology'].tolist()
        string_tree_list=df['string_tree'].tolist()
        if options.estimate_topologies:
            cleaned_topology_list=[s.split('=')[-1].split(';')[0] for s in options.estimate_topologies]
        else:
            c = Counter(topologies_list)
            to_plots = c.most_common(options.top_trees_to_estimate)
            cleaned_topology_list=[d[0] for d in to_plots]
        no_leaves = len(topologies_list[0].split('-')[0].split('.'))

        if True:
            nodes=df['pops'].tolist()[0].split('-')
            leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
            if len(leaves)==no_leaves:
                pass #everything is good
            elif len(leaves)==no_leaves-1:
                leaves.append(options.outgroup)
            else:
                assert False, 'The number of leaves could not be obtained'
            leaves=sorted(leaves)

        from tree_plotting import plot_as_directed_graph
        for i, to_plot in enumerate(cleaned_topology_list):
            relevant_string_trees=[]
            for string_tree, topology in zip(string_tree_list, topologies_list):
                if topology==to_plot:
                    relevant_string_trees.append(string_tree)
            branches_intervals, admixture_proportion_intervals=branch_and_proportion_quantiles(relevant_string_trees)
            branch_names=[branches_interval[0] for branches_interval in branches_intervals]

            admixture_names=[ad[0] for ad in admixture_proportion_intervals]
            tree = identifier_to_tree(to_plot,
                                      leaves=generate_predefined_list_string(deepcopy(leaves)),
                                      branch_lengths=generate_predefined_list_string(deepcopy(branch_names)),
                                      admixture_proportions=generate_predefined_list_string(deepcopy(admixture_names)))

            org_keys = list(tree.keys())
            for key in org_keys:
                node = tree[key]
                if node_is_admixture(node):
                    new_name = get_admixture_proportion_from_key(tree, key)
                    tree = rename_key(tree, key, new_name)
            adms=get_all_admixture_origins(tree)
            adm_interpretation={}
            for key, (branch_name, node_destination) in list(adms.items()):
                adm_interpretation[key]='For the lineages that pass through {}, this is the proportion that follows branch {} to node {}'.format(key, branch_name,node_destination)
            plot_as_directed_graph(tree, drawing_name=options.prefix+'topology_labels_' + str(i + 1) + '.png', plot_edge_lengths=True,  popup=options.popup)
            if options.write_estimates_to_file:
                branch_file=options.write_estimates_to_file[i*2+0]
                admixtures_file=options.write_estimates_to_file[i*2+1]
            else:
                branch_file=options.prefix+'branch_estimates_'+str(i+1)+'.txt'
                admixtures_file=options.prefix+'admixture_estimates_'+str(i+1)+'.txt'
            with open(branch_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%','mean','upper 95%'])+'\n')
                for v in branches_intervals:
                    f.write(','.join(map(str,v))+'\n')
            with open(admixtures_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%', 'mean', 'upper 95%','interpretation'])+'\n')
                for v in admixture_proportion_intervals:
                    f.write(','.join(map(str,list(v)+[adm_interpretation[v[0]]]))+'\n')

if __name__=='__main__':
    main(sys.argv[1:])
