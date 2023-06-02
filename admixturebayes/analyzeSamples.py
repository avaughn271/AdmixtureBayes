from argparse import ArgumentParser
from copy import deepcopy
import sys
import pandas as pd

from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, admixture_sorted_unique_identifier, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _add_to_waiting, get_populations

from Rtree_operations import (get_leaf_keys, get_real_parents, get_real_children, rename_root, change_admixture, get_categories,
                              screen_and_prune_one_in_one_out, remove_non_mixing_admixtures, node_is_non_admixture)

def leave_node(key, node, population, target_nodes, follow_branch):
    if node_is_non_admixture(node):
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch=1, population=new_pop, target_nodes=target_nodes, child_key=key, dependent='none')]

class follow_branch_class(object):
    
    def __init__(self, sub_graph_nodes):
        self.sub_graph_nodes=sub_graph_nodes
        self.seen_merging=False
        
    def __call__(self, parent_key, branch, population, target_nodes, child_key, dependent='none'):
        if self.seen_merging:
            return parent_key, population, dependent
        subset=population.subset_of_the_candidates(self.sub_graph_nodes)
        if subset=='partly':
            target_nodes.append((child_key, branch))
        elif subset=='all':
            self.seen_merging=True
        return parent_key, population, dependent

def get_branches_to_keep(tree, subgraph_keys):
    node_keys=get_leaf_keys(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    follow_branch=follow_branch_class(subgraph_keys)
    ready_nodes=list(zip(node_keys,pops))
    waiting_nodes={}
    taken_nodes=[]
    target_nodes=[]
    while True:
        for key,pop in ready_nodes:
        
            upds=leave_node(key, tree[key], pop, target_nodes, follow_branch)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return target_nodes

def find_root_name(tree):
    parents_seen=set()
    for k in tree:
        ps=get_real_parents(tree[k])
        for p in ps:
            parents_seen.add(p)
    rootset=parents_seen-set(tree.keys())
    return next(iter(rootset))

def prune_to_subtree(tree,branches_to_keep):
    sub_tree={}
    for key,b in branches_to_keep:
        sub_tree[key]=tree[key]
    root_name=find_root_name(sub_tree)
    sub_tree=rename_root(sub_tree, root_name)
    sub_tree=remove_empty_children(sub_tree)
    sub_tree=screen_and_prune_one_in_one_out(sub_tree)
    return sub_tree

def get_subtree(tree, subgraph_keys):
    tree2=remove_non_mixing_admixtures(deepcopy(tree))
    branches_to_keep=get_branches_to_keep(tree2, subgraph_keys)
    return prune_to_subtree(tree2, branches_to_keep)

def remove_empty_children(tree):
    for k in tree:
        child_keys=get_real_children(tree[k])
        children_to_keep=[]
        for child_key in child_keys:
            if child_key in tree:
                children_to_keep.append(child_key)
        if len(child_keys)!=len(children_to_keep):
            for n,ch in enumerate(children_to_keep):
                tree[k][5+n]=ch
            for n in range(n+1,2):
                tree[k][5+n]=None
    return tree

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             pre_thin_data_set_function=identity, 
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols, dtype={'no_admixes':object})
    df = df[cols]
    df= pre_thin_data_set_function(df)
    all_results=[]
    
    for n,(i,r) in enumerate(df.iterrows()):
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)

        for row_summarize_function in row_summarize_functions:
            add_dic, skip=row_summarize_function(**d_dic)
            d_dic.update(add_dic)
        all_results.append(thinned_d_dic(d_dic))
    return all_results

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, subnodes=[], outgroup_name=''):
        self.nodes=sorted(nodes_to_be_sorted)
        self.subnodes=subnodes
        self.outgroup_name=outgroup_name
        
    def __call__(self, tree, **not_needed):
        first_level=tree.split('-')[0]
        no_pops=len(first_level.split('.'))
        if no_pops==len(self.nodes): #checking if
            Rtree=identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        elif no_pops==len(self.nodes)+1 and self.outgroup_name:
            self.nodes=sorted(self.nodes+[self.outgroup_name])
            Rtree = identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        else:
            assert False, 'Either the outgroup name was not specified or something is seriously wrong because' \
                          ' the number of nodes did not match the size of the trees'

        if self.subnodes:#DETTE TAGER IKKE ORDENTLIG HOJDE FOR KOVARIANSMATRICERNE SOM BLIVER FORKERTE
            try:
                Rtree=get_subtree(Rtree, self.subnodes)
            except AssertionError:
                print('input_tree', tree)
                print('nodes', self.nodes)
                print('subnodes', self.subnodes)
                assert False
        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, outgroup_name='out', subnodes=[]):
        self.outgroup_name=outgroup_name
        self.subnodes=subnodes
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        full_tree=deepcopy(Rtree)
        if self.subnodes:
            full_tree=get_subtree(full_tree, self.subnodes)
        return {'full_tree':full_tree}, False

def get_subpops(pops, sub_graph_keys):
    ss_subgraph_keys=set(sub_graph_keys)
    new_pops=[]
    for pop in pops:
        new_pop=ss_subgraph_keys.intersection(pop.split('.'))
        new_pops.append('.'.join(new_pop))
    if '' in new_pops:
        new_pops.remove('')
    return '_'.join(sorted(list(set(new_pops))))

class topology(object):
    
    def __init__(self, nodes):
        self.nodes=nodes
        
    def __call__(self, Rtree=None, **kwargs):
        if 'string_tree' in kwargs:
            topology=kwargs['string_tree'].split('=')[-1].split(';')[0]
            return {'topology':topology},False
        if Rtree is None:
            Rtree=kwargs['full_tree']
        top=admixture_sorted_unique_identifier(Rtree, leaf_order=self.nodes, not_opposite=True)
        return {'topology':top}, False
    
class get_pops(object):
    
    def __init__(self, keys_to_include=None):
        self.keys_to_include=keys_to_include
            
    def __call__(self, full_tree=None, **kwargs):
        if full_tree is None:
            tree=kwargs['Rtree']
        else:
            tree=full_tree
        pops=get_populations(tree, keys_to_include=self.keys_to_include)
        return {'pops':'-'.join(pops)}, False

class thinning(object):
    
    def __init__(self, burn_in_fraction=None, total=None, **values_to_filter_by):
        self.burn_in_fraction=burn_in_fraction
        self.total=total
        self.values_to_filter_by=values_to_filter_by
        
    def __call__(self, df):
        #first_removing burn-in
        n=len(df)
        print('Dataframe read with ' + str(n) + ' samples.')
        if self.burn_in_fraction is not None:
            df=df[int(n*self.burn_in_fraction):]
        print('Burn-in of ' + str(n-len(df)) + ' samples removed. There are now ' + str(len(df)) + ' samples.')
        for column,value in list(self.values_to_filter_by.items()):
            print('filtering on ', column,'==',value)
            df=df.loc[df[column]==value,:]
        if self.total is not None:
            n=len(df)
            stepsize = self.total
            df=df[::stepsize]
            print('Thinning every ' + str(stepsize) + ' samples complete. There are now ' + str(len(df)) + ' samples')
        return df

def run_posterior_main(args):

    parser = ArgumentParser(usage='pipeline for post analysis')

    parser.add_argument('--mcmc_results', required=True, type=str, help='The output file from an AdmixtureBayes run.')
    parser.add_argument('--subnodes', default=[], type=str, nargs='+',
                        help='The subset of populations to perform the analysis on. If not declared, the analysis will be done on the full dataset.')
    parser.add_argument('--result_file', default='thinned_samples.csv', type=str,
                        help='The resulting file. It will be comma-separated and contain one column per summary plus a header.')

    parser.add_argument('--thinning_rate', default=10, type=int,
                        help='thinning rate')
    parser.add_argument('--burn_in_fraction', default=0.5, type=float,
                        help='the proportion of the rows that are discarded as burn in period')
    parser.add_argument('--save_summaries', default=['no_admixes', 'topology', 'pops','string_tree'], nargs='*', type=str,
                        help='The list of summaries to save')

    options= parser.parse_args(args)

    if not options.subnodes:
        df = pd.read_csv(options.mcmc_results)

        n=len(df)
        print('Dataframe read with ' + str(n) + ' samples.')
        if options.burn_in_fraction is not None:
            df=df[int(n*options.burn_in_fraction):]
        print('Burn-in of ' + str(n-len(df)) + ' samples removed. There are now ' + str(len(df)) + ' samples.')
        stepsize = options.thinning_rate
        df=df[::stepsize]
        print('Thinning every ' + str(stepsize) + ' samples complete. There are now ' + str(len(df)) + ' samples')

        treee = df['tree'].values.tolist()
        no_admixes = df['no_admixes'].values.tolist()
        descendantsets = df['descendant_sets'].values.tolist()
        pops = descendantsets[0]
        pops = pops.replace("-", " ")
        pops = pops.replace(".", " ")
        pops = pops.split(" ")
        uniquepops = []
        for x in pops:
            if x not in uniquepops:
                uniquepops.append(x)
        uniquepops.sort()

        uniquepopstring = uniquepops[0]
        for j in range(len(uniquepops) - 1):
            uniquepopstring = uniquepopstring +"=" + uniquepops[j + 1]
        topologyy = ["temp"] * len(treee)

        string_tree = ["temp"] * len(treee)
        for i in range(len(topologyy)):
            topologyy[i] = treee[i].split(";")[0]
            string_tree[i] = uniquepopstring + "=" + treee[i]
        d = {'topology': topologyy, 'pops': descendantsets, 'no_admixes':no_admixes, 'string_tree': string_tree}
        df = pd.DataFrame(data=d)
        df.to_csv(options.result_file, index=False)
    else:
        if options.subnodes:
            subnodes_wo_outgroup=options.subnodes
            subnodes_with_outgroup=options.subnodes
        else:
            subnodes_with_outgroup=[]
            subnodes_wo_outgroup=[]

        thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.thinning_rate)

        totallist = []
        a = pd.read_csv(options.mcmc_results, nrows=3)
        stringg = (a.loc[0,["descendant_sets"]])[0]
        stringg = (stringg.split('-'))
        for i in stringg:
            totallist.extend(i.split('.'))
        nodes = []
        for i in totallist:
            if i not in nodes:
                nodes.append(i)
        nodes.sort()

        nodes_wo_outgroup = deepcopy(nodes)
        nodes_with_outgroup = deepcopy(nodes_wo_outgroup)
        nodes_with_outgroup.sort()
        nodes_wo_outgroup.sort()
        subnodes_with_outgroup.sort()
        subnodes_wo_outgroup.sort()

        row_sums=[]

        class pointers(object):

            def __init__(self):
                self.count=0
                self.dic={}

            def __call__(self, name):
                self.dic[name]=self.count
                self.count+=1

            def __getitem__(self, key):
                return self.dic[key]

        name_to_rowsum_index=pointers()

        row_sums.append(make_Rtree(deepcopy(nodes_wo_outgroup), subnodes=subnodes_wo_outgroup, outgroup_name=''))
        name_to_rowsum_index('Rtree')
        row_sums.append(make_full_tree(outgroup_name='', subnodes=options.subnodes))
        name_to_rowsum_index('full_tree')

        if options.subnodes:
            nodes=options.subnodes
            nodes_with_outgroup=subnodes_with_outgroup
            nodes_wo_outgroup=subnodes_wo_outgroup
        else:
            nodes=nodes_with_outgroup
        options.save_summaries.remove('string_tree')
        row_sums.append(topology(nodes=nodes))
        name_to_rowsum_index('topology')
        row_sums.append(get_pops(keys_to_include=nodes))
        name_to_rowsum_index('pops')

        def save_thin_columns(d_dic):
            return {summ:d_dic[summ] for summ in list(set(options.save_summaries+[]))}
        all_results=iterate_over_output_file(options.mcmc_results,
                                                cols=['tree', 'add', 'layer', 'no_admixes'],
                                                pre_thin_data_set_function=thinner,
                                                row_summarize_functions=row_sums,
                                                thinned_d_dic=save_thin_columns)

        if True:
            summaries=list(all_results[0].keys())
            with open(options.result_file, 'w') as f:
                f.write(','.join(summaries)+'\n')
                for row in all_results:
                    s_summs=[str(row[summ]) for summ in summaries]
                    f.write(','.join(s_summs)+ '\n')
    
if __name__=='__main__':
    import sys
    run_posterior_main(sys.argv[1:])