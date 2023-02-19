from copy import deepcopy
import pandas

def rename_root(tree, old_name):
    for _,node in list(tree.items()):
        if node[0]==old_name:
            node[0]='r'
        if (node[1] is not None and node[1]==old_name):
            node[1]='r'
    return tree

def update_parent_and_branch_length(tree, child_key, child_branch, new_parent, new_branch_length):
    tree[child_key][child_branch]=new_parent
    tree[child_key][child_branch+3]=new_branch_length
    return tree

def insert_children_in_tree(tree):
    children={key:[] for key in tree}
    for key in tree:
        parents = get_real_parents(tree[key])
        for parent in parents:
            if parent!='r':
                children[parent].append(key)
    for key in tree:
        tree[key]=_update_parents(tree[key], children[key])
    return tree

def _update_parents(node, new_parents):
    if len(new_parents)==1:
        res=node[:5]+[new_parents[0],None]
        return res
    if len(new_parents)==2:
        res=node[:5]+new_parents
        return res
    if len(new_parents)==0:
        res=node[:5]+[None]*2
        return res

def get_real_parents(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

class generate_numbered_nodes(object):

    def __init__(self, prefix='n', admixture_prefix='a'):
        self.node_count=0
        self.admixture_count=0
        self.prefix='n'
        self.admixture_prefix='a'

    def __call__(self, admixture=False):
        if admixture:
            self.admixture_count+=1
            return self.admixture_prefix+str(self.admixture_count)
        self.node_count+=1
        return self.prefix+str(self.node_count)

class generate_predefined_list_string(object):

    def __init__(self, listi):
        self.listi=listi

    def __call__(self):
        return self.listi.pop(0)

def identifier_to_tree(identifier, leaves=None, inner_nodes=None, branch_lengths=None, admixture_proportions=None):
    '''
    Transforms an identifier of the form qwert-uio-asdfg-jk into a dictionary tree using the generators of leaves, inner_nodes, branch_lengths and admixture_proportions.
    '''
    levels=identifier.split('-')
    n_leaves=len(levels[0].split('.'))
    if leaves is None:
        leaf_values=sorted(['s'+str(n+1) for n in range(n_leaves)])
    else:
        leaf_values=[leaves() for _ in range(n_leaves)]
    tree={leaf:[None]*5 for leaf in leaf_values}
    trace_lineages=[(leaf,0) for leaf in leaf_values]
    if inner_nodes is None:
        inner_nodes=generate_numbered_nodes('n')
    if branch_lengths is None:
        def f():
            return 1.0
        branch_lengths= f
    if admixture_proportions is None:
        def g():
            return 0.4
        admixture_proportions=g
    for level in levels:
        identifier_lineages=level.split('.')
        parent_index={}
        indexes_to_be_removed=[]
        for n,identifier_lineage in enumerate(identifier_lineages):
            if identifier_lineage=='c':
                ##there is a coalecence for the n'th lineage, and it should be replaced by a new lineage
                new_key=inner_nodes()
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                tree[new_key]=[None]*5
                parent_index[n]=new_key
                trace_lineages[n]=(new_key,0)
            elif identifier_lineage=='w':
                pass
            elif identifier_lineage=='a':
                new_key=inner_nodes(admixture=True)
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                new_admixture_proportion=admixture_proportions()
                tree[new_key]=[None,None,new_admixture_proportion,None,None]
                trace_lineages[n]=(new_key,0)
                trace_lineages.append((new_key,1))
            else:
                try:
                    new_key=parent_index[int(identifier_lineage)]
                except KeyError as e:
                    print(e)

                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                indexes_to_be_removed.append(n)

        ##remove lineages
        trace_lineages=[trace_lineage for n,trace_lineage in enumerate(trace_lineages) if n not in indexes_to_be_removed]
    root_key=new_key
    del tree[root_key]
    tree=rename_root(tree, new_key)

    return insert_children_in_tree(tree)

df = pandas.read_csv('mcmc_samples.csv')
posteriorr = df['posterior'].values.tolist()
modetree = ((df[["tree"]]).iloc[posteriorr.index(max(posteriorr))])[0]
ff = open("MAPadd.txt", "w")
ff.write(str((df[["add"]]).iloc[posteriorr.index(max(posteriorr))][0] ) +  "\n")
ff.close
leaves = sorted(list(set([leaf for node in (df['descendant_sets'].tolist()[0].split('-')) for leaf in node.split('.')])))

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
    bdat=pandas.DataFrame.from_records(branches)
    adat=pandas.DataFrame.from_records(admixtures)
    bresults=[]
    aresults=[]
    for n,(mean) in enumerate(bdat.mean()):
        bresults.append(('c'+str(n+1), mean))
    if len(admixtures[0])>0:
        for n,(mean) in enumerate(adat.mean()):
            aresults.append(('ax'+str(n+1), mean))
    return bresults, aresults

chosen_topology = modetree
temptopology = chosen_topology.split(';')[0]
for i in leaves[::-1]:
    chosen_topology = i + "=" + chosen_topology
branches_intervals, admixture_proportion_intervals=branch_and_proportion_quantiles([chosen_topology])

branch_names=[branches_interval[0] for branches_interval in branches_intervals]
admixture_names=[ad[0] for ad in admixture_proportion_intervals]

tree = identifier_to_tree(temptopology, leaves=generate_predefined_list_string(deepcopy(leaves)), branch_lengths=generate_predefined_list_string(deepcopy(branch_names)), admixture_proportions=generate_predefined_list_string(deepcopy(admixture_names)))

f = open("MAPtree.txt", "w")

for node in tree:
    if (tree[node][2] is not None):
        for adds in admixture_proportion_intervals:
            for branchh in branches_intervals:
                if tree[node][2] == adds[0]:
                    if tree[node][3] == branchh[0]:
                        f.write( node + " " + tree[node][0] + " "+ str(branchh[1]) + " "+ str(adds[1])+"\n")
                    if tree[node][4] == branchh[0]:
                        f.write( node + " " + tree[node][1] + " "+ str(branchh[1]) + " "+ str(1-adds[1])+"\n")
    else:
        for branchh in branches_intervals:
            if tree[node][3] == branchh[0]:
                        f.write( node + " " + tree[node][0] + " "+ str(branchh[1]) + " 1.00"+"\n")
f.close