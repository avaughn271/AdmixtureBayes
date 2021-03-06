from Rtree_operations import update_specific_branch_lengths
from copy import deepcopy
from numpy.random import normal
from Rtree_operations import get_leaf_keys, get_all_branches, node_is_non_admixture
from numpy import zeros, insert
from Rtree_to_covariance_matrix import Population, _add_to_waiting, _thin_out_dic
from scipy.linalg import svd

class Coefficient_Matrix():
    
    def __init__(self, nodes_to_index, branch_to_index, custom_list={}):
        self.ni=nodes_to_index
        self.bi=branch_to_index
        self.custom_list=custom_list
        if custom_list:
            self.cofmat=zeros((len(custom_list), len(branch_to_index)))
        else:
            self.cofmat=zeros((len(nodes_to_index), len(branch_to_index)))
        
    def update(self, branch, population):
        j=self.bi[branch]
        if self.custom_list:
            for pop1_index in range(len(population.members)):
                for pop2_index in range(pop1_index, len(population.members)):
                    pop1=population.members[pop1_index]
                    pop2=population.members[pop2_index]
                    key=None
                    if (pop1,pop2) in self.custom_list:
                        key=(pop1,pop2)
                    elif (pop2,pop1) in self.custom_list:
                        key=(pop2,pop1)
                    if key is not None:
                        self.cofmat[self.custom_list[key],j]=population.weights[pop1_index]*population.weights[pop2_index]
        else:
            for pop, w in zip(population.members, population.weights):
                self.cofmat[self.ni[pop],j]=w**2
    
    def get_matrix(self):
        return self.cofmat


def nullspace(A, atol=1e-13, rtol=0):
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

def get_orthogonal_branch_space(tree, add_one_column=True):
    cof,_, bi= make_coefficient_matrix(tree)
    if add_one_column:
        cof=insert(cof, cof.shape[1], 1, axis=1)
    ad=nullspace(cof)
    return ad, bi, cof.T
    
def make_coefficient_matrix(tree, node_keys=None, branch_keys=None):
    '''
    Instead of constructing the covariance matrix, this function calculates the coefficient matrix, C, to solve
    
    w=Cx
    
    where w is the diagonal of the covariance matrix and x is the vector of branch lengths. Hence, C depends on the admixture proportions and the topology.
    '''
    if node_keys is None:
        node_keys=sorted(get_leaf_keys(tree))
    if branch_keys is None:
        branch_keys=get_all_branches(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=list(zip(node_keys,pops))
    ni={node_key:n for n,node_key in enumerate(node_keys)}
    bi={branch:n for n,branch in enumerate(branch_keys)}
    cofmat=Coefficient_Matrix(ni,bi, get_all_pairs(node_keys))
    waiting_nodes={}
    taken_nodes=[]
    while True:
        for key,pop in ready_nodes:
            upds=leave_node(key, tree[key], pop, cofmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return cofmat.get_matrix(), ni,bi

def leave_node(key, node, population, cofmat):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch=(key,0), population=population, cofmat=cofmat)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch=(key,0), population=population, cofmat=cofmat, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch=(key,1), population=new_pop, cofmat=cofmat, dependent='none')]

def follow_branch(parent_key, branch, population, cofmat, dependent="none"):
    cofmat.update(branch, population)
    return parent_key, population, dependent

def get_all_pairs(nodes):
    all_nodes={}
    count=0
    for node_index in range(len(nodes)):
        for node_index2 in range(node_index, len(nodes)):
            node1=nodes[node_index]
            node2=nodes[node_index2]
            all_nodes[(node1, node2)]=count
            count+=1
    return all_nodes
    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in list(dic.items()):
        l_determined[i]=e
    return l_determined

def get_added_branch_pieces(org, param):
    return org.dot(normal(scale=param, size=org.shape[1]))

def rescale(x, sigma=0.01, pks={}, update_add=True):
    tree, add=x
    pks['rescale_constrained_adap_param']=sigma
    new_tree=deepcopy(tree)
    orgs, bi,_= get_orthogonal_branch_space(new_tree, add_one_column=update_add)
    #removedprin norm(orgs, axis=0)
    branches=reverse_dic_to_list(bi)
    branch_pieces= get_added_branch_pieces(orgs, sigma)
    if update_add:
        b=branch_pieces[:-1]
        #removedprin 'ADD!'
    else:
        #removedprin 'NO ADD!'
        b=branch_pieces
    new_tree= update_specific_branch_lengths(new_tree, branches, b, add=True)
    if new_tree is None:
        return x,1,0 #rejecting by setting backward jump probability to 0.
    if update_add:
        new_add=add+branch_pieces[-1]
    else:
        new_add=0
    if new_add<0:
        return x,1,0 #rejecting by setting backward jump probability to 0.
    return  (new_tree, new_add),1,1

class rescale_constrained_class(object):
    new_nodes=0
    proposal_name='rescale_constrained'
    adaption=True
    input='both'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=0
    reverse='rescale_constrained'
    
    
    def __init__(self, **kwargs):
        self.kwargs=kwargs
    
    def __call__(self,*args, **kwargs):
        kwargs.update(self.kwargs)
        #removedprin kwargs
        return rescale(*args, **kwargs)

