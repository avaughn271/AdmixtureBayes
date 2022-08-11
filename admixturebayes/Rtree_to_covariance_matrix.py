from copy import deepcopy
from numpy import zeros, ix_, outer

from Rtree_operations import node_is_non_admixture, get_leaf_keys

class Covariance_Matrix():
    
    def __init__(self, nodes_to_index):
        self.ni=nodes_to_index
        self.covmat=zeros((len(nodes_to_index), len(nodes_to_index)))
    
    def get_indices(self, nodes):
        return [self.ni[n] for n in nodes]
    
    def get_addon(self, branch_length, weights):
        return branch_length*outer(weights, weights)
    
    def update(self, branch_length, population):
        indices=self.get_indices(population.members)
        self.covmat[ix_(indices,indices)]+=self.get_addon(branch_length, population.weights)
        
    def get_matrix(self):
        return self.covmat
    
Covariance_Matrix2=Covariance_Matrix

class Population:
    
    def __init__(self, weights,members):
        self.weights=weights
        self.members=members
        
    def get_weight(self, member):
        for m,w in zip(self.members, self.weights):
            if m==member:
                return w
        
    def subset_of_the_candidates(self,candidates):
        if any(( (cand in self.members) and (self.get_weight(cand)>1e-7) for cand in candidates)):
            if any((   (cand not in self.members) or (self.get_weight(cand)<1.0-1e-7) for cand in candidates)):
                return 'partly'
            else:
                return 'all'
        return 'none'
        
    def get_population_string(self, min_w, keys_to_remove=[]):
        return '.'.join(sorted([m for m,w in zip(self.members,self.weights) if w>0.0 and m not in keys_to_remove]))
        
    def remove_partition(self, weight):
        n_w=[w*weight for w in self.weights]
        self.weights=[w*(1-weight) for w in self.weights]
        return Population(n_w, deepcopy(self.members))
    
    def merge_with_other(self, other):
     
        self.weights=[w+other.weights[other.members.index(m)] if m in other.members else w for m,w in zip(self.members,self.weights) ]
        tmpl=[(w,m) for w,m in zip(other.weights, other.members) if m not in self.members]
        if tmpl:
            x_w,x_m=list(zip(*tmpl))
            self.weights.extend(x_w)
            self.members.extend(x_m)
        return self

def leave_node(key, node, population, covmat):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch_length=node[4], population=new_pop, covmat=covmat, dependent='none')]

def leave_node_and_check_admixtures(key, node, population, covmat):
    if node_is_non_admixture(node):
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat)],[]
    else:
        members=population.members[:]
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch_length=node[4], population=new_pop, covmat=covmat, dependent='none')], members

def follow_branch(parent_key, branch_length, population, covmat, dependent="none"):
    covmat.update(branch_length, population)
    return parent_key, population, dependent

def _add_to_waiting(dic,add,tree):
    key,pop,dep=add
    if key in dic:#this means that the node is a coalescence node
        dic[key][0][1]=pop
        dic[key][1][1]=dep
    else:
        if key=='r' or node_is_non_admixture(tree[key]):
            dic[key]=[[pop,None],[dep,"empty"]]
        else:
            dic[key]=[[pop],[dep]]
    return dic

def _full_node(key,dic):
    if key in dic:
        for dep in dic[key][1]:
            if dep=="empty":
                return False
        return True
    return False

def _merge_pops(pops):
    if len(pops)==1:
        return pops[0]
    else:
        return pops[0].merge_with_other(pops[1])

def _thin_out_dic(dic, taken):
    ready_nodes=[]
    for key,[pops, deps] in list(dic.items()):
        full_node=True
        for dep in deps:
            if dep is None or not (dep=="none" or _full_node(dep,dic) or dep in taken):
                full_node=False
            else:
                pass
        if full_node:
            taken.append(key)
            ready_nodes.append((key,_merge_pops(pops)))
            del dic[key]
    return dic,ready_nodes
                
def make_covariance(tree, node_keys=None, old_cov=False):
    if node_keys is None:
        node_keys=sorted(get_leaf_keys(tree))
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=list(zip(node_keys,pops))
    covmat=Covariance_Matrix2({node_key:n for n,node_key in enumerate(node_keys)})
    if old_cov:
        covmat=Covariance_Matrix2({node_key:n for n,node_key in enumerate(node_keys)})
    waiting_nodes={}
    taken_nodes=[]
    while True:
        for key,pop in ready_nodes:
            upds=leave_node(key, tree[key], pop, covmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return covmat.get_matrix()

class dummy_covmat(object):
    
    def update(self, *args, **kwargs):
        pass
        

def get_populations(tree, min_w=0.0, keys_to_include=None):
    
    node_keys=sorted(get_leaf_keys(tree))
    if keys_to_include is None:
        keys_to_remove=[]
    else:
        keys_to_remove=list(set(node_keys)-set(keys_to_include))
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=list(zip(node_keys,pops))
    waiting_nodes={}
    taken_nodes=[]
    covmat=dummy_covmat()
    pop_strings=[]
    while True:
        for key,pop in ready_nodes:
            pop_strings.append(pop.get_population_string(0.0, keys_to_remove))
            upds=leave_node(key, tree[key], pop, covmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            big_pop=ready_nodes[0][1]
            pop_strings.append(big_pop.get_population_string(0.0, keys_to_remove))
            break
    if '' in pop_strings:
        pop_strings.remove('')
    return sorted(list(set(pop_strings)))

def get_populations_string(tree, min_w=0.0, keys_to_include=None):
    return '-'.join(get_populations(tree, 0.0, keys_to_include))
