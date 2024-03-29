from copy import deepcopy
from numpy import zeros, ix_, outer
import pandas
import sys

def node_is_non_admixture(node):
    return (node[1] is None)

def get_leaf_keys(tree):
    res=[]
    for key, node in list(tree.items()):
        if (node[1] is None and node[5] is None):
            res.append(key)
    return res

class Covariance_Matrix():
    
    def __init__(self, nodes_to_index):
        self.ni=nodes_to_index
        self.covmatr=zeros((len(nodes_to_index), len(nodes_to_index)))
    
    def get_indices(self, nodes):
        return [self.ni[n] for n in nodes]
    
    def get_addon(self, branch_length, weights):
        return branch_length*outer(weights, weights)
    
    def update(self, branch_length, population):
        indices=self.get_indices(population.members)
        self.covmatr[ix_(indices,indices)]+=self.get_addon(branch_length, population.weights)

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
        
    def get_population_string(self, keys_to_remove=[]):
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
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, dependent='none'), follow_branch(parent_key=node[1],branch_length=node[4], population=new_pop, covmat=covmat, dependent='none')]

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
                
def make_covariance(tree):
    node_keys=sorted(get_leaf_keys(tree))
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=list(zip(node_keys,pops))
    covmat=Covariance_Matrix({node_key:n for n,node_key in enumerate(node_keys)})
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
    return covmat.covmatr

def getcovariance(inputt, outputt, addfile, add_Corr, covmatrixx):
    df = pandas.read_csv(inputt, header = None)
    Column1 = []
    Column2 = []
    Column3 = []
    Column4 = []

    treee = {}
    for i in range(df.shape[0]):
        Column1.append(((df.loc[i])[0]).split(" ")[0])
        Column2.append((df.loc[i])[0].split(" ")[1])
        Column3.append(((df.loc[i])[0]).split(" ")[2])
        Column4.append((df.loc[i])[0].split(" ")[3])
    UniqueNodes = pandas.unique(Column1 + Column2)
    rootname = "r"

    for j in UniqueNodes:
        if not(j in Column1):
            rootname = j

    for i in range(len(Column1)):
        if Column1[i] == rootname:
            Column1[i] = "r"
        if Column2[i] == rootname:
            Column2[i] = "r"
    
    UniqueNodes = pandas.unique(Column1 + Column2)

    for j in UniqueNodes:
        if not (j in Column2):
            treee[j] = [ Column2[Column1.index(j)], None, None,  float(Column3[Column1.index(j)]), None, None, None ]
        elif not(j in Column1):
            pass # this is root
        elif Column1.count(j) > 1:
            indic = []
            for idnum, listval in enumerate(Column1):
                if listval == j:
                    indic.append(idnum)
            treee[j] =  [ Column2[indic[0]],
            Column2[indic[1]],
            float(Column4[indic[0]]), 
            float(Column3[indic[0]]), 
            float(Column3[indic[1]]),
            Column1[Column2.index(j)],  None ]
        else:
            indic = []
            for idnum, listval in enumerate(Column2):
                if listval == j:
                    indic.append(idnum)
            treee[j] =  [ Column2[Column1.index(j)], None, None,  float(Column3[Column1.index(j)]), None, Column1[indic[0]],  Column1[indic[1]] ]

    df = pandas.read_csv(addfile, header = None)
    b  = df[0][0]

    file1 = open(covmatrixx, 'r')
    Lines = file1.readlines()
    a = 1.0
    if add_Corr:
        a = float(((Lines[len(Lines)-1])[11:]))

    (pandas.DataFrame((make_covariance(treee) +  b/a))).to_csv(outputt,  header=False, index = False)

for j in range(0, int(sys.argv[1])):
    getcovariance("TemporaryFiles/MAPtree" + str(j+1) + ".txt", "TemporaryFiles/MAPCov" + str(j+1) + ".txt", "TemporaryFiles/MAPadd" + str(j+1) + ".txt", True,  "TemporaryFiles/covariance_matrix" + str(j+1) + ".txt")

    for iiiiiii in range(100):
        getcovariance('TemporaryFiles/Tree' + str(1+iiiiiii) + "_" + str(j+1) + '.txt', 
        'TemporaryFiles/Cov' + str(1+iiiiiii) + "_" + str(j+1) + '.txt', 
        'TemporaryFiles/add' + str(1+iiiiiii) + "_" + str(j+1) + '.txt', True, "TemporaryFiles/covariance_matrix" + str(j+1) + ".txt")