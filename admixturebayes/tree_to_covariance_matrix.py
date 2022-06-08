
from copy import deepcopy
from numpy import zeros, ix_, outer

class Population:
    
    def __init__(self, weights,members):
        self.weights=weights
        self.members=members
        
    def remove_partition(self, weight):
        #removedprin "self.weight",self.weights
        n_w=[w*weight for w in self.weights]
        self.weights=[w*(1-weight) for w in self.weights]
        #removedprin "weight",weight
        #removedprin "self.weight",self.weights
        return Population(n_w, deepcopy(self.members))
    
    def merge_with_other(self, other):
        #removedprin "merge", self.pop, other.pop
        #removedprin "self",self.members, self.weights
        #removedprin "other", other.members, other.weights
     
        self.weights=[w+other.weights[other.members.index(m)] if m in other.members else w for m,w in zip(self.members,self.weights) ]
        tmpl=[(w,m) for w,m in zip(other.weights, other.members) if m not in self.members]
        if tmpl:
            x_w,x_m=list(zip(*tmpl))
            self.weights.extend(x_w)
            self.members.extend(x_m)
        
        #elf.pop={member:(self.pop.get(member,0.0)+other.pop.get(member,0.0)) for member in set(self.pop.keys()+other.pop.keys())}
        #removedprin self.pop
        
    def get_contributions_as_iterable(self, branch_length):
        #removedprin "calculating for the piece:"
        #removedprin self.pop
        #removedprin branch_length
        #list_of_res=[]
        for s1,w1 in self.pop.items():
            for s2,w2 in self.pop.items():
                #list_of_res.append(((s1,s2),w1*w2*branch_length))
                yield ((s1,s2),w1*w2*branch_length)
                #yield ((s1,s2),w1*w2*branch_length)
        #return list_of_res

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
        #self.covmat[ix_(indices,indices)]+=branch_length*outer(weights, weights)
    
def follow_branch(branch, population, covmat):
    branch_length=branch[1]
    if len(branch)==2:
        covmat.update(branch_length, population)
        return None, (branch[0], population), "coalesce"
    else:
        a_event=branch[2]
        a_direction=a_event[0]
        if a_direction==">":
            if len(a_event)==4:
                return branch, ((a_event[1],a_event[3]),None), "waiting for admixture in"
            else:
                covmat.update(branch_length, population)
                population.merge_with_other(a_event[4])
                return branch[:1]+branch[3:], (None, population), "received admixture in"
        else:
            covmat.update(branch_length, population)
            a_dest=a_event[1]
            a_proportion=a_event[2]
            a_identifier=a_event[3]
            new_pop=population.remove_partition(a_proportion)
            return branch[:1]+branch[3:], ((a_dest,a_identifier),(population,new_pop)), "moving admixture out"
        
def calc_covariance_from_generators(generators, node_to_index):
    res=zeros((len(node_to_index), len(node_to_index)))
    for generator in generators:
        for (a,b), x in generator:
            res[node_to_index[a],node_to_index[b]]+=x
    return res

            
            
        
def make_covariance(tree, nodes):

    node_to_index={node:i for i,node in enumerate(nodes)}
    pops=[Population([1.0],[node]) for node in nodes]
    covmat=Covariance_Matrix(node_to_index)
    #tree=deepcopy(tree)
    tree_dict={b[1]:deepcopy(b[:1]+b[2:]) for b in tree}
    tree_dict["r"]="stop"
    attaches=dict(  (branch_key,(pops[n],"ready")) for n,branch_key in enumerate(nodes) if (branch_key in nodes))

    waiting_in_list={}
    waiting_out_list={}

    #all_back_at_root=False
    while len(attaches)>=2:
        moves=0
        for branch_key, branch in tree_dict.items():
            if branch_key in attaches and attaches[branch_key][1]=="ready" and len(attaches)>=2:
                moves+=1
                if branch[1]<0:
                    return None
                new_branch, attachments, code= follow_branch(branch, attaches[branch_key][0], covmat)
                if code == "coalesce": #the branch has coalesced with another branch
                    del attaches[branch_key]
                    new_branch_key=attachments[0]
                    new_population=attachments[1]
                    if new_branch_key in attaches:
                        new_population.merge_with_other(attaches[new_branch_key][0])
                        attaches[new_branch_key]=(new_population, "ready")
                    else:
                        attaches[new_branch_key]=(new_population, "waiting for coalescer")
                elif code == "waiting for admixture in": 
                    identifier_key=attachments[0][1]
                    travel_from_branch=attachments[0][0]
                    if identifier_key in waiting_out_list:
                        new_branch[2].append(waiting_out_list[identifier_key])
                        attaches[branch_key]=(attaches[branch_key][0], "ready")
                        attaches[travel_from_branch]=(attaches[travel_from_branch][0], "ready")
                        del waiting_out_list[identifier_key]
                    else:
                        waiting_in_list[identifier_key]=travel_from_branch
                        attaches[branch_key]=(attaches[branch_key][0], "waiting for admixturer")
                elif code == "moving admixture out":
                    
                    travel_to_branch_key=attachments[0][0]
                    identifier_key=attachments[0][1]
                    new_population_stay=attachments[1][0]
                    new_population_move=attachments[1][1]
                    
                    if identifier_key in waiting_in_list:
                        #here there is already a branch waiting for the population that is being emitted.
                        destination_branch=tree_dict[travel_to_branch_key]
                        destination_branch[2].append(new_population_move)
                        del waiting_in_list[identifier_key]
                        attaches[branch_key]=(new_population_stay, "ready")
                        attaches[travel_to_branch_key]=(attaches[travel_to_branch_key][0], "ready")  
                    else:
                        #no branch is waiting and we will put this on hold waiting for it to happen
                        attaches[branch_key]=(new_population_stay, "waiting for admixture out")
                        waiting_out_list[identifier_key]=new_population_move
                elif code == "received admixture in":
                    _,new_population=attachments
                    attaches[branch_key]=(new_population, "ready")
                    
                tree_dict[branch_key]=new_branch
        if moves==0:
            return None #this means that the tree is illegal, because it does not progress.
    return covmat.covmat