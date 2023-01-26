from mcmc_proposals import addadmix_class, deladmix_class, sliding_regraft_class, rescale_class, rescale_admixtures_class, rescale_constrained_class, rescale_add_class
from numpy.random import choice
from Rtree_operations import get_number_of_admixes
#ALLGOOD
class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        if no_nodes==2:
            self.n+=1
            return ['x'+str(self.n)+a for a in ['a','b']]
        elif no_nodes==1:
            self.n+=1
            return 'x'+str(self.n)
        else:
            return ''
    
def initialize_proposals(proposals):
    all_props=[addadmix_class, deladmix_class, rescale_class, sliding_regraft_class, rescale_add_class, rescale_constrained_class,  rescale_admixtures_class]
    all_props_dic={cl.proposal_name:cl for cl in all_props}
    res=[]
    for proposal in proposals:
        res.append(all_props_dic[proposal]())
    return res
    
def draw_proposal(props, k, proportions):
    
    legal_indices=[i for i,prop in enumerate(props) if prop.require_admixture<=k]
    normaliser=sum([proportion for n,proportion in enumerate(proportions) if n in legal_indices])
    new_proportions=[float(proportion)/normaliser for n,proportion in enumerate(proportions) if n in legal_indices]
    
    chosen_index_i= choice(len(legal_indices), 1, p=new_proportions)[0]
    chosen_index=legal_indices[chosen_index_i]
    
    return chosen_index
    
def get_args2(names):
    args=[]
    if names:
        args.append(names)
    return args

class simple_adaptive_proposal(object):
    
    def __init__(self, proposals, proportions):
        self.props=initialize_proposals(proposals)
        self.proportions=proportions
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, x, pks={}):
        tree,add=x
        k=get_number_of_admixes(tree)
        index = draw_proposal(self.props, k, self.proportions)
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        proposal_input= self.props[index].input
        args=get_args2(names)
        if proposal_input=='add':
            new_add, forward, backward = self.props[index](add, *args, pks=pks)
            return (tree, new_add), forward, backward, 1.0
        if proposal_input=='tree':
            new_tree, forward, backward = self.props[index](tree, *args, pks=pks)
            return (new_tree, add), forward, backward, 1.0
        else:
            new_x, forward, backward = self.props[index](x, *args, pks=pks)
            return new_x, forward, backward, 1.0
        
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information
        