from mcmc_proposals import addadmix_class, deladmix_class, sliding_regraft_class, rescale_class, rescale_admixtures_class, rescale_constrained_class, rescale_add_class
from numpy.random import choice
from Rtree_operations import get_number_of_admixes
from math import exp

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

class simple_adaption(object):
    
    def __init__(self, start_value=0.1, name='adap'):
        self.value=start_value
        self.count=10
        self.name=name
    
    def adapt(self, mhr):
        self.count+=1
        self.value=self.value*exp((10/self.count**0.9)*(min(1.0,mhr) - 0.234))
    
def initialize_proposals(proposals):
    all_props=[addadmix_class, deladmix_class, rescale_class, sliding_regraft_class, rescale_add_class, rescale_constrained_class,  rescale_admixtures_class]
    all_props_dic={cl.proposal_name:cl for cl in all_props}
    res=[]
    for proposal in proposals:
        res.append(all_props_dic[proposal]())
    return res
    
def draw_proposal(props, k):
    
    if k == 0: # which indices are legal will depend on the whether there is admixture or not
        legal_indices = [1,2,3,5,6]
        new_proportions = 1/5
    else:
        legal_indices = [0,1,2,3,4,5,6]
        new_proportions = 1/7

    chosen_index = choice(legal_indices, 1)[0]
    
    effect_of_chosen_index=props[chosen_index].admixture_change
    if effect_of_chosen_index != 0:
        if k + effect_of_chosen_index > 0:
            new_proportions2 = 1/7
        else:
            new_proportions2 = 1/5

        return chosen_index, new_proportions, new_proportions2
    else:
        return chosen_index, 1.96,1.96 #it is not really 1.96 and 1.96 but only the ratio between them matters and I like 1.96
    
def get_args2(names, adap_object):
    args=[]
    if names:
        args.append(names)
    if adap_object is not None:
        args.append(adap_object.value)
    return args 

class simple_adaptive_proposal(object):
    
    def __init__(self, proposals):
        self.props=initialize_proposals(proposals)
        self.adaps=[simple_adaption() if prop.adaption else None for prop in self.props]
        self.node_naming=new_node_naming_policy()
        self.recently_called_index=None
        
    def __call__(self, x, pks={}):
        tree,add=x
        k=get_number_of_admixes(tree)
        index, jforward, jbackward = draw_proposal(self.props, k)
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        self.recently_called_index=index
        proposal_input= self.props[index].input
        args=get_args2(names, self.adaps[index])
        
        if proposal_input=='add':
            new_add, forward, backward = self.props[index](add, *args, pks=pks)
            return (tree, new_add), forward, backward, 1.0, jforward, jbackward
        if proposal_input=='tree':
            new_tree, forward, backward = self.props[index](tree, *args, pks=pks)
            return (new_tree, add), forward, backward, 1.0, jforward, jbackward
        else:
            new_x, forward, backward = self.props[index](x, *args, pks=pks)
            return new_x, forward, backward, 1.0, jforward, jbackward
        
    def adapt(self, mhr):
        if self.props[self.recently_called_index].adaption:
            self.adaps[self.recently_called_index].adapt(mhr)
        
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information
        