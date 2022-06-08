from Rproposal_admix import addadmix_class, deladmix_class
from Rproposal_regraft import regraft_class
from Rproposal_rescale import rescale_class
from Rproposal_sliding_regraft import sliding_regraft_class, sliding_regraft_class_resimulate
from Rproposal_rescale_marginally import rescale_marginally_class
from Rproposal_sliding_rescale import sliding_rescale_class
from Rproposal_rescale_add import rescale_add_class
from Rproposal_rescale_constrained import rescale_constrained_class
from Rproposal_rescale_admix import rescale_admixtures_class
from Rproposal_rescale_admix_correction import rescale_admix_correction_class
from numpy.random import choice
from Rtree_operations import get_number_of_admixes
from math import exp

class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        #removedprin self.n
        if no_nodes==2:
            self.n+=1 
            return ['x'+str(self.n)+a for a in ['a','b']]
        elif no_nodes==1:
            self.n+=1  
            return 'x'+str(self.n)
              
        else:
            return ''

class basic_meta_proposal(object):
    
    def __init__(self):
        self.props=[addadmix_class(), deladmix_class(), regraft_class(), rescale_class(), rescale_marginally_class()]
        self.params=[None, None, None, [0.1], [0.1]]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, x, pks={}):
        tree, add=x
        index=choice(len(self.props),1)[0]
        if get_number_of_admixes(tree)==0 and index<=1:
            index=0
            backj=0.2
            forwj=0.4
        elif get_number_of_admixes(tree)==1 and index==1:
            backj=0.4
            forwj=0.2
        else:
            backj=0.2
            forwj=0.2
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj

    def adapt(self,mhr, u, post_new, post, temperature):
        pass
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
    
class no_admix_proposal(object):
    
    def __init__(self):
        self.props=[regraft_class(), rescale_class()]
        self.params=[None, [0.5]]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,0.5,0.5
    
    def adapt(self,mhr, u, post_new, post, temperature):
        pass
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']

class adaptive_proposal_no_admix(object):
    
    def __init__(self):
        self.props=[sliding_regraft_class(), rescale_class()]
        start_value_of_sigma=0.1
        start_value_of_slider=0.1
        self.node_naming=new_node_naming_policy()
        self.recently_called_type=None
        self.regraft_count=10
        self.rescale_count=10
        self.multiplier=10
        self.desired_mhr=0.0534
        self.alpha=0.9
        self.params=[[start_value_of_slider], [start_value_of_sigma]]

    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        backj=0.5
        forwj=0.5
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        self.recently_called_type=self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj
    
    def adapt(self,mhr, u, post_new, post, temperature):
        if self.recently_called_type=='rescale':
            new_val, self.rescale_count= standard_update(self.rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[1][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale')
            self.params[1]=[new_val]
        if self.recently_called_type=='sliding_regraft':
            new_val, self.regraft_count= standard_update(self.regraft_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[0][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='regraft_slider',
                                                         max_val=15.0)
            self.params[0]=[new_val]
            
            
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        #information['params']=self.params
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
        #self.params=information['params']

def get_random_proposal_without_deleting_empty(no_props, no_admixes):
    index=choice(no_props,1)[0]
    if no_admixes==0 and index<=1:
        index=0
        backj=1.0/no_props
        forwj=2.0/no_props
    elif no_admixes==1 and index==1:
        backj=2.0/no_props
        forwj=1.0/no_props
    else:
        backj=1.0/no_props
        forwj=1.0/no_props
    return index,backj, forwj

def get_args(names, params):
    args=[]
    if names:
        args.append(names)
    if params is not None:
        args.extend(params)
    return args

class simple_adaption(object):
    
    def __init__(self, start_value=0.1, count=10, multiplier=10, desired_mhr=0.234, alpha=0.9, maxval=15, name='adap'):
        self.value=start_value
        self.count=count
        self.multiplier=multiplier
        self.desired_mhr=desired_mhr
        self.alpha=alpha
        self.maxval=15
        self.name=name
        
    def get_value(self):
        return self.value
    
    def adapt(self, mhr):
        self.value, self.count= standard_update(self.count, 
                                                self.multiplier, 
                                                self.alpha, 
                                                self.value, 
                                                mhr, 
                                                desired_mhr=self.desired_mhr, 
                                                verbose=False,
                                                name=self.name)
    
def initialize_proposals(proposals, extras={}):
    all_props=[addadmix_class, deladmix_class, regraft_class, 
               rescale_class, sliding_regraft_class, sliding_regraft_class_resimulate,
               rescale_marginally_class, sliding_rescale_class, rescale_add_class,
               rescale_constrained_class,  rescale_admixtures_class, rescale_admix_correction_class]
    all_props_dic={cl.proposal_name:cl for cl in all_props}
    #removedprin all_props_dic
    res=[]
    for proposal in proposals:
        if proposal in extras:
            res.append(all_props_dic[proposal](**extras[proposal]))
        else:
            res.append(all_props_dic[proposal]())
    return res
    
def draw_proposal(props, k, proportions):
    
    legal_indices=[i for i,prop in enumerate(props) if prop.require_admixture<=k]    
    normaliser=sum([proportion for n,proportion in enumerate(proportions) if n in legal_indices])
    new_proportions=[float(proportion)/normaliser for n,proportion in enumerate(proportions) if n in legal_indices]
    
    chosen_index_i= choice(len(legal_indices), 1, p=new_proportions)[0]
    chosen_index=legal_indices[chosen_index_i]
    
    effect_of_chosen_index=props[chosen_index].admixture_change
    if effect_of_chosen_index!=0:
        legal_indices2=[i for i,prop in enumerate(props) if prop.require_admixture <= k+effect_of_chosen_index]    
        normaliser2=sum([proportion for n,proportion in enumerate(proportions) if n in legal_indices2])
        new_proportions2=[float(proportion)/normaliser2 for n,proportion in enumerate(proportions) if n in legal_indices2]
        reverse_type= props[chosen_index].reverse
        reverse_index= next((index for index, prop in enumerate(props) if prop.proposal_name==reverse_type))
        reverse_index_i= next((index_i for index_i, index in enumerate(legal_indices2) if index==reverse_index))
        return chosen_index, new_proportions[chosen_index_i], new_proportions2[reverse_index_i]
    else:
        return chosen_index, 1.96,1.96 #it is not really 1.96 and 1.96 but only the ratio between them matters and I like 1.96
    
def get_args2(names, adap_object):
    args=[]
    if names:
        args.append(names)
    if adap_object is not None:
        args.append(adap_object.get_value())
    return args    

class simple_adaptive_proposal(object):
    
    def __init__(self, proposals, proportions, extras={}):
        '''
        extras is of the form {proposal_name: {parameter1:argument1, parameter2:argument2, ...,},...} and is used
        to set any extra parameters in the proposals.
        '''
        self.props=initialize_proposals(proposals, extras)
        self.proportions=proportions
        self.adaps=[simple_adaption() if prop.adaption else None for prop in self.props]
        self.node_naming=new_node_naming_policy()
        self.recently_called_type=None
        self.recently_called_index=None
        
    def __call__(self, x, pks={}):
        tree,add=x
        k=get_number_of_admixes(tree)
        index, jforward, jbackward = draw_proposal(self.props, k, self.proportions)
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        self.recently_called_type=self.props[index].proposal_name
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
        
    def adapt(self, mhr, u, post_new, post, temperature):
        if self.props[self.recently_called_index].adaption:
            self.adaps[self.recently_called_index].adapt(mhr)
        
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        #information['params']=self.params
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
        #self.params=information['params']


class adaptive_proposal(object):
    
    def __init__(self, resimulate_regrafted_branch_length=False):
        self.props=[addadmix_class(), deladmix_class(), sliding_regraft_class(), rescale_constrained_class(), sliding_rescale_class(), rescale_add_class(), rescale_admixtures_class()]
        if resimulate_regrafted_branch_length:
            self.props[2]=sliding_regraft_class_resimulate(resimulate_regrafted_branch_length)
        start_value_of_sigma=0.1
        start_value_of_slider=0.1
        start_value_of_sliding_rescales=0.1
        start_value_of_sigma_add=0.1
        start_value_admix_rescale=0.1
        self.node_naming=new_node_naming_policy()
        self.recently_called_type=None
        self.regraft_count=10
        self.rescale_count=10
        self.sliding_rescale_count=10
        self.admix_rescale_count=10
        self.rescale_add_count=10
        self.multiplier=10
        self.desired_mhr=0.234
        self.alpha=0.9
        self.params=[None, None, [start_value_of_slider], [start_value_of_sigma], [start_value_of_sliding_rescales], [start_value_of_sigma_add], [start_value_admix_rescale]]
        
    def __call__(self, x, pks={}):
        tree, add=x
        index, backj, forwj = get_random_proposal_without_deleting_empty(len(self.props), get_number_of_admixes(tree))
        
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        self.recently_called_type=self.props[index].proposal_name
        args=get_args(names, self.params[index])
        
        if self.recently_called_type[-3:] == 'add':
            new_add, forward, backward =self.props[index](add, *args, pks=pks)
            return (tree,new_add),forward,backward,1,forwj,backj
        else:
            new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
            return (new_tree,add),forward,backward,1,forwj,backj

    def adapt(self, mhr, u, post_new, post, temperature):
        if self.recently_called_type=='rescale_constrained':
            new_val, self.rescale_count= standard_update(self.rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[3][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale')
            self.params[3]=[new_val]
        if self.recently_called_type=='sliding_regraft':
            new_val, self.regraft_count= standard_update(self.regraft_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[2][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='regraft_slider',
                                                         max_val=15.0)
            self.params[2]=[new_val]
            
        if self.recently_called_type=='sliding_rescale':
            new_val, self.sliding_rescale_count= standard_update(self.sliding_rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[4][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='sliding_rescale',
                                                         max_val=15.0)
            self.params[4]=[new_val]
        if self.recently_called_type == 'rescale_add':
            new_val, self.rescale_add_count= standard_update(self.rescale_add_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[5][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale_add',
                                                         max_val=15.0)
            self.params[5]=[new_val]
        if self.recently_called_type == 'rescale_admixtures':
            new_val, self.rescale_add_count= standard_update(self.admix_rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[6][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale_admix',
                                                         max_val=15.0)
            self.params[6]=[new_val]
            
            
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        #information['params']=self.params
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
        #self.params=information['params']
        
        
def standard_update(count, multiplier, alpha, old_value, mhr, desired_mhr=0.234, verbose=False, max_val=float('inf'), name='value'):
    count+=1
    gamma=multiplier/count**alpha
    change=exp(gamma*(min(1.0,mhr)-desired_mhr))
    value=old_value*change
    value=min(value, max_val)
    if verbose:
        print('old_'+name+'=',old_value)
        print('mhr=',mhr)
        print('count=',count)
        print('gamma=', gamma)
        print('multiplier=', change)
        print('new_'+name+'=', value)
    return value,count