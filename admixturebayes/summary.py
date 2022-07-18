#import matplotlib.pyplot as plt
from Rtree_operations import get_number_of_admixes, get_all_branch_lengths
from tree_statistics import unique_identifier

from numpy import isfinite

def values_to_numbers(list_of_objects):
    res=[]
    seen={}
    count=0
    for element in list_of_objects:
        if element not in seen:
            count+=1
            seen[element]=count
        res.append(seen[element])
    return res

def thin_out_nans(x, index):
    return list(zip(*[(y,i) for y,i in zip(x,index) if isfinite(y)]))


class Summary(object):
       
    def __init__(self, name, pandable=True, output='double'):
        self.name=name
        self.pandable=pandable
        self.output=output
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return self.name+'= '+str(output)
    
    def summary_of_phylogeny(self, tree):
        return None
    
        
class s_basic_tree_statistics(Summary):
    
    def __init__(self, function, name, output='double', args_extra=[]):
        super(s_basic_tree_statistics, self).__init__(name, output=output)
        self.function=function
        self.args_extra=args_extra
        
    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return self.function(tree, *self.args_extra)

    def summary_of_phylogeny(self, tree):
        return self.function(tree, *self.args_extra)
        
    
class s_no_admixes(Summary):
    
    def __init__(self):
        super(s_no_admixes,self).__init__('no_admixes', output='integer')

    def __call__(self, **kwargs):
        old_tree=kwargs['tree']
        return get_number_of_admixes(old_tree)
    
    def summary_of_phylogeny(self, tree):
        return get_number_of_admixes(tree)
    
    

class s_total_branch_length(Summary):

    def __init__(self):
        super(s_total_branch_length,self).__init__('total_branch_length')

    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return sum(get_all_branch_lengths(tree))
    
    def summary_of_phylogeny(self, tree):
        return sum(get_all_branch_lengths(tree))
    
class s_average_branch_length(Summary):
    
    def __init__(self):
        super(s_average_branch_length,self).__init__('average_branch_length')

    def __call__(self, **kwargs):
        tree=kwargs['tree']
        all_branch=get_all_branch_lengths(tree)
        return sum(all_branch)/len(all_branch)
    
    def summary_of_phylogeny(self, tree):
        all_branch=get_all_branch_lengths(tree)
        return sum(all_branch)/len(all_branch)
    
class s_variable(Summary):
    
    def __init__(self, variable, pandable=True, output='double'):
        super(s_variable, self).__init__(variable, pandable, output)

    def __call__(self, **kwargs):
        if self.name not in kwargs:
            return None
        return kwargs[self.name]
    
class s_variable_recalculated(Summary):
    
    def __init__(self, variable, pandable=True, output='double', pks_function=None, args=[]):
        super(s_variable_recalculated, self).__init__(variable, pandable, output)
        self.pks_function=pks_function
        self.args=args

    def __call__(self, **kwargs):
        if self.name not in kwargs:
            return None
        return kwargs[self.name]
    
    def summary_of_phylogeny(self, tree):
        pks={}
        ad=self.pks_function(tree, *self.args, pks=pks)
        return pks[self.name]
    
class s_posterior(Summary):
    
    def __init__(self):
        super(s_posterior, self).__init__('posterior', output='double')

    def __call__(self, **kwargs):
        return sum(kwargs['posterior'][:2])
    
class s_likelihood(Summary):
    
    def __init__(self):
        super(s_likelihood, self).__init__('likelihood', output='double')

    def __call__(self, **kwargs):
        return kwargs['posterior'][0]
    
class s_prior(Summary):
    
    def __init__(self):
        super(s_prior, self).__init__('prior', output='double')

    def __call__(self, **kwargs):
        return kwargs['posterior'][1]
    
def _id(x):
    return x
    
class s_bposterior_difference(Summary):
    
    def __init__(self, summary=_id, summary_name=''):
        super(s_bposterior_difference, self).__init__(summary_name, output='double')
        self.summary=summary
        self.summary_name=summary_name
        
    def __call__(self, **kwargs):
        return self.summary(kwargs['proposed_posterior'])-self.summary(kwargs['old_post'])
    
class s_covariance_difference(Summary):
    
    def __init__(self, subset_of_matrix=_id, comparison=_id, summary_name=''):
        super(s_bposterior_difference, self).__init__(summary_name, output='double')
        self.subset_of_matrix=summary
        self.comparison=comparison
        
    def __call__(self, **kwargs):
        return self.comparison(self.subset_of_matrix(kwargs['proposed_posterior'][3]),self.subset_of_matrix(kwargs['old_post'][3]))
    
    
class s_tree_identifier(Summary):
    
    def __init__(self):
        super(s_tree_identifier,self).__init__('tree_identifier', output='string')
        
    def __call__(self, **kwargs):
        old_tree=kwargs['tree']
        return unique_identifier(old_tree)
    
    def summary_of_phylogeny(self, tree):
        return unique_identifier(tree)
    
    

class s_tree_identifier_new_tree(s_tree_identifier):
    
    def __init__(self):
        super(s_tree_identifier,self).__init__('tree_identifier_new_tree', output='string')
        
    def __call__(self, **kwargs):
        tree=kwargs['proposed_tree']
        return unique_identifier(tree)
        
        
 