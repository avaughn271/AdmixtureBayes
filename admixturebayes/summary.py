from Rtree_operations import get_number_of_admixes, get_all_branch_lengths

class Summary(object):
       
    def __init__(self, name, pandable=True, output='double'):
        self.name=name
        self.pandable=pandable
        self.output=output
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return self.name+'= '+str(output)
        
class s_basic_tree_statistics(Summary):
    
    def __init__(self, function, name, output='double', args_extra=[]):
        super(s_basic_tree_statistics, self).__init__(name, output=output)
        self.function=function
        self.args_extra=args_extra
        
    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return self.function(tree, *self.args_extra)
    
class s_no_admixes(Summary):
    
    def __init__(self):
        super(s_no_admixes,self).__init__('no_admixes', output='integer')

    def __call__(self, **kwargs):
        old_tree=kwargs['tree']
        return get_number_of_admixes(old_tree)

class s_total_branch_length(Summary):

    def __init__(self):
        super(s_total_branch_length,self).__init__('total_branch_length')

    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return sum(get_all_branch_lengths(tree))

class s_variable(Summary):
    
    def __init__(self, variable, pandable=True, output='double'):
        super(s_variable, self).__init__(variable, pandable, output)

    def __call__(self, **kwargs):
        if self.name not in kwargs:
            return None
        return kwargs[self.name]
    
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
