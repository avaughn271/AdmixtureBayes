from numpy import mean, array, delete

def fixed(num):
    return within_belt(num, 1e-6)

def within_belt(num, belt_size):
    if num<belt_size or num>1-belt_size:
        return True
    return False

def apply_filter(freqs,ns,names, locus_filter):
    new_freqs=[]
    new_ns=[]
    for i in range(ns.shape[1]):
        x=freqs[:,i]
        n=ns[:,i]
        if locus_filter(x,n, names):
            new_freqs.append(x)
            new_ns.append(n)
    return array(new_freqs).T, array(new_ns).T

class filter(object):
    
    def __init__(self, outgroup_name=''):
        self.outgroup_name=outgroup_name
        
    def __call__(self, freqs, pop_sizes, names=None):
        return True
        '''
        Return either True or false based on the numbers here. The default here returns true
        '''
    
    def apply_filter(self, xs, ns, names=None):
        new_freqs=[]
        new_ns=[]
        for i in range(ns.shape[1]):
            x=xs[:,i]
            n=ns[:,i]
            if self.__call__(x,n, names):
                new_freqs.append(x)
                new_ns.append(n)
        return array(new_freqs).T, array(new_ns).T
    
class snp_filter(filter):
    
    def __init__(self, outgroup_name=''):
        
        super(snp_filter, self).__init__(outgroup_name)
        
    def __call__(self, freqs, pop_sizes, names=None):
        if fixed(mean(freqs)):
            return False
        return True
        
        
class outgroup_other_filter(filter):
    
    def __init__(self, outgroup_name=''):
        
        super(outgroup_other_filter, self).__init__(outgroup_name)
        self.outgroup_index=None
        
    def __call__(self, freqs, pop_sizes, names=None):
        if self.outgroup_index is None:
            self.outgroup_index=self.get_outgroup_index(names)
        if fixed(freqs[self.outgroup_index]):
            return False
        others=mean(delete(freqs, self.outgroup_index))
        if fixed(others):
            return False
        return True
    
    def get_outgroup_index(self, names):
        if names is None:
            return 0
        assert self.outgroup_name, 'Unspecified outgroup'
        return next((n for n, e in enumerate(names) if e==self.outgroup_name))
    
class all_pops(filter):
    
    def __init__(self, outgroup_name=''):
        super(all_pops, self).__init__(outgroup_name)
        
    def __call__(self, freqs, pop_sizes, names=None):
        return all((not fixed(freq) for freq in freqs))
    
    
        
class outgroup_filter(filter):
    
    def __init__(self, outgroup_name=''):
        
        super(outgroup_filter, self).__init__(outgroup_name)
        self.outgroup_index=None
        
    def __call__(self, freqs, pop_sizes, names):
        if self.outgroup_index is None:
            self.outgroup_index=self.get_outgroup_index(names)
        if fixed(freqs[self.outgroup_index]):
            return False
        return True
    
    def get_outgroup_index(self, names):
        if names is None:
            return 0
        assert self.outgroup_name, 'Unspecified outgroup'
        return next((n for n, e in enumerate(names) if e==self.outgroup_name))
    
    