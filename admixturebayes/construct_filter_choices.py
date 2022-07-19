class filter(object):
    def __init__(self, outgroup_name=''):
        self.outgroup_name=outgroup_name
        
    def __call__(self, freqs, pop_sizes, names=None):
        return True
    
def make_filter(filter_type=['none','snp', 'outgroup_other','outgroup', 'all_pops'],
                outgroup_name='', covariance_pipeline=[]):
    return filter()