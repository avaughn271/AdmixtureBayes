from filters import outgroup_other_filter, filter, snp_filter, outgroup_filter, all_pops


def initor(a):
    if not isinstance(a, str):
        return a[0]
    else:
        return a
    
    
def make_filter(filter_type=['none','snp', 'outgroup_other','outgroup', 'all_pops'],
                outgroup_name='', covariance_pipeline=[]):
    filter_type=initor(filter_type)
    
    #taking care of special cases
    if 23 in covariance_pipeline and filter_type=='snp':
        filter_type='none'
    if 6 in covariance_pipeline and (7 in covariance_pipeline or 8 in covariance_pipeline) and filter_type=='snp':
        filter_type='none'
        
    if filter_type=='outgroup':
        return outgroup_filter(outgroup_name)
    elif filter_type=='outgroup_other':
        return outgroup_other_filter(outgroup_name)
    elif filter_type=='snp':
        return snp_filter()
    elif filter_type=='none':
        return filter()
    elif filter_type=='all_pops':
        return all_pops()
    else:
        assert False, 'Unexpected behaviour of make_filter because filter type is: '+str(filter_type)
    
    
    