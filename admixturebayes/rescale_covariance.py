from math import log

def rescale_empirical_covariance(m, normalizer= ['min', 'max']):
    '''
    It is allowed to rescale the empirical covariance matrix such that the inferred covariance matrix takes values that are closer to the mean of the prior.
    
    In the nicest phylogeny with branches on each level 8+4+2, we expect n*log_2(n)
    In the worst phylogeny we get n+(n-1)+..+2 =(n+1)*n/2-1. 
    
    '''
    
    if not isinstance(normalizer, str):
        normalizer=normalizer[0]
    
    n=m.shape[0]
    actual_trace=m.trace()
    min_expected_trace=log(n)/log(2)*n
    max_expected_trace=n*(n+1)/2-1
    
    if normalizer=='min':
        multiplier= min_expected_trace/actual_trace
    elif normalizer=='max':
        multiplier= max_expected_trace/actual_trace
    else:
        assert False, 'normalizer not set properly. normalizer= '+str(normalizer)
    
    return m*multiplier, multiplier