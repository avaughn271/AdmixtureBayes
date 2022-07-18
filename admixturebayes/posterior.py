from prior import prior
from likelihood import likelihood, likelihood_from_matrix
from math import log
from numpy import loadtxt

def zero_likelihood(*args, **kwargs):
    return 0

class posterior_class(object):
    
    def __init__(self, 
                 emp_cov, 
                 M=10, 
                 p=0.5, 
                 use_skewed_distr=False, 
                 multiplier=None, 
                 nodes=None, 
                 use_uniform_prior=False, 
                 treemix=False,
                 add_variance_correction_to_graph=False,
                 prefix='',
                 variance_correction_file='',
                 prior_run=False,
                 unadmixed_populations=[],
                 r=0,
                 collapse_row='',
                 ):
        '''
        M can either be a float - the degrees of freedom in the wishart distribution or the constant variance in the treemix normal approximation of the covariance matrix.
        or M can be a matrix - the same size of emp_cov where each entry is the variance of that entry. 
        '''
        self.emp_cov=emp_cov
        self.M=M
        self.p=p
        self.base_r=r
        if treemix:
            print("error")
        else:
            if prior_run:
                self.lik=zero_likelihood
                self.likmat=zero_likelihood
            else:
                self.lik=likelihood
                self.likmat=likelihood_from_matrix
            
        self.use_skewed_distr=use_skewed_distr
        self.multiplier=multiplier
        self.nodes=nodes
        self.collapse_row=collapse_row
        self.use_uniform_prior=use_uniform_prior
        self.unadmixed_populations=unadmixed_populations
        
        if add_variance_correction_to_graph:
            if variance_correction_file:
                self.b=loadtxt(variance_correction_file)
            else:
                self.b=loadtxt(prefix+'variance_correction.txt')
            if multiplier:
                self.b*=multiplier
        else:
            self.b=None


        
    def __call__(self, x, pks={}, verbose=False,r=None):
        if r is None:
            r=self.base_r
        prior_value = prior(x,p=self.p,
                                              use_skewed_distr=self.use_skewed_distr,pks=pks,
                                              use_uniform_prior=self.use_uniform_prior,
                                              unadmixed_populations=self.unadmixed_populations,
                                              r=r)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        
        likelihood_value=self.lik(x, self.emp_cov,self.b, self.M, nodes=self.nodes, collapse_row=self.collapse_row, pks=pks)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value

def rescale_empirical_covariance(m):
    '''
    It is allowed to rescale the empirical covariance matrix such that the inferred covariance matrix takes values that are closer to the mean of the prior.
    '''
    
    n=m.shape[0]
    actual_trace=m.trace()
    expected_trace=log(n)/log(2)*n
    multiplier= expected_trace/actual_trace
    return m*multiplier, multiplier
