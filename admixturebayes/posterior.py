from prior import prior
from numpy import loadtxt
from Rtree_to_covariance_matrix import make_covariance
from covariance_scaled import reduce_covariance
from scipy.stats import wishart
from numpy.linalg import LinAlgError

def likelihood(x, emp_cov, b, M=12,nodes=None, collapse_row='', pks={}):
    tree, add= x
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    par_cov=make_covariance(tree, nodes)
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    if collapse_row:
        n=len(nodes)-1
        par_cov=reduce_covariance(par_cov, n)
    if b is not None:
        par_cov+=b
    pks['covariance']=par_cov
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    try:
        d=wishart.logpdf(emp_cov, df=M, scale= (par_cov+add)/M)
    except (ValueError, LinAlgError) as e:
        return -float("inf")
    return d

class posterior_class(object):
    
    def __init__(self, 
                 emp_cov, 
                 M=10, 
                 p=0.5, 
                 multiplier=None, 
                 nodes=None, 
                 prefix='',
                 variance_correction_file='',
                 r=0,
                 collapse_row='',
                 ):
        '''
        M can either be a float - the degrees of freedom in the wishart distribution or the constant variance in the normal approximation of the covariance matrix.
        or M can be a matrix - the same size of emp_cov where each entry is the variance of that entry. 
        '''
        self.emp_cov=emp_cov
        self.M=M
        self.p=p
        self.base_r=r
        self.lik=likelihood
            
        self.multiplier=multiplier
        self.nodes=nodes
        self.collapse_row=collapse_row
        
        if variance_correction_file:
            self.b=loadtxt(variance_correction_file)
        else:
            self.b=loadtxt(prefix+'variance_correction.txt')
        if multiplier:
            self.b*=multiplier


        
    def __call__(self, x, pks={}, verbose=False,r=None):
        if r is None:
            r=self.base_r
        prior_value = prior(x,p=self.p, r=r)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        
        likelihood_value=self.lik(x, self.emp_cov,self.b, self.M, nodes=self.nodes, collapse_row=self.collapse_row, pks=pks)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        return likelihood_value, prior_value

