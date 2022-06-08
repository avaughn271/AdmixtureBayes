from prior import prior
from likelihood import likelihood, n_mark, likelihood_from_matrix, likelihood_treemix, likelihood_treemix_from_matrix
from scipy.stats import norm
from math import log
from generate_prior_trees import generate_phylogeny
from tree_statistics import identifier_to_tree_clean
from Rtree_operations import get_number_of_leaves, get_number_of_admixes
from Rtree_to_covariance_matrix import make_covariance
from numpy import median, amin, amax, loadtxt
from numpy.linalg import norm


def initialize_posterior2(emp_cov=None, 
                         true_tree=None, 
                         M=None, 
                         use_skewed_distr=False, 
                         p=0.5, 
                         rescale=False, 
                         model_choice=['empirical covariance',
                                       'true tree covariance',
                                       'wishart on true tree covariance',
                                       'empirical covariance on true tree',
                                       'no likelihood'],
                         simulate_true_tree=False,
                         true_tree_no_leaves=None,
                         true_tree_no_admixes=None,
                         nodes=None,
                         simulate_true_tree_with_skewed_prior=False,
                         reduce_cov=None,
                         add_outgroup_to_true_tree=False,
                         reduce_true_tree=False):
    
    if not isinstance(model_choice, str):
        model_choice=model_choice[0]
        
    if model_choice == 'no likelihood':
        return initialize_prior_as_posterior(), {}
        
    if (model_choice == 'true tree covariance' or 
        model_choice == 'wishart on true tree covariance' or
        model_choice == 'empirical covariance on true tree'):
        
        if simulate_true_tree:
            true_tree= generate_phylogeny(true_tree_no_leaves,
                               true_tree_no_admixes,
                               nodes,
                               simulate_true_tree_with_skewed_prior)
            
        elif isinstance(true_tree, str):
            if ';' in true_tree: #this means that the true tree is a s_tree
                true_tree_s=true_tree
                true_tree=identifier_to_tree_clean(true_tree_s)
            else:
                with open(true_tree, 'r') as f:
                    true_tree_s=f.readline().rstrip()
                true_tree=identifier_to_tree_clean(true_tree_s)
        
        
        
                
        true_tree=Rtree_operations.simple_reorder_the_leaves_after_removal_of_s1(true_tree)
                
        
        no_leaves = get_number_of_leaves(true_tree)
        no_admixes = get_number_of_admixes(true_tree)
        
        
        
        cov=make_covariance(true_tree)
        
        if reduce_cov is not None:
            pass 
        if reduce_true_tree is not None:
            true_tree=Rtree_operations.remove_outgroup(true_tree, reduce_true_tree)
            if reduce_true_tree=='s1' or reduce_true_tree==0:
                pass
        if emp_cov is not None:
            if isinstance(emp_cov, str):
                pass
    
    if M is None:
        M=n_mark(emp_cov)
    if rescale:
        emp_cov, multiplier = rescale_empirical_covariance(emp_cov)
        print('multiplier is', multiplier)
    def posterior(x,pks={}):
        #removedprin tot_branch_length
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood(x, emp_cov,M=M)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value
    if rescale:
        return posterior, multiplier
    return posterior

def initialize_treemix_posterior(emp_cov, variances, p=0.5, use_skewed_distr=False, multiplier=None, nodes=None, use_uniform_prior=False):
    def posterior(x,pks={}):
        #removedprin tot_branch_length
        #removedprin get_number_of_leaves(x[0]), emp_cov.shape[0]
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks, use_uniform_prior=use_uniform_prior)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood_treemix(x, emp_cov, variances=variances, nodes=nodes)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value
    if multiplier is not None:
        return posterior, multiplier
    return posterior


def initialize_posterior(emp_cov, M=10, p=0.5, use_skewed_distr=False, multiplier=None, nodes=None, use_uniform_prior=False):
    def posterior(x,pks={}):
        #removedprin tot_branch_length
        #removedprin get_number_of_leaves(x[0]), emp_cov.shape[0]
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks, use_uniform_prior=use_uniform_prior)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood(x, emp_cov,M=M, nodes=nodes)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value
    if multiplier is not None:
        return posterior, multiplier
    return posterior

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
                 collapse_row_rearrange=False
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
            self.lik=likelihood_treemix
            self.likmat=likelihood_treemix_from_matrix
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
        #print('posterior_nodes', self.nodes) #ANDREWDEBUG
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
        if verbose:
            print('empirical_matrix=', self.emp_cov)
            print('input_matrix=', pks['covariance']+x[1])
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value
    
    def get_likelihood_from_matrix(self, matrix, pks={}, verbose=False):
        val=self.likmat(matrix, self.emp_cov, self.b, self.M, pks=pks)
        if verbose:
            print('empirical_matrix=', self.emp_cov)
            print('input_matrix=', matrix)
        return val
    
    def get_max_likelihood(self, pks={}, verbose=False):
        val=self.likmat(self.emp_cov, self.emp_cov,None, self.M, pks=pks)
        if verbose:
            print('empirical_matrix=', self.emp_cov)
            print('input_matrix=', self.emp_cov)
        return val
    
    def alternative_emp_cov_likelihood(self, alternative_emp_cov, x, pks={}, verbose=False):
        val=self.lik(x,alternative_emp_cov,self.b, self.M, nodes=self.nodes, pks=pks)
        if verbose: 
            print('empirical_matrix=', self.emp_cov)
            print('input_matrix=', pks['covariance']+x[1])
        return val
    
    def get_non_empirical_max_likelihood(self, x, pks={}, verbose=False):
        p_cov=make_covariance(x[0])+x[1]
        if self.b is not None:
            p_cov+=self.b
        val=self.likmat(p_cov, p_cov, None, self.M, pks=pks)
        if verbose:
            print('empirical_matrix=', p_cov)
            print('input_matrix=', p_cov)
        return val
    
    def get_size_diff(self, x):
        t,add=x
        p_cov=make_covariance(t)+add
        if self.b is not None:
            p_cov+=self.b
        diffs=p_cov-self.emp_cov
        max_dif=amax(diffs)
        min_dif=amin(diffs)
        return median(p_cov/self.emp_cov), (max_dif+min_dif)/(abs(max_dif)+abs(min_dif)), norm(diffs)
        
def initialize_big_posterior(emp_cov, M=None, use_skewed_distr=False, p=0.5):
    if M is None:
        M=n_mark(emp_cov)
    def posterior(x,pks={}):
        #removedprin tot_branch_length
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood(x, emp_cov,M=M, pks=pks)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        prior_values=(pks['branch_prior'], pks['no_admix_prior'], pks['admix_prop_prior'], pks['top_prior'])
        covariance=pks['covariance']
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value, prior_values, covariance
    return posterior
        

def initialize_prior_as_posterior(p=0.5):
    def posterior(x,pks={}):
        #removedprin tot_branch_length
        prior_value=prior(x,p=p,pks=pks)
        if prior_value==-float('inf'):
            return prior_value
        pks['prior']=prior_value
        pks['likelihood']=0
        return 0,prior_value
    return posterior

def initialize_trivial_posterior():
    def posterior(x, pks={}):
        if isinstance(x, float):
            res= norm.logpdf(x)
            pks['prior']= res
            return res
#         elif isinstance(x, list) and isinstance(x[0], float):
#             res= multivariate_normal.logpdf(x)
#             pks['prior']= res
#             return res
        else:
            assert False, 'input in posterior was not recognizable.'
    return posterior

def print_pks(pks):
    for key, value in list(pks.items()):
        print('\t',key,':',value)

def call_post(posterior_function,tree, posterior_function_name='posterior', tree_name='tree'):
    pks={}
    print(posterior_function_name+'('+tree_name+')=', posterior_function(tree, pks=pks))
    print_pks(pks)
    
def rescale_empirical_covariance(m):
    '''
    It is allowed to rescale the empirical covariance matrix such that the inferred covariance matrix takes values that are closer to the mean of the prior.
    '''
    
    n=m.shape[0]
    actual_trace=m.trace()
    expected_trace=log(n)/log(2)*n
    multiplier= expected_trace/actual_trace
    return m*multiplier, multiplier