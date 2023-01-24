from numpy import random
from math import exp, log
from Rtree_operations import scale_tree_copy, get_number_of_admixes, get_all_branch_lengths, get_number_of_ghost_populations
from tree_statistics import unique_identifier_and_branch_lengths, get_admixture_proportion_string
from Rtree_to_covariance_matrix import get_populations_string
from multiprocessing import Queue, Process

class basic_chain_class_as_process(object):
    
    def __init__(self, basic_chain_class):
        self.chain=basic_chain_class
        self.process= Process(target=self.chain)
        self.process.start()
    
class basic_chain_class(object):
    
    def __init__(self, posterior_function, proposal, resxeed):
        if resxeed is None: #if not applied the different chains will run with the same seexd, putting the assumptions of the model in danger.
            random.seed()
        elif isinstance(resxeed, int): #SxEEDINGDEBUG
            random.seed(resxeed)
        self.posterior_function=posterior_function
        self.proposal=proposal
        self.task_queue= Queue()
        self.response_queue = Queue()

    def __call__(self):
        while True:
            input = self.task_queue.get()
            self.response_queue.put(self.run_chain(input))
            
    def run_chain(self, p):
        start_tree, post, N, sample_verbose_scheme, overall_thinning, i_start_from, temperature, proposal_update, multiplier = p
        return basic_chain(start_tree, self.posterior_function, self.proposal,  post,  N,  sample_verbose_scheme, i_start_from, temperature, proposal_update, multiplier)
        
class basic_chain_pool(object):
    
    def __init__(self, posterior_function, proposals, seeds, posterior_function_list=[]): #SxEEDDEBUG
        if seeds is None:
                seeds=[None]*len(proposals)
        if posterior_function_list:
            self.group = [basic_chain_class_as_process(
                basic_chain_class( post_function, proposal, seed)) for proposal, seed, post_function in
                zip(proposals, seeds, posterior_function_list)]
        else:
            self.group=[basic_chain_class_as_process(
                basic_chain_class(posterior_function, proposal,seed)) for proposal,seed in zip(proposals,seeds)]

    def order_calculation(self, list_of_lists_of_arguments):
        '''
        The list of arguments should math that of p in basic_chain_class.run_chain()
        '''
        counter=0
        for chainn, list_of_arguments in zip(self.group, list_of_lists_of_arguments):
            chainn.chain.task_queue.put(list_of_arguments)
            counter+=1
        return [chainn.chain.response_queue.get() for chainn in self.group]

def one_jump(x, post, temperature, posterior_function, proposal):
    newx,g1,g2,Jh,j1,j2=proposal(x, {})
    post_new=posterior_function(newx)

    likelihood_old, prior_old = post[:2]
    likelihood_new, prior_new = post_new[:2]
    
    if g2<=0 or j2<=0:
        logmhr=-float('inf')
    else:
        logmhr=(likelihood_new-likelihood_old)/temperature+(prior_new-prior_old)+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh)
    if logmhr>100:
        mhr=float('inf')
    else:
        mhr=exp(logmhr)
            
    u=random.random()
    proposal.adapt(mhr)
    if u<mhr:
        return newx,post_new
    return x,post

def basic_chain(start_x, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, i_start_from=0, 
                temperature=1.0, proposal_update=None, multiplier=None):
    
    proposal.node_naming.n=proposal_update['n']
    
    x=start_x
    
    iteration_summary=[]
        
    for i in range(i_start_from,i_start_from+N):
        new_x,new_post=one_jump(x, post, temperature, posterior_function, proposal)
        if i%40==0:
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme, i,
                                                               tree=scale_tree_copy(new_x[0], 1.0/multiplier),
                                                               add=new_x[1], posterior=new_post ))
        x=new_x
        post=new_post
    return x, post, list(zip(*iteration_summary)),proposal.get_exportable_state()

def _calc_and_print_summaries(sample_verbose_scheme,iteration, tree, add, posterior):
    res=[iteration]
    if len(sample_verbose_scheme) != 0:
        res.extend([sum(posterior[:2]),posterior[0],posterior[1],get_number_of_admixes(tree),add,
        sum(get_all_branch_lengths(tree)), get_number_of_ghost_populations(tree), get_populations_string(tree),
        unique_identifier_and_branch_lengths(tree), get_admixture_proportion_string(tree)])
    else:
        res.extend([None,None,None,None,None,None,None,None,None,None])
    return res
