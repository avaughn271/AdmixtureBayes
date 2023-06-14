from numpy import random
from math import exp, log
from Rtree_operations import scale_tree_copy

from multiprocessing import Queue, Process

class basic_chain_class_as_process(object):
    
    def __init__(self, basic_chain_class):
        self.chain=basic_chain_class
        self.process= Process(target=self.chain)
        self.process.start()
    
class basic_chain_class(object):
    
    def __init__(self, summaries, posterior_function, proposal, resxeed):
        if resxeed is None: #if not applied the different chains will run with the same seexd, putting the assumptions of the model in danger.
            random.seed()
        elif isinstance(resxeed, int): #SxEEDINGDEBUG
            random.seed(resxeed)
        self.summaries=summaries
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
        return basic_chain(start_tree,  self.summaries,  self.posterior_function, self.proposal,  post,  N,  sample_verbose_scheme, i_start_from, temperature, proposal_update, multiplier)
        
class basic_chain_pool(object):
    
    def __init__(self, summaries, posterior_function, proposals, seeds, posterior_function_list=[]): #SxEEDDEBUG
        if seeds is None:
                seeds=[None]*len(proposals)
        if posterior_function_list:
            self.group = [basic_chain_class_as_process(
                basic_chain_class(summaries, post_function, proposal, seed)) for proposal, seed, post_function in
                zip(proposals, seeds, posterior_function_list)]
        else:
            self.group=[basic_chain_class_as_process(
                basic_chain_class(summaries, posterior_function, proposal,seed)) for proposal,seed in zip(proposals,seeds)]

    def order_calculation(self, list_of_lists_of_arguments):
        '''
        The list of arguments should math that of p in basic_chain_class.run_chain()
        '''
        counter=0
        for chainn, list_of_arguments in zip(self.group, list_of_lists_of_arguments):
            chainn.chain.task_queue.put(list_of_arguments)
            counter+=1
        return [chainn.chain.response_queue.get() for chainn in self.group]

def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    post_new=posterior_function(newx)

    likelihood_old, prior_old, desiredprior_old = post[:3]
    likelihood_new, prior_new, desiredprior_new = post_new[:3]
    #print(post)
    #print(post_new)
    
    if g2<=0 or j2<=0 or likelihood_new == -float('inf') or prior_new == -float('inf') or desiredprior_new == -float('inf'):
        logmhr=-float('inf')
    else:
        #logmhr=(likelihood_new-likelihood_old)/temperature+(prior_new-prior_old)+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh) #original
        #logmhr=(likelihood_new-likelihood_old+prior_new-prior_old)/temperature + log(g2)+log(j2)-log(j1)-log(g1)+log(Jh) #heat all
        logmhr=(likelihood_new-likelihood_old + desiredprior_old - desiredprior_new + prior_new - prior_old)/temperature+(desiredprior_new-desiredprior_old)+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh)
    if logmhr>100:
        mhr=float('inf')
    else:
        mhr=exp(logmhr)
            
    u=random.random()
    proposal.adapt(mhr)
    if u<mhr:
        return newx,post_new
    return x,post

def basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, i_start_from=0, 
                temperature=1.0, proposal_update=None, multiplier=None):
    if proposal_update is not None:
        proposal.node_naming.n=proposal_update['n']
    
    x=start_x
    if post is None:
        post=posterior_function(x)
    
    iteration_summary=[]
        
    for i in range(i_start_from,i_start_from+N):
        proposal_knowledge_scraper={}
        new_x,new_post=one_jump(x, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        if i%40==0:
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               tree=scale_tree_copy(new_x[0], 1.0/multiplier),
                                                               add=new_x[1],
                                                               posterior=new_post,
                                                               old_post=post,
                                                               old_tree=scale_tree_copy(x[0],1.0/multiplier),
                                                               iteration_number=i,**proposal_knowledge_scraper))
        x=new_x
        post=new_post
    
    return x, post, list(zip(*iteration_summary)),proposal.get_exportable_state()

def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    iteration=kwargs['iteration_number']
    res=[iteration]
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        save_bool = (save_num!=0) and (iteration % save_num==0)
        if save_bool:
            val=s(**kwargs)
            res.append(val)
        else:
            res.append(None)
    return res
    