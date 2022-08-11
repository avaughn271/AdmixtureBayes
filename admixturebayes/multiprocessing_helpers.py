
from multiprocessing import Queue, Process
from MCMC import basic_chain
from numpy import random

class basic_chain_class_as_process(object):
    
    def __init__(self, basic_chain_class):
        self.chain=basic_chain_class
        self.process= Process(target=self.chain)
        self.process.start()
    
    def start(self, p):
        self.chain.task_queue.put(p)
    
    def terminate(self):
        self.process.terminate()
        
    def complete(self):
        return self.chain.response_queue.get()
    
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
        return basic_chain(start_tree,
                           self.summaries, 
                           self.posterior_function, 
                           self.proposal, 
                           post, 
                           N, 
                           sample_verbose_scheme, 
                           overall_thinning, 
                           i_start_from, 
                           temperature, 
                           proposal_update,
                           multiplier)
        
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
        for chain, list_of_arguments in zip(self.group, list_of_lists_of_arguments):
            chain.start(list_of_arguments)
            counter+=1
        assert counter==len(self.group)
        return [chain.complete() for chain in self.group]
    
    def terminate(self):
        for chain in self.group:
            chain.terminate()