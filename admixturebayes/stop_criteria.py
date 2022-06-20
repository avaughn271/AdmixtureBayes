from subprocess import call
import os


class stop_criteria(object):
    

    
    def __init__(self, 
                        frequency=20000, 
                        summaries=['no_admixes','average_branch_length','add','descendant_sets'], 
                        outfile='tmp_stop_criteria.txt',
                        verbose_level='normal',
                        burn_in=0.5,
                        topological_threshold=200,
                        continuous_threshold=200,
                        Rscript_path='Rscript',
                 ):
        self.counter=0
        self.frequency=frequency
        self.summaries=summaries
        self.non_topological_summaries=len(summaries)
        self.outfile=outfile
        self.topological_threshold = topological_threshold
        self.continuous_threshold = continuous_threshold
        self.Rscript_path = Rscript_path
        self.burn_in = burn_in
        self.verbose_level=verbose_level
        if not topological_threshold<0:
            self.summaries.extend(['Zero_Ntree','random_Ntree','mode_Ntree'])
        
        
    def __call__(self, cum_iterations, filename):
        if self.check_yet(cum_iterations):
            return self.stop_yet(filename)
        return False
        
    def check_yet(self, cum_iterations):
        if int(cum_iterations/self.frequency)>self.counter:
            self.counter+=1
            return True
        return False
    
    def stop_yet(self, filename):
        dir = os.path.dirname(__file__)
        command=[self.Rscript_path,os.path.join(dir, 'ESS.R'), filename, str(self.burn_in), self.outfile,  str(self.verbose_level)]+self.summaries
        #if self.verbose_level=='normal':
        #    print(command)
        call(command)
        return self.check_outfile()
    
    def check_outfile(self):
        with open(self.outfile, 'r') as f:
            f.readline()
            for n,lin in enumerate(f.readlines()):
                ess=float(lin.split()[-1])
                name=lin.split()[0]
                #if self.verbose_level=='normal':
                #    print(n, name, ess)
                if n<self.non_topological_summaries:
                    if ess < self.continuous_threshold:
                        return False
                else:
                    if ess < self.topological_threshold:
                        return False
        return True
