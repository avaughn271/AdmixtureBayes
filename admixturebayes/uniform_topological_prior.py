from Rtree_operations import get_number_of_admixes, get_number_of_leaves
from math import log

class uniform_prior(object):
    
    def __init__(self, leaves):
        self.leaves=leaves
        self.admixture_probabilities=[]
        self.p_a_counts={}
        self.dic_of_counts={(1,0,0,0):1}
        
    def probability(self,tree=None, admixtures=None):
        if admixtures is None:
            A=get_number_of_admixes(tree)
        else:
            A=admixtures
        if len(self.admixture_probabilities)>A:
            return log(self.admixture_probabilities[A])-log(2)*A
        else:
            val=self.get_new_probability(A)
            return log(val)-log(2)*A
        
    def get_new_probability(self, A):
        sum=0
        for p in range(0, int(self.leaves/2)+1):
            if (p,A) not in self.p_a_counts:
                res=self.get_specific_count((self.leaves,p,A,0))
                self.p_a_counts[(p,A)]=res
            else:
                res=self.p_a_counts[(p,A)]
            sum+=res
        return 1.0/sum
        
    def get_specific_count(self,tuple_of_interest):
        if tuple_of_interest in self.dic_of_counts:
            return self.dic_of_counts[tuple_of_interest]
        L,P,A,E=tuple_of_interest
        if L<=0 or A<0 or E<0 or P<0 or E>A or 2*P>L:
            return 0
        term1=(E+1)*self.get_specific_count((L-1,P,A,E+1))
        term2=(L-2*P+1)*self.get_specific_count((L-1,P-1,A,E))
        term3=(L+2*P+3*A-2*E-2)*self.get_specific_count((L-1,P,A,E))
        term4=2*(P+1)/float(L*(L+1))*self.get_specific_count((L+1,P+1,A-1,E-1))
        term5=2*(P+1)*(P+2)/float(L*(L+1))*self.get_specific_count((L+1,P+2,A-1,E))
        term6=2*(P+1)*(L-2*P-1)/float(L*(L+1))*self.get_specific_count((L+1,P+1,A-1,E))
        term7=(L-2*P)*(L-2*P+1)/float(2*L*(L+1))*self.get_specific_count((L+1,P,A-1,E))
        s=term1+term2+term3+term4+term5+term6+term7
        self.dic_of_counts[tuple_of_interest]=s
        return s

 
def uniform_topological_prior_function(tree):
    up=uniform_prior(get_number_of_leaves(tree))
    return up.probability(tree=tree)
