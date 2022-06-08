from numpy.random import normal

def rescale(add, sigma=0.01, pks={}):
    pks['rescale_add_adap_param']=sigma
    new_add=add+normal()*sigma
    if new_add<0:
        return add,1,0 #rejecting by setting backward jump probability to 0.
    return new_add,1,1

class rescale_add_class(object):
    new_nodes=0
    proposal_name='rescale_add'
    adaption=True
    input='add'
    require_admixture=0
    reverse_require_admixture=0
    reverse='rescale_add'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)
