import numpy as np

#maximizes the function, function in the interval [lower_limit,\infty).
def I_cant_believe_I_have_to_write_this_function_myself(function, lower_limit):
    old_x=lower_limit
    new_x=lower_limit*2
    old_y=function(old_x)
    lower_y=old_y
    new_y=function(new_x)
    sgn=1
    max_step_size_increases=20
    step_size_increases=0
    step_size=lower_limit
    for _ in range(20):
        c=0
        while new_y<old_y and c<100:
            old_x=new_x
            new_x+=sgn*step_size
            old_y=new_y
            if new_x<lower_limit:
                new_x=lower_limit
                new_y=lower_y
                break
            new_y=function(new_x)
            c+=1
            if c>10 and max_step_size_increases>step_size_increases:
                c=0
                step_size*=2
                step_size_increases+=1
        sgn=-sgn
        step_size*=0.5
        old_x=new_x
        new_x=old_x+sgn*step_size
        old_y=new_y
        new_y=function(new_x)
    return new_x

def variance_mean_based(sample_of_matrices, divisor=None, verbose_level='normal'):
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    var_wishart= np.var(sample_of_matrices, axis=0)
    r=mean_wishart.shape[0]
    var_rom_mean_wishart=np.square(mean_wishart)+np.outer(np.diag(mean_wishart),np.diag(mean_wishart))
    def penalty_function(df_l):
        df=df_l
        val=np.linalg.norm(var_wishart-var_rom_mean_wishart/df)**2
        return np.log(val)
    
    rval=I_cant_believe_I_have_to_write_this_function_myself(penalty_function, r)
    return rval

def initor(a):
    if not isinstance(a, str):
        return a[0]
    else:
        return a
    
class Estimator(object):
    
    def __init__(self,
                 reduce_also=True):
        self.reduce_also=reduce_also
        self.fitted_value=None
        self.initial_Sigma_generator=None
        self.initial_Sigma=None
        
    def get_reduce_index(self):
        return 0
    
    def __call__(self,xs,ns, names=None, extra_info={}):
        '''
        Should return the covariance estimate and possibly set the variable fitted_value
        '''
        pass
