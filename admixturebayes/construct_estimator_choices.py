from estimators import initor
from covariance_scaled import ScaledEstimator

def make_estimator(reduce_method, 
                   variance_correction, 
                   nodes, 
                   reducer,
                   add_variance_correction_to_graph=False,
                   prefix='',
                   save_variance_correction=True):

    variance_correction=initor(variance_correction)
    
    est=ScaledEstimator(reduce_method=reduce_method,
                 scaling='average_sum',
                 variance_correction=variance_correction,
                 add_variance_correction_to_graph=add_variance_correction_to_graph,
                 prefix_for_saving_variance_correction=prefix,
                 save_variance_correction=save_variance_correction)
    return(est)
