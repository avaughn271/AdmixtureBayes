from meta_proposal import simple_adaptive_proposal
import warnings

def get_proposals(all_proposals, all_proportions):
    
    thinned_proportions=[]
    thinned_proposals=[]
    
    for proposal, proportion in zip(all_proposals, all_proportions):
        if proportion> 1e-8:
            thinned_proportions.append(proportion)
            thinned_proposals.append(proposal)
    return thinned_proportions, thinned_proposals

def make_proposal(deladmix,
                  addadmix,
                  rescale,
                  regraft,
                  rescale_add,
                  rescale_admix,
                  rescale_admix_correction,
                  rescale_constrained,
                  rescale_marginally,
                  sliding_regraft,
                  sliding_rescale,
                  MCMC_chains,
                  cancel_preserve_root_distance=False,
                  no_add=False,
                  ):
    extras={}
    if cancel_preserve_root_distance:
        extras.update({'deladmix':{'preserve_root_distance':False}, 'addadmix':{'preserve_root_distance':False}})
    if no_add:
        extras['rescale_constrained']={'update_add':False}
        if rescale_add>0:
            rescale_add=0
            warnings.warn('removing the proposal rescale_add because no_add was used')

    all_proposals = ['deladmix', 'addadmix', 'rescale',
                     'regraft', 'rescale_add', 'rescale_admixtures',
                     'rescale_admix_correction',
                     'rescale_constrained', 'rescale_marginally',
                     'sliding_regraft', 'sliding_rescale']
    all_proportions = [deladmix, addadmix, rescale,
                       regraft, rescale_add, rescale_admix,
                       rescale_admix_correction,
                       rescale_constrained, rescale_marginally,
                       sliding_regraft, sliding_rescale]

    proportions, proposals = get_proposals(all_proposals, all_proportions)

    mp= [simple_adaptive_proposal(proposals, proportions, extras) for _ in range(MCMC_chains)]
    return mp
