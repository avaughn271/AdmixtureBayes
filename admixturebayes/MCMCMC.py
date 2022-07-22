#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
from MCMC import basic_chain
from pathos.multiprocessing import freeze_support
from Rtree_operations import get_number_of_admixes
from posterior import no_admixes
import pandas as pd
from multiprocessing_helpers import basic_chain_pool
from numpy.random import choice, random
from math import exp
from itertools import chain
import time
import os

def _basic_chain_unpacker(args):
    return basic_chain(*args)

def MCMCMC(starting_trees, 
           posterior_function,
           summaries,
           temperature_scheme, 
           printing_schemes, 
           iteration_scheme, 
           overall_thinnings, 
           proposal_scheme, 
         n_arg, m_arg, verboseee,  cores=4,
           no_chains=None,
           numpy_seeds=None,
           multiplier= None,
           result_file=None,
           stop_criteria=None,
           posterior_function_list=[]):
    '''
    this function runs a MC3 using the basic_chain_unpacker. Let no_chains=number of chains. The inputs are
        starting_trees: a list of one or more trees that the chains should be started with
        posterior_function: one 'initialized'(see MCMC.py) unnormalized posterior function, that should be simulated from.
        summaries: a list of instances of realization of concrete classes of the superclass Summary. It is closely linked to printing_scheme, 
                   because they are only meaningful if specified in that.
        temperature_scheme: one instance of a class that has the functions:
                                    - get_temp(i): returns the the temperature for the i'th chain
                                    - update_temp(permut): updates the temperatures for each chain using the permutation, permut
        printing_scheme: a list of either one or no_chains dictionaries, where each dictionary gives the sample_verbose_scheme of basic_chain
        iteration_scheme: a list that sums to the total number of iterations, where each entry is the number of MH-steps before the next flipping
        overall_thinnings: an integer indicating how much should be skipped before each summary statistics is calculated
        proposal_scheme: a list of instances of classes that handles proposals and updates of said proposals. It has the basic functions
                                    - prop(x,pks): proposes the next tree(and returns proposal densities and statistics in pks)
                                    - adapt(mhr): updates the proposals based on the mhr ratio, mhr
                                    - extract_new_values(): after two chains have changed position in the flipping, this function gets the information that should be changed_later
                                    - wear_new_values(information): the new values from extract_new_values replaces the old values.
        #ps:     - a list of doubles in the interval [0,1] indicating the parameter of the geometrically distributed prior on the number of admixture events for each chain. Or:
        #        - a double in the interval [0,1] indicating the same parameter for all geometrically distributed prior on the number of admixture events.
                                    
    '''
    
    
    if no_chains is None:
        no_chains=cores
        
    if len(printing_schemes)==1:
        printing_schemes=[printing_schemes[0]]*no_chains
        
    df_result=None
    total_permutation=list(range(no_chains))
    xs = starting_trees

    freeze_support()
    
    #if numpy_sxeeds is None: #Sxeeddebug
    #    numpy_sxeeds=[None]*no_chains


    rs=[]
    ps=[posterior_function.p]
    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme, numpy_seeds)
    posteriors = [posterior_function(x) for x in xs]

    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]
    
    cum_iterations=0
    start_time=time.time()
    for no_iterations in iteration_scheme:
        if cum_iterations % 1000 == 0 and verboseee != "silent":
            print("Currently on iteration " +  str(cum_iterations) + " out of " + str(n_arg * m_arg))
        #letting each chain run for no_iterations:
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, printing_schemes, overall_thinnings, no_iterations, cum_iterations, proposal_updates, multiplier)
        if no_chains==1:#debugging purposes
            new_state=[_basic_chain_unpacker(next(iteration_object))]
        else:
            new_state = pool.order_calculation(iteration_object)
        xs, posteriors, df_add,proposal_updates = _unpack_everything(new_state, summaries, total_permutation)
        df_result=_update_results(df_result, df_add)
        if result_file is not None:
            if cum_iterations==0:
                start_data_frame(df_result, result_file)
            elif df_result.shape[0]>1000:
                add_to_data_frame(df_result, result_file)
                df_result=df_result[0:0]
        #making the mc3 flips and updating:
        if not posterior_function_list:
            xs, posteriors, permut, proposal_updates = flipping(xs, posteriors, temperature_scheme, proposal_updates,
                                                                rs, ps,
                                                                [posterior_function])  # trees, posteriors, range(len(trees)),[None]*len(trees)#
        total_permutation=_update_permutation(total_permutation, permut)
        cum_iterations+=no_iterations
        if stop_criteria is not None:
            if stop_criteria(cum_iterations, result_file):
                print("Stop criteria reached. MCMC terminating now.")
                if os.path.exists("trees_tmp.txt"):
                    os.remove("trees_tmp.txt")
                if os.path.exists("stop_criteria.txt"):
                    os.remove("stop_criteria.txt")
                break
            
    pool.terminate()
    if result_file is None:
        return df_result
        
def _update_permutation(config, permut):
    return [config[n] for n in permut]

def start_data_frame(df, result_file):
    df=df.loc[df.layer==0,:]
    df.to_csv(result_file, header=True)

def add_to_data_frame(df_add, result_file):
    df_add=df_add.loc[df_add.layer==0,:]
    with open(result_file, 'a') as f:
        df_add.to_csv(f, header=False)


def r_correction(x1,x2, r1,r2,p1,p2):
    (tree1,_),(tree2,_)=x1,x2
    n1=get_number_of_admixes(tree1)
    n2=get_number_of_admixes(tree2)
    cadmix_prior11= no_admixes(p=p1, admixes=n1, r=r1)
    cadmix_prior12 = no_admixes(p=p2, admixes=n1, r=r2)
    cadmix_prior21= no_admixes(p=p1, admixes=n2, r=r1)
    cadmix_prior22= no_admixes(p= p2, admixes=n2, r=r2)

    return cadmix_prior12-cadmix_prior11, cadmix_prior21-cadmix_prior22
    
def flipping(xs, posteriors, temperature_scheme, proposal_updates, rs=[],ps=[], posterior_function_list=[]):
    n=len(xs)
    step_permutation=list(range(n))
    count=0
    for _ in range(40):
        i,j = choice(n,2,False)
        post_i,post_j=posteriors[i],posteriors[j]
        temp_i,temp_j=temperature_scheme.get_temp(i), temperature_scheme.get_temp(j)
        logalpha=-(post_i[0]-post_j[0])*(1.0/temp_i-1.0/temp_j)
        if rs:
            i_correction, j_correction=r_correction(xs[i],xs[j], rs[i],rs[j],ps[i],ps[j])
            logalpha+=i_correction+j_correction
        else:
            i_correction, j_correction=0,0
        if logalpha>0 or random() < exp(logalpha):
            count+=1
            step_permutation[i], step_permutation[j]= step_permutation[j], step_permutation[i]
            posteriors[j],posteriors[i]=(post_i[0],post_i[1]+i_correction),(post_j[0], post_j[1]+j_correction)
            xs[i], xs[j] = xs[j], xs[i]

            proposal_updates[i], proposal_updates[j]=proposal_updates[j], proposal_updates[i]
    return xs, posteriors, step_permutation, proposal_updates

def _update_results(df_result, df_add):
    if df_result is None:
        df_result = df_add
    else:
        df_result = pd.concat([df_result, df_add])
    return df_result

##should match basic_chain_class.run_chain(start_tree, post, N, sample_verbose_scheme, overall_thinning, i_start_from, temperature, proposal_update)
def _pack_everything(xs, posteriors, temperature_scheme,printing_schemes,overall_thinnings,no_iterations,cum_iterations, proposal_updates=None, multiplier=None):
    return ([x, 
             posterior,
             no_iterations,
             printing_scheme,
             overall_thinnings,
             cum_iterations,
             temperature_scheme.get_temp(i),
             proposal_update,
             multiplier] for i,(x,posterior,printing_scheme,proposal_update) in enumerate(zip(xs,posteriors,printing_schemes,proposal_updates)))

def _unpack_everything(new_state, summaries, total_permutation):
    xs,posteriors, summs, proposal_updates = list(zip(*new_state))
    list_of_smaller_data_frames=[]
    for summ_data, n, i in zip(summs, total_permutation, list(range(len(total_permutation)))):
        iter_chain=chain((('iteration', summ_data[0]),),
                         ((summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,summ_data[1:])))
        df=pd.DataFrame.from_dict(dict(iter_chain))
        df['origin']=n
        df['layer']=i
        list_of_smaller_data_frames.append(df)
    df=pd.concat(list_of_smaller_data_frames)
    return list(xs), list(posteriors), df, list(proposal_updates)
    
