from pathos.multiprocessing import freeze_support
import pandas as pd
from MCMC import basic_chain_pool
from numpy.random import choice, random
from math import exp, log
from itertools import chain

def MCMCMC(starting_trees,    posterior_function, summaries, temperature_scheme,  printing_schemes, 
           iteration_scheme,  proposal_scheme, n_arg, verboseee,
           no_chains=None, numpy_seeds=None, multiplier= None, result_file=None):
    '''
    this function runs a MC3 using the basic_chain_unpacker. Let no_chains=number of chains. The inputs are
        starting_trees: a list of one or more trees that the chains should be started with
        proposal_scheme: a list of instances of classes that handles proposals and updates of said proposals. It has the basic functions
                                    - prop(x,pks): proposes the next tree(and returns proposal densities and statistics in pks)
                                    - adapt(mhr): updates the proposals based on the mhr ratio, mhr                                    
    '''
    df_result=None
    total_permutation=list(range(no_chains))
    xs = starting_trees

    freeze_support()
    
    #if numpy_sxeeds is None: #Sxeeddebug
    #    numpy_sxeeds=[None]*no_chains

    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme, numpy_seeds)
    posteriors = [posterior_function(x) for x in xs]

    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]
    
    cum_iterations=0
    totalMC3Acceptance = 0
    totalMC3FlipPropososals = 0
    chain0acceptance = 0
    chain0proposals = 0
    for no_iterations in iteration_scheme:
        if cum_iterations % 1000 == 0 and verboseee != "silent":
            print("Currently on iteration " +  str(cum_iterations) + " out of " + str(n_arg * 50))
        #letting each chain run for no_iterations:
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, printing_schemes, no_iterations, cum_iterations, proposal_updates, multiplier)
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
        xs, posteriors, permut, proposal_updates, totalMC3Acceptance, chain0acceptance, chain0proposals = flipping(xs, posteriors, temperature_scheme, proposal_updates, totalMC3Acceptance,chain0acceptance, chain0proposals)
        total_permutation=[total_permutation[n] for n in permut]
        totalMC3FlipPropososals = totalMC3FlipPropososals + 1
        if cum_iterations % 1000 == 0 and verboseee != "silent":
            print("Total MC3 Acceptance is: ", totalMC3Acceptance / totalMC3FlipPropososals, chain0acceptance/chain0proposals )
        cum_iterations+=no_iterations
    for chain in pool.group:
            chain.process.terminate()

def start_data_frame(df, result_file):
    df=df.loc[df.layer==0,:]
    df.to_csv(result_file, header=True)

def add_to_data_frame(df_add, result_file):
    df_add=df_add.loc[df_add.layer==0,:]
    with open(result_file, 'a') as f:
        df_add.to_csv(f, header=False)
    
def flipping(xs, posteriors, temperature_scheme, proposal_updates, totalMC3Acceptance,chain0acceptance, chain0proposals):
    n=len(xs)
    step_permutation=list(range(n))
    count=0
    for _ in range(1):
        i,j = choice(n,2,False)
        if (i == 0 or j == 0):
            chain0proposals = chain0proposals + 1
        post_i,post_j=posteriors[i],posteriors[j]
        temp_i,temp_j=temperature_scheme[i], temperature_scheme[j]
        logalpha = (1.0/temp_j - 1.0/temp_i) * (post_i[0] + post_i[1])  + (1.0/temp_i - 1.0/temp_j) * (post_j[0] + post_j[1])
        if log(random()) < logalpha:
            count+=1
            step_permutation[i], step_permutation[j]= step_permutation[j], step_permutation[i]
            posteriors[j],posteriors[i]=(post_i[0],post_i[1]),(post_j[0], post_j[1])
            xs[i], xs[j] = xs[j], xs[i]
            if (i == 0 or j == 0):
                chain0acceptance = chain0acceptance + 1
                print("switching 0")

            proposal_updates[i], proposal_updates[j]=proposal_updates[j], proposal_updates[i]
    totalMC3Acceptancenew = totalMC3Acceptance + count
    return xs, posteriors, step_permutation, proposal_updates, totalMC3Acceptancenew, chain0acceptance, chain0proposals

def _update_results(df_result, df_add):
    if df_result is None:
        df_result = df_add
    else:
        df_result = pd.concat([df_result, df_add])
    return df_result

def _pack_everything(xs, posteriors, temperature_scheme,printing_schemes,no_iterations,cum_iterations, proposal_updates=None, multiplier=None):
    return ([x, posterior,no_iterations, printing_scheme, 40, cum_iterations, temperature_scheme[i], proposal_update, multiplier] for i,(x,posterior,printing_scheme,proposal_update) in enumerate(zip(xs,posteriors,printing_schemes,proposal_updates)))

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
