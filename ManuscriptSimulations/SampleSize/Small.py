import msprime, tskit, numpy, random, sys, math

demography = msprime.Demography()
demography.add_population(name="pop1", initial_size=20000)
demography.add_population(name="pop2", initial_size=20000)
demography.add_population(name="pop3", initial_size=20000)
demography.add_population(name="outgroup", initial_size=20000)
demography.add_population(name="anc1", initial_size=20000)
demography.add_population(name="anc2", initial_size=20000)
demography.add_population(name="anc3", initial_size=20000)

demography.add_population_split(time=800, derived=["pop1", "pop2"], ancestral="anc1")
demography.add_population_split(time=1000, derived=["anc1", "pop3"], ancestral="anc2")
demography.add_population_split(time=1300, derived=["anc2", "outgroup"], ancestral="anc3")

ts = msprime.sim_ancestry({"pop1": 2, "pop2": 2, "pop3": 2,  "outgroup": 2}, demography=demography, sequence_length=3e7, recombination_rate = 1e-8)

mts = msprime.sim_mutations(ts, rate = 5e-8)

with open("Data.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)
