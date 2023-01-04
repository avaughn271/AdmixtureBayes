import msprime
from math import exp, log
demography = msprime.Demography()
outgroupsize = -(5000/2)/log(1-.08) #PLACE1
demography.add_population(name="outgroup", initial_size=outgroupsize)   #edit this.
PopulationsToSizes = [   #PLACE2
    ["A", -(2000/2)/log(1-.05)],    #FINAL
    ["B", -(2000/2)/log(1-.01) ],  #FINAL
    ["C", -(1000/2)/log(1-.04)],     #FINAL
    ["D",  -(1000/2)/log(1-.08)],  #FINAL
    ["E",   10000],  #FINAL
    ["pop1",  -(2000/2)/log(1-.015)],
    ["pop2", -(2000/2)/log(1-.01)],    #FINAL
    ["pop3", -(3000/2)/log(1-.03)],    #FINAL
    ["pop4", -(1000/2)/log(1-.005)  ],    #FINAL
    ["pop5", -(1000/2)/log(1-.002)  ]    #FINAL
]
####################################################################
for i in PopulationsToSizes:
    demography.add_population(name=i[0], initial_size=i[1])

AbsoluteTimes = []
InitialList = []

def getabsolutetime(val):
    for i in AbsoluteTimes:
        if val == i[0]:
            return(i[1])
    return 0.0

def addrelevantedges(absolutetime, child, ancestral):
    for i in PopulationsToSizes:
        if child == i[0]:
            timee = 1-exp(-(absolutetime - getabsolutetime(child))/(2*i[1]))
            break
    admixture = 1.0
    AbsoluteTimes.append([ancestral, absolutetime])
    return(child + " " + ancestral + " " + str(timee) + " "  + str(admixture))

def causeevent(timee, der1, der2, anc):
    demography.add_population_split(time=timee, derived=[der1, der2], ancestral=anc)
    InitialList.append(addrelevantedges(timee, der1, anc))
    InitialList.append(addrelevantedges(timee, der2, anc))

####################################################################

causeevent(1000, "pop5", "pop4", "B")
causeevent(2000, "pop1", "pop2", "A")
causeevent(3000, "B", "pop3", "C")
causeevent(4000, "C", "A", "D")

outgrouptime = 5000  # chnaged from 80000
demography.add_population_split(time=outgrouptime, derived=["D", "outgroup"], ancestral="E")

demography.sort_events()
ts = msprime.sim_ancestry(
    {"pop1": 5,
     "pop2": 5, 
     "pop3": 5,
     "pop4": 5,  "pop5": 5, 
     "outgroup": 5},
      demography=demography, sequence_length=20000000, recombination_rate = 1e-8)      #PLACE3 ENDS
######################################################################
outgroup = 1 - exp(-outgrouptime/(outgroupsize + outgroupsize)) +   1 - exp(-1000/(-(1000/2)/log(1-.08) + -(1000/2)/log(1-.08)))
mts = msprime.sim_mutations(ts, rate = 1e-8)

with open("TemporaryFiles/Data.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)

with open('TemporaryFiles/TrueTree.txt', 'w') as f:
    for line in InitialList:
        f.write(line + "\n")
        print(line)

with open('TemporaryFiles/TrueAdd.txt', 'w') as f:
    f.write(str(outgroup) + "\n")
print(outgroup)