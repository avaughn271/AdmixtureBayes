import msprime
from math import exp, log
demography = msprime.Demography()
outgroupsize =    -(8000/2)/log(1-.184)
demography.add_population(name="outgroup", initial_size=outgroupsize)   #edit this.
PopulationsToSizes = [   #PLACE2
    ["A1",   -(1500/2)/log(1-.003)],    #FINAL
    ["A2",    -(3500/2)/log(1-.018)],   #FINAL
    ["B",   -(1000/2)/log(1-.044)],  #FINAL
    ["C1",   -(1499/2)/log(1-.013)],     #FINAL
    ["C2",   -(499/2)/log(1-.002)],       #FINAL
    ["D",    -(3000/2)/log(1-.018)],  #FINAL
    ["E",   -(2000/2)/log(1-.032)],  #FINAL
    ["F",    -(2000/2)/log(1-.041)],     #FINAL
    ["G", 5000],   
    ["H",    -(1000/2)/log(1-.110)],       
    ["I",    -(1000/2)/log(1-.044)],
    ["J",   -(2000/2)/log(1-.028)], 
    ["pop1",    -(5000/2)/log(1-.082)],
    ["pop2",  -(3000/2)/log(1-.097)],    #FINAL
    ["pop3",  -(500/2)/log(1-.019)  ],    #FINAL
    ["pop4",  -(501/2)/log(1-.016) ],    #FINAL
    ["pop5",  -(1000/2)/log(1-.030) ],    #FINAL
    ["pop6",    -(6000/2)/log(1-.039)]
]

for i in PopulationsToSizes: print(i)
####################################################################
for i in PopulationsToSizes:
    demography.add_population(name=i[0], initial_size=i[1])

DivPairs = []
AbsoluteTimes = []
AdmixProp = []
InitialList = []

def getalphabetical(val):
    for i in DivPairs:
        if val in i and i[0] <  i[1]:
            return(i[0])
        if val in i and i[1] <  i[0]:
            return(i[1])
    return(val)

def getadmixtureproportions(val):
    for i in AdmixProp:
        if val == i[0]:
            return(i[1])
    return 1.0

def getabsolutetime(val):
    for i in AbsoluteTimes:
        if val == i[0]:
            return(i[1])
    return 0.0

def addrelevantedges(absolutetime, child, ancestral):
    childreturn = getalphabetical(child)
    for i in PopulationsToSizes:
        if child == i[0]:
            timee = 1-exp(-(absolutetime - getabsolutetime(childreturn))/(2*i[1]))
            break
    admixture = getadmixtureproportions(child)
    AbsoluteTimes.append([ancestral, absolutetime])
    return(childreturn + " " + ancestral + " " + str(timee) + " "  + str(admixture))

def addadmixture(absolutetime, child, ancestral1, ancestral2, prop1, prop2):
    childreturn = ancestral1
    if ancestral2 < ancestral1:
        childreturn = ancestral2
    for i in PopulationsToSizes:
        if child == i[0]:
            timee = 1-exp(-(absolutetime - getabsolutetime(child))/(2*i[1]))
            break
    admixture = getadmixtureproportions(child)
    AbsoluteTimes.append([childreturn, absolutetime])
    AdmixProp.append([ancestral1, prop1])
    AdmixProp.append([ancestral2, prop2])
    DivPairs.append([ancestral1, ancestral2])
    return(child + " " + childreturn + " " + str(timee) + " "  + str(admixture))

def causeevent(timee, der1, der2, anc):
    demography.add_population_split(time=timee, derived=[der1, der2], ancestral=anc)
    InitialList.append(addrelevantedges(timee, der1, anc))
    InitialList.append(addrelevantedges(timee, der2, anc))

####################################################################
demography.add_admixture(time=500, derived="pop3", ancestral=["A1", "A2"], proportions=[0.39, 0.61])       #PLACE3
demography.add_admixture(time=501, derived="pop4", ancestral=["C1", "C2"], proportions=[0.6, 0.4])       #PLACE3

InitialList.append(addadmixture(500, "pop3",  "A1", "A2", 0.39, 0.61))
InitialList.append(addadmixture(501, "pop4",  "C1", "C2", 0.6, 0.4))

causeevent(1000, "C2", "pop5", "D")
causeevent(2000, "A1", "C1", "B")
causeevent(3000, "pop2", "B", "E")
causeevent(4000, "A2", "D", "J")
causeevent(5000, "pop1", "E", "F")
causeevent(6000, "J", "pop6", "I")
causeevent(7000, "I", "F", "H")

demography.add_population_split(time=8000, derived=["H", "outgroup"], ancestral="G")

demography.sort_events()
ts = msprime.sim_ancestry(
    {"pop1": 5,"pop2": 5, "pop3": 5, "pop4": 5, "pop5": 5, "pop6": 5,"outgroup": 5},  
      demography=demography, sequence_length=80000000, recombination_rate = 1e-8)      #PLACE3 ENDS
######################################################################
outgroup = 1 - exp(-8000/(outgroupsize + outgroupsize)) + 1 - exp(-1000/( -(1000/2)/log(1-.110) -(1000/2)/log(1-.110)))
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