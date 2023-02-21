import msprime
from math import exp, log
demography = msprime.Demography()
outgroupsize =  -(6000/2)/log(1-.080)
demography.add_population(name="outgroup", initial_size=outgroupsize)   #edit this.
PopulationsToSizes = [   #PLACE2
    ["A1", -(2000/2)/log(1-.015)],    #FINAL
    ["A2",  -(1000/2)/log(1-.030)],   #FINAL
    ["B", -(1000/2)/log(1-.055)],  #FINAL
    ["C", -(1000/2)/log(1-.07)],     #FINAL
    ["D", -(3000/2)/log(1-.06)],  #FINAL
    ["E",  -(1000/2)/log(1-.08)],  #FINAL
    ["F",  10000],     #FINAL
    ["pop1", -(4000/2)/log(1-.025)],
    ["pop2", -(3000/2)/log(1-.01)],    #FINAL
    ["pop3", -(1000/2)/log(1-.005)],    #FINAL
    ["pop4", -(2000/2)/log(1-.02)]    #FINAL
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
demography.add_admixture(time=1000, derived="pop3", ancestral=["A1", "A2"], proportions=[0.4, 0.6])       #PLACE3

InitialList.append(addadmixture(1000, "pop3",  "A1", "A2", 0.4, 0.6))

causeevent(2000, "A2", "pop4", "D")
causeevent(3000, "A1", "pop2", "B")
causeevent(4000, "pop1", "B", "C")
causeevent(5000, "C", "D", "E")

demography.add_population_split(time=6000, derived=["E", "outgroup"], ancestral="F")

demography.sort_events()
ts = msprime.sim_ancestry(
    {"pop1": 15,"pop2": 15, "pop3": 15, "pop4": 15,  "outgroup": 15}, 
      demography=demography, sequence_length=200000000, recombination_rate = 1e-8)      #PLACE3 ENDS
######################################################################
outgroup = 1 - exp(-6000/(outgroupsize + outgroupsize)) + 1 - exp(-1000/(-(1000/2)/log(1-.139)-(1000/2)/log(1-.139) ))
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