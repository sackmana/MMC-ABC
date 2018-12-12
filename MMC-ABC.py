import random
import sys
import math
import itertools
import scipy.stats
import os
import numpy
from multiprocessing import Pool
from multiprocessing import cpu_count

#########################################################################################################################################################
#########################################################################################################################################################
######### Please carefully read the MMC-ABC manual for information on how to use this program, as well as information on the proper formating of data ###  
######### files. Depending on the commands used, MMC-ABC is capable of inferring N, psi, and site-specific s and h from time-series polymorphism ########
######### data.##########################################################################################################################################
#########################################################################################################################################################



s_min=-0.2
s_max=0.6
numsims=10000
recomb=1e-8
psi_min=0
psi_max=0.3
numthreads=cpu_count()
n_min=500
n_max=10000
n_apriori=False
psi_apriori=False
output_file_number=1
Part2Replicates=1000
het=0.5
ploidy=2
h_min=0
h_max=1

#########################################################################################################################################################
#########################################################################################################################################################
######### Here we read in the command line arguments used. For a detailed description of available arguments, please see the manual. ####################  
#########################################################################################################################################################
#########################################################################################################################################################

list_of_possible_args=["--h_yes","--ploidy","--h","--s_min","--part2sims","--s_max","--part2N","--part2psi","--part1sims","--numthreads","--n_min","--n_max","--n_apriori","--recomb","--psi_min","--psi_max","--output_file_number","--psi_apriori","--h_min","--h_max"]
for thing in sys.argv[1:-1:2]:
	if str(thing) not in list_of_possible_args:
		print "Invalid variable entry"
		sys.exit()
		
if("--s_min" in sys.argv):
	s_min=float(sys.argv[sys.argv.index("--s_min")+1])
if("--s_max" in sys.argv):
	s_max=float(sys.argv[sys.argv.index("--s_max")+1])
if("--part1sims" in sys.argv):
	numsims=int(sys.argv[sys.argv.index("--part1sims")+1])
if("--numthreads" in sys.argv):
	numthreads=int(sys.argv[sys.argv.index("--numthreads")+1])
if("--n_max" in sys.argv):
	n_max=int(sys.argv[sys.argv.index("--n_max")+1])
if("--n_min" in sys.argv):
	n_min=int(sys.argv[sys.argv.index("--n_min")+1])
if("--n_apriori" in sys.argv):
	n_apriori=int(sys.argv[sys.argv.index("--n_apriori")+1])
if("--recomb" in sys.argv):
	recomb=sys.argv[sys.argv.index("--recomb")+1]
if("--psi_min" in sys.argv):
	psi_min=float(sys.argv[sys.argv.index("--psi_min")+1])
if("--psi_max" in sys.argv):
	psi_max=float(sys.argv[sys.argv.index("--psi_max")+1])
if("--psi_apriori" in sys.argv):
	psi_apriori=float(sys.argv[sys.argv.index("--psi_apriori")+1])
if("--output_file_number" in sys.argv):
	output_file_number=sys.argv[sys.argv.index("--output_file_number")+1]
if("--part2sims" in sys.argv):
	Part2Replicates=int(sys.argv[sys.argv.index("--part2sims")+1])
if("--ploidy" in sys.argv):
	ploidy=int(sys.argv[sys.argv.index("--ploidy")+1])
if("--h" in sys.argv):
	het=int(sys.argv[sys.argv.index("--h")+1])
if("--h_min" in sys.argv):
	h_min=float(sys.argv[sys.argv.index("--h_min")+1])
if("--h_max" in sys.argv):
	h_max=float(sys.argv[sys.argv.index("--h_max")+1])

print "Beginning Part 1 of MMC-ABC: Estimation of N and Psi"

######### If the user specifies a haploid population, we set h=1 for the slim simulations. #########
if ploidy==1:
	het=1

filenumber=sys.argv[-1]

#########################################################################################################################################################
#########################################################################################################################################################
######### Here we read in the input WF_ABC formatted file, and extract the frequencies at each sampled timepoint.  Note that currently this program only  
######### allows for an equal sample size at each timepoint of all sites. ###############################################################################
#########################################################################################################################################################

input=open(filenumber,"r")
text=input.readlines()
input.close()
numsites=int(text[0].split()[0])
timepoints=int(text[0].split()[1])
gens=[int(s) for s in text[1].split(',')[0:-1]]

samplesizes=[]
counts=[]
freqs=[]
for x in range(numsites):
	freqs.append([])
	samplesizes.append([])
	counts.append([])
	for y in range(timepoints):
		samplesizes[x].append(float(text[(2+x*2)].split(',')[y]))
		counts[x].append(float(text[(3+x*2)].split(',')[y]))
		freqs[x].append(float(text[(3+x*2)].split(',')[y])/float(text[(2+x*2)].split(',')[y]))

######### We create empty vectors to store values of F and variance in frequency changes. #########

Fvalues=[]

######### Now we measure F across each site in the input file. #########

for x in range(numsites):
	Fvalues.append([])
	
	######### We only measure frequency changes for time point pairs that have frequencies above 0.01 and below 0.99 on both ends. #########
	
	for y in range(timepoints-1):
		a=freqs[x][y]
		b=freqs[x][y+1]
		if a>0.01 and b>0.01:
			if a<0.99 and b<0.99:
				sample1=float(samplesizes[x][y])
				sample2=float(samplesizes[x][y+1])
				time1=float(gens[y])
				time2=float(gens[y+1])
				c=(a+b)/2
				Fs=((a-b)**2)/(c*(1-c))
				n=scipy.stats.hmean([sample1,sample2])
				Fsprime=(Fs*(1-(1/(2*n)))-(2/n))/((1+(Fs/4))*(1-(1/sample2))*(time2-time1))
				Fvalues[x].append(Fsprime)


######### Now we measure mean F across all sites. #########
		
total=0.0
count=0.0
for thing in Fvalues:
	for thingy in thing:
		total+=thingy
		count+=1
if ploidy==2:
	#print 1/(total/count)
	EstimatedSampleNe=1/(2*(total/count))
elif ploidy==1:
	#print 1/(total/count)
	EstimatedSampleNe=1/(total/count)


################################################################################################################################################
################################################################################################################################################
######### Now we set up the parameters for the ABC simulations. ################################################################################
################################################################################################################################################
################################################################################################################################################


gensrescaled=[s-gens[0] for s in gens]
firstgens=[]
startfreqers1=[]
full_trajectories=[]
for x in range(len(freqs)):
	countbit=0
	for y in range(len(freqs[x])):
		if freqs[x][y]!=0 and countbit==0:
			firstgens.append(gensrescaled[y])
			countbit+=1
			startfreqers1.append(freqs[x][y])
			full_trajectories.append(counts[x][y:len(counts[x])])
startfreqers2=startfreqers1
full_trajectories2=full_trajectories

firstgens2=firstgens


ziplist_startfreqs=zip(firstgens,startfreqers1,full_trajectories)
ziplist_startfreqs.sort()
firstgens,startfreqers1,full_trajectories=zip(*ziplist_startfreqs)

numberpergen=[]
for x in range(len(gensrescaled)-1):
	numberpergen.append(firstgens.count(gensrescaled[x]))
	
numberpergeninput='c('
for thing in numberpergen[0:-1]:
	numberpergeninput+=str(thing)+','
numberpergeninput+=str(numberpergen[-1])+')'

gensrescaled=gensrescaled[1:len(gensrescaled)]
outputtimes='c('
for thing in range(len(gensrescaled)-1):
	outputtimes+=str(gensrescaled[thing])+','
outputtimes+=str(gensrescaled[-1])+')'


#########################################################################################################################################################
############# We create vectors of N and psi drawn randomly from their priors to use for the ABC simulations ############################################
############# If either parameter was specified instead of using a prior, that is handled here as well. #################################################  
#########################################################################################################################################################
#########################################################################################################################################################

SimNs=[]
SimPsis=[]
startfreqers=[]
startfreqs=[]
if n_apriori!=False:
	for x in range(numsims):
		SimNs.append(n_apriori)
		NE=str(n_apriori)
		startfreqers.append([])
		if ploidy==2:
			startfreqers[x]=[round(int(NE)*2*s) for s in startfreqers1 if s!=0]
		elif ploidy==1:
			startfreqers[x]=[round(int(NE)*s) for s in startfreqers1 if s!=0]
		startfreqs.append('')
		startfreqs[x]='c('
		for thing in range(len(startfreqers[x])-1):
			startfreqs[x]+=str(int(startfreqers[x][thing]))+','
		startfreqs[x]+=str(int(startfreqers[x][-1]))+')'
elif n_apriori==False:
	for x in range(numsims):
		NE=str(random.randint(n_min,n_max))
		SimNs.append(int(NE))
		startfreqers.append([])
		if ploidy==2:
			startfreqers[x]=[round(int(NE)*2*s) for s in startfreqers1 if s!=0]
		elif ploidy==1:
			startfreqers[x]=[round(int(NE)*s) for s in startfreqers1 if s!=0]
		startfreqs.append('')
		startfreqs[x]='c('
		for thing in range(len(startfreqers[x])-1):
			startfreqs[x]+=str(int(startfreqers[x][thing]))+','
		startfreqs[x]+=str(int(startfreqers[x][-1]))+')'
if psi_apriori!=False:
	for x in range(numsims):
		SimPsis.append(psi_apriori)
elif psi_apriori==False:
	for x in range(numsims):
		SimPsis.append(round(numpy.random.uniform(psi_min,psi_max),4))
SimNeEstimates=[]
for x in range(numsims):
	SimNeEstimates.append(0)
abcsim1range=range(numsims)
choicevector=[0,1]

######### We define the simulation function #########

def abcsim1(x):
		######### Each simulation run will generate a temporary file name 'output_x.txt' where x=simulation #. The file will be read, analyzed, and then deleted. ########
		#print x
		if ploidy==2:
			os.system('slim -d N='+str(SimNs[x])+' -d psi='+str(SimPsis[x])+' -d startfreqs="'+startfreqs[x]+'" -d s=0 -d numberpergeninput="'+numberpergeninput+'" -d numsites='+str(len(startfreqers[x]))+' -d recomb='+str(recomb)+' -d outputtimes="'+outputtimes+'" part1_diploid.txt > output'+str(x)+'.txt')
		elif ploidy==1:
			os.system('slim -d N='+str(SimNs[x])+' -d psi='+str(SimPsis[x])+' -d startfreqs="'+startfreqs[x]+'" -d s=0 -d numberpergeninput="'+numberpergeninput+'" -d numsites='+str(len(startfreqers[x]))+' -d recomb='+str(recomb)+' -d outputtimes="'+outputtimes+'" part1_haploid.txt > output'+str(x)+'.txt')

		######### We analyze the simulation data, estimating F for each simulation. #########

		with open('output'+str(x)+'.txt','r') as siminputpart1:
			timeliners=list(siminputpart1)[(-1*(timepoints-1)):]
		timelines=[]
		######### timeliners now contains a list with the text describing allele counts at each timepoint #########
		######### Below we will convert this into a list of counts and a list of frequencies stored in simfreqtimes #########	
	
		for thing in timeliners:
			timelineline=thing.split(',')
			del(timelineline[-1])
			timelines.append([float(s) for s in timelineline])
		simfreqtimes=[]

		######### Here we are attempting to replicate the filtering criteria used to generate the sampled data wherein we remove all sites at #########
		######### which we have fewer than three informative timepoints and sites for which all time points are below 2.5% #########
		indices=[]
		for h in range(len(timelines[0])):
			countingtimes=0
			for j in range(timepoints-1):
				if ploidy==2:
					if timelines[j][h] != SimNs[x]*2 and timelines[j][h]!=0:
						countingtimes+=1
				elif ploidy==1:
					if timelines[j][h] != SimNs[x] and timelines[j][h]!=0:
						countingtimes+=1
			if countingtimes<2:
				indices.append(h)
		indices.sort(reverse=True)
		startfreqsim1=startfreqers[x]
		for thing in indices:
			for j in range(timepoints-1):
				del timelines[j][thing]
			del startfreqsim1[thing]

		########## We down-sample the data to match the sample sizes from our sample data ##########

		for thing in range(len(startfreqsim1)):
			simfreqtimes.append([])
			for quabber in range(timepoints):
				simfreqtimes[thing].append(0.0)
		for thing in range(len(timelines)):
			for thingy in range(len(timelines[thing])):
				simfreqtimes[thingy][thing+1]=(timelines[thing][thingy])
		for thing in range(len(simfreqtimes)):
			simfreqtimes[thing][([0]+gensrescaled).index(firstgens[thing])]=(startfreqsim1[thing])
			for thingy in range(timepoints):
				if ploidy==2:
					simfreqtimes[thing][thingy]=sum(numpy.random.choice(choicevector,int(samplesizes[thing][thingy]),p=[(((2*SimNs[x])-float(simfreqtimes[thing][thingy]))/(2*SimNs[x])),(float(simfreqtimes[thing][thingy])/(2*SimNs[x]))]))/float(samplesizes[thing][thingy]) #Note that this is another place where uneven sample sizes do not work.
				elif ploidy==1:
					simfreqtimes[thing][thingy]=sum(numpy.random.choice(choicevector,int(samplesizes[thing][thingy]),p=[(((SimNs[x])-float(simfreqtimes[thing][thingy]))/(SimNs[x])),(float(simfreqtimes[thing][thingy])/(SimNs[x]))]))/float(samplesizes[thing][thingy]) #Note that this is another place where uneven sample sizes do not work.

		########## We calculate Ne across all simulated sites ##########
	
		SimFvalues=[]
		outputtimes2=[0]+gensrescaled
		for z in range(len(simfreqtimes)):
			SimFvalues.append([])
			for y in range(timepoints-1):
				sample1=float(samplesizes[z][y])
				sample2=float(samplesizes[z][y+1])
				a=simfreqtimes[z][y]
				b=simfreqtimes[z][y+1]
				if a>0.01 and b>0.01:
					if a<0.99 and b<0.99:
						if ((a+b)/2)!=0 and ((a+b)/2)!=1:
							time1=float(outputtimes2[y])
							time2=float(outputtimes2[y+1])
							c=(a+b)/2
							Fs=((a-b)**2)/(c*(1-c))
							n=scipy.stats.hmean([sample1,sample2]) #Note yet another place where uneven sample sizes do not work.
							Fsprime=0
							if ((1+(Fs/4))*(1-(1/n))*(time2-time1))!=0:
								Fsprime=(Fs*(1-(1/(2*n)))-(2/n))/((1+(Fs/4))*(1-(1/n))*(time2-time1))
							SimFvalues[z].append(Fsprime)
		
		total=0.0
		count=0.0
		for thing in SimFvalues:
			for thingy in thing:
				total+=thingy
				count+=1
		SimEstimatedSampleNe=1000000
		if ploidy==2:
			if count!=0:
				SimEstimatedSampleNe=1/(2*(total/count))
		elif ploidy==1:
			if count!=0:
				SimEstimatedSampleNe=1/(total/count)
		
		########## We delete the temporary data file ##########
		
		os.system('rm output'+str(x)+'.txt')

		return(x,SimEstimatedSampleNe)
		os.system('rm output'+str(x)+'.txt')


######### We set up multithreading of the ABC simulations #########

p=Pool(numthreads)
valuelists=p.map(abcsim1,abcsim1range)
p.close()
p.join()
valuelists.sort()
for x in abcsim1range:
	SimNeEstimates[x]=float(valuelists[x][1])

NeEstimateDiffs=[abs(s-EstimatedSampleNe) for s in SimNeEstimates]
ziplist=zip(NeEstimateDiffs,SimNs,SimPsis)
ziplist.sort()
NeEstimateDiffs,SimNs1,SimPsis1 =zip(*ziplist)

print "Mean of the posterior of N: "+str(sum(SimNs1[0:(int((round(numsims*0.01))))])/float(round(numsims*0.01)))
print "Mean of the posterior of Psi: " +str((sum(SimPsis1[0:(int((round(numsims*0.01))))])/float(round(numsims*0.01))))

########## We output the joint posteriors of N and Psi to a file named 'N_psi_posterior.csv' ##########

Npsi_posterior=open(filenumber.split('.')[0]+'_N_psi_joint_posterior.csv','w')
Npsi_posterior.write('N,Psi,\n')
for x in range(int(round(numsims*0.01))):
	Npsi_posterior.write(str(SimNs1[x])+','+str(SimPsis1[x])+'\n')
Npsi_posterior.close()

################################################################################################################################################
############ Now we beging Part 2 of MMC-ABC ###################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

print "Beginning part two of MMC-ABC: Estimation of site-specific selection coefficients using the joint posterior of N and Psi."


######### We generate vectors of N and Psi to be used for the ABC simulations, drawn randomly and jointly from their joint posterior. ##########
######### If values of N or Psi were specified instead of a prior for part 1, those values are used here. ##########

joint_posterior=zip(SimNs1[0:(int((round(numsims*0.01))))],SimPsis1[0:(int((round(numsims*0.01))))])
randomlysampled=numpy.random.randint(0,int(round(numsims*0.01)),Part2Replicates)
part2_N_psi_vector=[joint_posterior[s] for s in randomlysampled]
Part2N=[]
for x in range(Part2Replicates):
	Part2N.append(part2_N_psi_vector[x][0])
if("--n_apriori" in sys.argv):
	Part2N=[int(sys.argv[sys.argv.index("--n_apriori")+1])]*Part2Replicates
Part2Psi=[]
for x in range(Part2Replicates):
	Part2Psi.append(part2_N_psi_vector[x][1])
if("--psi_apriori" in sys.argv):
	Part2Psi=[float(sys.argv[sys.argv.index("--psi_apriori")+1])]*Part2Replicates


########## We generate vectors of values to be used for the initial counts of each mutation, calculated from their observed frequencies ##########
########## in the data and the population size used (adjusted for haploid and diploid populations) ##########
startertracker=[]

Part2Startfreqs=[]
for y in range(Part2Replicates):
	if ploidy==2:
		Part2Startfreqs.append([int(s*Part2N[y]*2) for s in startfreqers2])
	elif ploidy==1:
		Part2Startfreqs.append([int(s*Part2N[y]) for s in startfreqers2])
Part2Svector=[]
Part2HetVector=[]
Part2OutputTimes=[]
for x in range(len(startfreqers2)):
	startertracker.append(([0]+gensrescaled).index(firstgens2[x]))
	outputtimes='c('
	for thing in range(([0]+gensrescaled).index(firstgens2[x]),len(gensrescaled)):
		outputtimes+=str(([0]+gensrescaled)[thing])+','
	outputtimes+=str(gensrescaled[-1])+')'
	Part2OutputTimes.append(outputtimes)
	Part2Svector.append([])
	Part2HetVector.append([])
	for y in range(Part2Replicates):
		Part2Svector[x].append(round(numpy.random.uniform(s_min,s_max),2))
		Part2HetVector[x].append(round(numpy.random.uniform(h_min,h_max),2))


########## We define the simulation function, which simulates a population of size N and parameter Psi chosen from their joint posterior ###########
########## with initial frequency at each site following the observed data. The simulations output a temporary file for each replicate ##########
########## at each site, containing output from SLiM and an output simulated trajectory for the mutation. Those trajectories are read and ##########
########## used to calculate the distance between the summary statistics for the data and the simulated trajectory. ##########


def abcsim2(x):
	#if x in range(0,Part2Replicates,100):
		#print x
	if ploidy==2:
		if ("--h_yes" in sys.argv):
			os.system('slim -d N='+str(Part2N[x])+' -d psi='+str(Part2Psi[x])+' -d h='+str(Part2HetVector[thing][x])+' -d startingfrequency='+str(Part2Startfreqs[x][thing])+' -d s=' + str(Part2Svector[thing][x])+' -d outputtimesinit="'+Part2OutputTimes[thing]+'" part2_diploid.txt > output'+str(x)+'.txt')
		else:
			os.system('slim -d N='+str(Part2N[x])+' -d psi='+str(Part2Psi[x])+' -d startingfrequency='+str(Part2Startfreqs[x][thing])+' -d s=' + str(Part2Svector[thing][x])+' -d h='+str(het)+' -d outputtimesinit="'+Part2OutputTimes[thing]+'" part2_diploid.txt > output'+str(x)+'.txt')
	elif ploidy==1:
		os.system('slim -d N='+str(Part2N[x])+' -d psi='+str(Part2Psi[x])+' -d startingfrequency='+str(Part2Startfreqs[x][thing])+' -d s=' + str(Part2Svector[thing][x])+' -d h='+str(het)+' -d outputtimesinit="'+Part2OutputTimes[thing]+'" part2_haploid.txt > output'+str(x)+'.txt')
	with open('output'+str(x)+'.txt','r') as filename2:
		temptraj=filename2.readlines()[-1]
	os.system('rm output'+str(x)+'.txt')
	SimFivalues2=[]
	SimFdvalues2=[]		
	simulated_trajectory=[int(s) for s in temptraj.split(',')]
	if ploidy==2:
		for alpha in range(len(simulated_trajectory)):
			simulated_trajectory[alpha]=sum(numpy.random.choice([0,1],int(samplesizes[thing][alpha+1+startertracker[thing]]),p=[((2*Part2N[x])-float(simulated_trajectory[alpha]))/(2*Part2N[x]),float(simulated_trajectory[alpha])/(2*Part2N[x])]))
		simulated_trajectory=[Part2Startfreqs[x][thing]/((2*Part2N[x])/samplesizes[thing][startertracker[thing]])]+simulated_trajectory
	elif ploidy==1:
		for alpha in range(len(simulated_trajectory)):
			simulated_trajectory[alpha]=sum(numpy.random.choice([0,1],int(samplesizes[thing][alpha+1+startertracker[thing]]),p=[((Part2N[x])-float(simulated_trajectory[alpha]))/(Part2N[x]),float(simulated_trajectory[alpha])/(Part2N[x])]))
		simulated_trajectory=[Part2Startfreqs[x][thing]/((Part2N[x])/samplesizes[thing][startertracker[thing]])]+simulated_trajectory
	for z in range(len(simulated_trajectory)-1):
		sample1=float(samplesizes[thing][z+startertracker[thing]])
		sample2=float(samplesizes[thing][z+startertracker[thing]+1])
		a=simulated_trajectory[z]/float(samplesizes[thing][z+startertracker[thing]])
		b=simulated_trajectory[z+1]/float(samplesizes[thing][z+1+startertracker[thing]])
		if a>0.01 or b>0.01:
			if a<0.99 or b<0.99:
				if ((a+b)/2)!=0 and ((a+b)/2)!=1:
					time1=float(outputtimessim2[z])
					time2=float(outputtimessim2[z+1])
					c=(a+b)/2
					Fs=((a-b)**2)/(c*(1-c))
					n=scipy.stats.hmean([sample1,sample2])
					Fsprime=0
					if ((1+(Fs/4))*(1-(1/n))*(time2-time1))!=0:
						Fsprime=(Fs*(1-(1/(2*n)))-(2/n))/((1+(Fs/4))*(1-(1/n))*(time2-time1))
					if a<b:
						SimFivalues2.append(Fsprime)
					elif a>b:
						SimFdvalues2.append(Fsprime)

	if SimFdvalues2!=[] and SimFivalues2!=[]:
		returnvalue=(numpy.sqrt((numpy.mean(SimFivalues2)-InputMeanFi2)**2+(numpy.mean(SimFdvalues2)-InputMeanFd2)**2))
	elif SimFivalues2==[] and SimFdvalues2==[]:
		returnvalue=(numpy.sqrt(InputMeanFi2**2+InputMeanFd2**2))
	elif SimFivalues2==[]:
		returnvalue=(numpy.sqrt((InputMeanFi2)**2+(numpy.mean(SimFdvalues2)-InputMeanFd2)**2))
	else:
		returnvalue=(numpy.sqrt((numpy.mean(SimFivalues2)-InputMeanFi2)**2+(InputMeanFd2)**2))
	if len(SimFivalues2+SimFdvalues2)<=1:
		returnvalue=100000000
	return(x,returnvalue)
	
gensrescaled2=[0]+gensrescaled
valuecums2=[]
valuecumsf2=[]
valuecums2het=[]


########## We measure Fsi and Fsd for each site and run the part 2 simulations, generating posteriors for each site. ##########

abcsim2range=range(len(Part2Startfreqs[0]))
for thing in abcsim2range:
	input_trajectory=full_trajectories2[thing]
	SimFValues3=[]
	InputFivalues2=[]
	InputFdvalues2=[]
	outputtimessim2=gensrescaled2[gensrescaled2.index(firstgens2[thing]):len(gensrescaled2)]
	for z in range(len(input_trajectory)-1):
		sample1=float(samplesizes[thing][z+startertracker[thing]])
		sample2=float(samplesizes[thing][z+startertracker[thing]+1])
		a=input_trajectory[z]/float(samplesizes[thing][z+startertracker[thing]])
		b=input_trajectory[z+1]/float(samplesizes[thing][z+1+startertracker[thing]])
		if a>0.01 or b>0.01:
			if a<0.99 or b<0.99:
				if ((a+b)/2)!=0 and ((a+b)/2)!=1:
					time1=float(outputtimessim2[z])
					time2=float(outputtimessim2[z+1])
					c=(a+b)/2
					Fs=((a-b)**2)/(c*(1-c))
					n=scipy.stats.hmean([sample1,sample2])
					Fsprime=0
					if ((1+(Fs/4))*(1-(1/n))*(time2-time1))!=0:
						Fsprime=(Fs*(1-(1/(2*n)))-(2/n))/((1+(Fs/4))*(1-(1/n))*(time2-time1))
					if a<b:
						InputFivalues2.append(Fsprime)
					elif a>b:
						InputFdvalues2.append(Fsprime)
	if InputFivalues2!=[]:
		InputMeanFi2=numpy.mean(InputFivalues2)
	else:
		InputMeanFi2=0
	if InputFdvalues2!=[]:
		InputMeanFd2=numpy.mean(InputFdvalues2)
	else:
		InputMeanFd2=0

	print "Site "+str(thing+1)
	abcsim2range2=range(Part2Replicates)
	p=Pool(numthreads)
	valuelists2=p.map(abcsim2,abcsim2range2)
	p.close()
	SimFValues3=[]
	for thingery in valuelists2:
		SimFValues3.append(thingery[1])
	s_output_file=open('s_outputs_'+filenumber,'a')
	########## We sort each simulation by distance and generate our posterior. ###########	
	if ("--h_yes" in sys.argv):
		ziplist_sim6=zip(SimFValues3,Part2Svector[thing],Part2HetVector[thing])
		ziplist_sim6.sort(key=lambda x: x[0])
		SimFValues3_sorted,SimFValues3_s_values,SimFValues3_het_values=zip(*ziplist_sim6)
		valuecumsf2.append(SimFValues3_sorted[0:int(Part2Replicates)])
		valuecums2.append(SimFValues3_s_values[0:int(Part2Replicates)])
		valuecums2het.append(SimFValues3_het_values[0:int(Part2Replicates)])
		s_output_file.write(str(numpy.mean(SimFValues3_s_values[0:int(Part2Replicates*0.01)]))+','+str(numpy.mean(SimFValues3_het_values[0:int(Part2Replicates*0.01)]))+'\n')
	else:
		ziplist_sim6=zip(SimFValues3,Part2Svector[thing])
		ziplist_sim6.sort(key=lambda x: x[0])
		SimFValues3_sorted,SimFValues3_s_values=zip(*ziplist_sim6)
		valuecumsf2.append(SimFValues3_sorted[0:int(Part2Replicates)])
		valuecums2.append(SimFValues3_s_values[0:int(Part2Replicates)])
		s_output_file.write(str(numpy.mean(SimFValues3_s_values[0:int(Part2Replicates*0.01)]))+'\n')
	s_output_file.close()
########## Now we output the posterior for each site to a csv file, with each column heading indicating the site number. ##########

s_posteriors=open(filenumber.split('.')[0]+'_s_posteriors.csv','w')
if ("--h_yes" in sys.argv):
	for x in range(len(startfreqers2)):
		s_posteriors.write('Site_'+str(x+1)+'_s,Site_'+str(x+1)+'_h,')
	s_posteriors.write('\n')
	for x in range(int(Part2Replicates*.01)):
		for y in range(len(startfreqers2)):
			s_posteriors.write(str(valuecums2[y][x])+','+str(valuecums2het[y][x])+',')
		s_posteriors.write('\n')
else:
	for x in range(len(startfreqers2)):
		s_posteriors.write('Site_'+str(x+1)+',')
	s_posteriors.write('\n')
	for x in range(int(Part2Replicates*.01)):
		for y in range(len(startfreqers2)):
			s_posteriors.write(str(valuecums2[y][x])+',')
		s_posteriors.write('\n')

s_posteriors.close()
s_output_file.close()

