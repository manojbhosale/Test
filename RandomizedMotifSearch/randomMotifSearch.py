import re
from random import randint

def getkmers(seq, k):
	return re.findall("(?=(\w{%s}))" % k,seq)

def getBestMotifs(dna, k):
	bestmotifs= list()
	for seq in dna:
		prefix = seq[:k]
		bestmotifs.append(prefix)
	return bestmotifs
	
def getKmerProbability(kmer, profilematrix):
	kmerDict = {"A":0,"C":1,"T":2,"G":3}
	prob = 1
	for i in range(len(kmer)):
		#print(kmer[i])
		prob *= profilematrix[kmerDict[kmer[i]]][i]
	return prob

def getMostProbKmer(s, k, profilematrix):
	matches = re.findall( r'(?=(\w{%s}))' % k,s)
	prob = 0
	reskmer = matches[0]
	#print(reskmer)
	for kmer in matches:
		#print(kmer," -> ",getKmerProbability(kmer, profilematrix))
		if prob < getKmerProbability(kmer, profilematrix):
			prob = getKmerProbability(kmer, profilematrix)
			reskmer = kmer
	return reskmer
	
def createProfileMatrix(motifs, k):
	#print(motifs)
	basesDict = {"A":0,"C":0,"T":0,"G":0}
	profilematrix=[[0 for x in range(k)] for x in range(4)]
	t = len(motifs)
	for i in range(k):
		#posCounter = 0
		basesDict = {"A":0,"C":0,"T":0,"G":0}
		#print(motifs)
		for seq in motifs:
			basesDict[seq[i]]+=1
		totalbasesIncolumn = (basesDict["A"]+basesDict["C"]+basesDict["G"]+basesDict["T"])
		profilematrix[0][i]=(basesDict["A"]+1)/(t+totalbasesIncolumn)
		profilematrix[1][i]=(basesDict["C"]+1)/(t+totalbasesIncolumn)
		profilematrix[2][i]=(basesDict["T"]+1)/(t+totalbasesIncolumn)
		profilematrix[3][i]=(basesDict["G"]+1)/(t+totalbasesIncolumn)
		#posCounter+=1
	return profilematrix
	
def findConsensus(motifs):
	basesDict = {"A":0,"C":0,"T":0,"G":0}
	k = len(motifs[0])
	selected = list()
	for i in range(k):
		basesDict = {"A":0,"C":0,"T":0,"G":0}
		for seq in motifs:
			basesDict[seq[i]]+=1
		maxval = max(basesDict.values())
		for key, value in basesDict.items():
			if value == maxval:
				selected.append(key)
				break
	return "".join(selected)

def HammingDistance(one, two):
	HammingDistance = 0
	i = 0;
	for base in one:
		if base is not two[i]:
			HammingDistance += 1
		i += 1
	return HammingDistance
	
def scoreMotifs(motifs):
	consensus = findConsensus(motifs)
	#print(motifs,"-->",consensus)
	score = 0
	for motif in motifs:
		score += HammingDistance(consensus, motif)
	return score
	
def createRandomMotifs(dna):
	randMotifs=[]
	for seq in dna:
		randomPos = randint(0,len(seq)-k-1)
		kmerlist = getkmers(seq,k)
		#print(randomPos)
		motif = kmerlist[randomPos]
		randMotifs.append(motif)
	return randMotifs
	
with open("input.txt","r") as infile:
	lines = infile.readlines()
	k, t = map(lambda x:int(x),lines[0].strip().split(" "))
	print("k = ",k," t = ",t)
	linecnt = 0
	dna =list()
	for line in lines:
		if linecnt == 0:
			linecnt+=1
			continue
		line.strip()
		dna.append(line)
#bestMotifs = getBestMotifs(dna, k)

bestMotifs = list()
baseseq = dna[0]
otherseq = dna[1:t]
kmerlist = getkmers(baseseq,k)
motifs = list()
flag = 0
round = 0
for times in range(1000):
	#randMotifs=[]
	randmotifs = createRandomMotifs(dna)
	#profilematrix = createProfileMatrix(motifs, k)
	#motifs=list()	
	if round == 0:
		round+=1	
		bestMotifs = randmotifs
	profilematrix = createProfileMatrix(randmotifs, k)
	while 1:
		
		motifs = list()
		for seq in dna:
			
			nextMotif = getMostProbKmer(seq, k, profilematrix)
			
			motifs.append(nextMotif)
			
		if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
			bestMotifs = motifs
			#print(scoreMotifs(bestMotifs))
		else:
			flag =1
			break
	if flag ==1:
		break
		
for	motif in bestMotifs:
	print(motif)


print(scoreMotifs(bestMotifs))

