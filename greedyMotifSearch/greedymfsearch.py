import re

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
		profilematrix[0][i]=basesDict["A"]/t
		profilematrix[1][i]=basesDict["C"]/t
		profilematrix[2][i]=basesDict["T"]/t
		profilematrix[3][i]=basesDict["G"]/t
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
bestMotifs = getBestMotifs(dna, k)


"""
print(nextMotif)
motifs = ['GGC', 'AAG', 'AAG', 'GTC'] 
string02 = 'CAATAATATTCG'
profilematrix = createProfileMatrix(motifs, k)
print(profilematrix)
nextMotif = getMostProbKmer(string02, k, profilematrix)

Motifs = ['GGC']  
- string02 = 'AAGAATCAGTCA'
- Profile01:
   {'A': [0.0, 0.0, 0.0], 'C': [0.0, 0.0, 1.0], 'T': [0.0, 0.0, 0.0], 'G': [1.0, 1.0, 0.0]}
- profileMostProbableKmer(string02, 3, Profile01): 'AAG'
- profileProbableKmers(string02, 3, Profile01)
   [('AAG', 0.0), ('GTC', 0.0), ('ATC', 0.0), ('AGA', 0.0), ('AAT', 0.0), ('AGT', 0.0), ('CAG', 0.0), ('GAA', 0.0), ('TCA', 0.0)]
"""

baseseq = dna[0]
otherseq = dna[1:t]
kmerlist = getkmers(baseseq,k)
motifs = list()
for kmer in kmerlist:
	motifs = list()
	motifs.append(kmer)
	for seq in otherseq:
		
		profilematrix = createProfileMatrix(motifs, k)
		#print(profilematrix)
		nextMotif = getMostProbKmer(seq, k, profilematrix)
		#print(nextMotif)
		motifs.append(nextMotif)
		
	if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
		bestMotifs = motifs
	#print(motifs," ",bestMotifs)
	#print(scoreMotifs(motifs)," ",scoreMotifs(bestMotifs))
#motif = ['GGC', 'AAG', 'CAA', 'CAC', 'CAA']
#print(motif)
#print(scoreMotifs(motif))
print(bestMotifs)
"""
GREEDYMOTIFSEARCH(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string
                      from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string
                          in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        output BestMotifs

Score(motifs)
  // note: I broke findConsensus ties by alphabetical order,
  // but not sure if it matters or not
  consensus = findConsensus(motifs)
  score = 0
  for each (motif in motifs)
    score += hammingDistance(consensus, motif)
  output score
"""