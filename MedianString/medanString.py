import re
from itertools import *


def generate(s,d):
	N = len(s)
	letters = "ACGT"
	pool = list(s)
	
	for indices in combinations(range(N), d):
		for replacements in product(letters, repeat=d):
			skip = False
			for i , a in zip(indices, replacements):
				if pool[i] == a: skip = True
			if skip: continue
		
			keys = dict(zip(indices, replacements))
			yield ''.join([pool[i] if i not in indices else keys[i] for i in range(N)])


def getkmers(seq, k):
	return re.findall("(?=(\w{%s}))" % k,seq)

def HammingDistance(one, two):
	HammingDistance = 0
	i = 0;
	for base in one:
		if base is not two[i]:
			#print(base," == ",two[i])
			HammingDistance += 1
		i += 1
	return HammingDistance
	
def DistanceBetweenPatternAndStrings(Pattern, Dna):
    k = len(Pattern)
    distance = 0
    for seq in Dna:
        HamDist = k
        for kmer in getkmers(seq,k):
            if HamDist > HammingDistance(Pattern, kmer):
                HamDist = HammingDistance(Pattern, kmer)
        distance = distance + HamDist
    return distance

def MedianString(kmers, Dna, k):
	distance = 1000000
	Median = ""
	for kmer in kmers:
		#print(DistanceBetweenPatternAndStrings(kmer, Dna))
		if distance > DistanceBetweenPatternAndStrings(kmer, Dna):
			distance = DistanceBetweenPatternAndStrings(kmer, Dna)
			Median = kmer
	return Median

def getPossibleKmers(Dna, k):
	kmerset = set()
	for seq in Dna:
		for kmer in getkmers(seq, k):
			kmerset.add(kmer)
			kmerset.union(generate(kmer, k))
	return kmerset
	
with open("input.txt","r") as infile:
	lines = infile.readlines()
	k = int(lines[0].strip())
	linecnt = 0
	i = 0
	Dna=[]
	for seq in lines:
		#print(seq)
		if linecnt == 0:
			linecnt+=1
			continue
		Dna.append(seq.strip())
		i+=1
#	print(Dna,k)
#print(DistanceBetweenPatternAndStrings(pattern, Dna))
#print(HammingDistance("ATC","ATG"))
#print(getkmers("ATCG",2))
kmers = getPossibleKmers(Dna, k)
#print(kmers)
print(MedianString(kmers, Dna, k ))