import re
import operator
from itertools import *


def getKmers(s, la):
	matches = re.findall( r'(?=(\w{%s}))' % la,s)
	return matches

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
			
			

def checkOccur(seq, k, d):
	kmerlist = getKmers(seq, k)
	kmerset = set()
	for kmer in kmerlist:
		kmerset.add(kmer)
		mismatchedKmers = generate(kmer,d)
		for mmkmer in mismatchedKmers:
			kmerset.add(mmkmer)
	return kmerset
			
with open("input.txt","r") as infile:
	
	alllines = infile.readlines();
	
linecount = 0
seqcount = 0;
occurlist=[]
newset = set()
oldset = set()
resultset = set()

for aline in alllines:
	if(linecount == 0):
		linecount+=1
		k,d = map(lambda x: int(x) ,aline.strip().split())
		continue
	#print(aline.strip())
	seq = aline.strip()
	newset = checkOccur(seq, k, d)	
	if(seqcount == 0):
		#oldset = newset
		resultset = newset
		seqcount+=1
		continue
	resultset = resultset & newset
	seqcount+=1

with open("result.txt","w") as outfile:
	print(resultset)	
	outfile.write(' '.join(resultset))