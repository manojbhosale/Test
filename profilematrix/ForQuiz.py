import re

def getKmerProbability(kmer):
	kmerDict = {"A":0,"C":1,"G":2,"T":3}
	prob = 1
	for i in range(len(kmer)):
		prob *= profilematrix[kmerDict[kmer[i]]][i]
	return prob

def getMostProbKmer(s, la):
	matches = re.findall( r'(?=(\w{%s}))' % la,s)
	prob = 0
	reskmer = ""
	for kmer in matches:
		if prob < getKmerProbability(kmer):
			prob = getKmerProbability(kmer)
			reskmer = kmer
	return reskmer
	
with open("input.txt","r") as infile:
	lines = infile.readlines()
	seq = lines[0].strip()
	k = int(lines[1].strip())
	linecnt = 0
	
	profilematrix=[[0 for i in range(k)]for i in range(k)]
	rowCnt = 0
	colCnt = 0
	for line in lines:
		if linecnt < 2:
			linecnt+=1
			continue
		values = map(lambda x: float(x), line.strip().split("  "))
		colCnt = 0
		for i in values:
			print(rowCnt," ",colCnt)
			profilematrix[rowCnt][colCnt] = i 
			colCnt+=1
		rowCnt+=1
		

	print(getMostProbKmer(seq, k))
		
print(getKmerProbability("TCGGTA"))
""" Read 2D matix
for i in range(rowCnt):
	for j in range(k):
		print(profilematrix[i][j]," ", end=" ")
	print()
"""







