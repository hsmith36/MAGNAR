#Python Script

import sys
from collections import defaultdict

print '\nRunning File: ', sys.argv[0]
print '\nUsing File: ', sys.argv[1]
outfile = open(sys.argv[2], "w")

#=====================================Define Functions Space====================================

#Gets all the name tags in the file
def setupNames(file):

	names = []
	for line in file:
		ln_lst = line.split(' ')
		for i in range(len(ln_lst)):
			if i > 0:
				nm_lst = ln_lst[i].split('|')
				nm = nm_lst[0]
				if nm not in names:
					names.append(nm)
	
	return names


#Creates a phenotype key to be put in the dictionary
def setupPheno(ln_lst):
	ln_lst.pop(0)#to remove the name 
	for i in range(len(ln_lst)):
		nm_lst = ln_lst[i].split('|')
		ln_lst[i] = nm_lst[0]

	pheno = [0]*len(names)

	for i in range(len(pheno)):
		if names[i] in ln_lst:
			pheno[i] = 1
	
	return pheno


#turns a list into a string separated by spaces
def listToString(pheno):
	pheno_str = ""
	for i in range(len(pheno)):
		pheno_str += str(pheno[i]) + " "
	
	return pheno_str
		

#Runs through the file and sets up the matrix line by line
def setupMatrix(file, names):

	matrix = defaultdict(list)
	pheno_name = [] #list of tuples containing the phenotype and name

	for line in file:
		ln_lst = line.split(' ')
		ln_name = ln_lst[0]
		ln_name = ln_name.replace(":", "")

		pheno = setupPheno(ln_lst)
		pheno_str = listToString(pheno)

		pheno_name.append((pheno_str, ln_name))

	for k, v in pheno_name:
		matrix[k].append(v)

	return matrix
		
#================================================================================================

file = open(sys.argv[1])
names = setupNames(file)

#prints 4-letter id names
all_names = ""
for nm in names:
	all_names += nm + " "
outfile.write(all_names + "|NAMES\n")

file = open(sys.argv[1])
matrix = setupMatrix(file, names)

#prints matrix
for pheno in matrix:
	println = ""
	println += pheno + "|"
	for name in matrix[pheno]:
		println += name + ","

	println = println[:len(println)-1] # to get rid of the final comma
	outfile.write(println + "\n")

outfile.close()
print '\nEnding the program\n'

