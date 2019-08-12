"""
	Specify the genes to look at, and look in the provided set of patients which have a mutation in that gene, and which do not. 

"""

import sys
import numpy as np
from os import listdir
from os.path import isfile, join

#Read the expression data file
#Get all the patients for which we have RNA-seq data
expressionDataFile = sys.argv[1]
patientsWithExprData = []
with open(expressionDataFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		
		splitLine = line.split("\t")
		
		#The SNP file uses only the first 4 parts of the ID
		for patientID in splitLine[1:len(splitLine)-1]: 
			splitID = patientID.split("-")
			newID = "-".join(splitID[0:4])
			patientsWithExprData.append(newID[:-1]) #Also remove the final character for the IDS to match
		
		break

#Get the expression data in numpy format
expressionData = []
samples = []
with open(expressionDataFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		if lineCount == 0:
			samples = line.split("\t")
			lineCount += 1
			continue
		if lineCount < 2:
			lineCount += 1
			continue
		splitLine = line.split("\t")
		fullGeneName = splitLine[0]
		geneName = fullGeneName.split("|")[0]

		data = splitLine[1:len(splitLine)-1] 
		fixedData = [geneName]
		fixedData += data
		expressionData.append(fixedData)

expressionData = np.array(expressionData, dtype="object")	
print expressionData

cosmicGeneFile = sys.argv[2]
cosmicGenes = []

with open(cosmicGeneFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			lineCount += 1
		cosmicGenes.append(line.split("\t")[0])
		lineCount += 1
	
#Get all genes instead of just cosmic genes
geneList = []
with open(sys.argv[5], 'r') as geneFile:
	
	lineCount = 0
	for line in geneFile:
		line = line.strip()
		splitLine = line.split("\t")

		#Obtain the name, chromosome and positions of the gene. 
		
		geneID = splitLine[3]
		
		geneList.append(geneID)
		

#use all genes instead of cosmic
cosmicGenes = geneList

#genes = ['CTCF', 'RAD21', 'REC8', 'SMC1A', 'SMC1B', 'SMC3', 'STAG1', 'STAG2', 'STAG3']
#Split into cohesin/CTCF only
genes = ['CTCF']
geneLocations = [['chr16', 67596310, 67673088]]
#genes = ['RAD21', 'REC8', 'SMC1A', 'SMC1B', 'SMC3', 'STAG1', 'STAG2', 'STAG3']

snvDir = sys.argv[3]

#Collect all patients with mutations in these genes, even though we do not have matching expression data for all of these patients

mutatedGenesCount = dict()
mutatedGenesCountWithExpression = dict()

patientsWithSNVMutations = []
patientsWithoutSNVMutations = []

allFiles = [f for f in listdir(snvDir) if isfile(join(snvDir, f))]
print len(allFiles)

patientsNoExpr = []
mutations = []
for currentFile in allFiles:
	
	if currentFile == "MANIFEST.txt":
		continue
	splitFileName = currentFile.split(".")
	patientID = splitFileName[0]

	if patientID not in patientsWithExprData:
		patientsNoExpr.append(patientID)
	#Load the contents of the file
	with open(snvDir + "/" + currentFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip() #remove newlines
			if lineCount < 1: #only read the line if it is not a header line
				lineCount += 1
				continue

			splitLine = line.split("\t")
			geneName = splitLine[0]
		
			if geneName in genes:
				if geneName not in mutatedGenesCount:
					mutatedGenesCount[geneName] = 0
				mutatedGenesCount[geneName] += 1
				
				
				
				if patientID in patientsWithExprData: #Check if we also have expression data for this patient
					mutations.append(splitLine) 
					if patientID not in patientsWithSNVMutations:
						patientsWithSNVMutations.append(patientID)
						if geneName not in mutatedGenesCountWithExpression:
							mutatedGenesCountWithExpression[geneName] = 0
						mutatedGenesCountWithExpression[geneName] += 1
			else: #if the patient does not have the mutation
				if patientID in patientsWithExprData:
					if patientID not in patientsWithoutSNVMutations:
						patientsWithoutSNVMutations.append(patientID)

mutations = np.array(mutations, dtype="object")
print len(patientsNoExpr)
#np.savetxt("mutations_ctcf.txt", mutations, delimiter="\t", fmt='%s')

#Get patients with CNVs (deletions primarily)
cnvFile = sys.argv[4]
patientsWithDeletions = []
patientsWithoutDeletions = []
with open(cnvFile, 'r') as f:
	
	lineCount = 0
	for line in f:
		if lineCount < 1:
			lineCount += 1
			continue

		line = line.strip()
		splitLine = line.split("\t")
		
		#convert the sample name to the same format as for the snvs
		splitSampleName = splitLine[0].split("-")
		
		code = int("".join(list(splitSampleName[3])[0:2])) #skip normal sample portions
		if code > 9:
			continue
		
		#convert the sample name to the same format as for the snvs
		firstSamplePart = "-".join(splitSampleName[0:3])
		shortSampleName = firstSamplePart + '-' + splitSampleName[6]
		
		if shortSampleName not in patientsWithExprData: #Check if we also have expression data for this patient
			continue
		#Check if there is a change at the location containing CTCF
		if 'chr' + splitLine[1] == geneLocations[0][0]:
			
			if int(splitLine[2]) <= geneLocations[0][2] and int(splitLine[3]) >= geneLocations[0][1]: #CTCF overlaps this segment
				if float(splitLine[5]) < 0: #assume that this is a deletion
					if shortSampleName not in patientsWithDeletions:
						patientsWithDeletions.append(shortSampleName)
						continue
					
		if shortSampleName not in patientsWithoutDeletions:
			patientsWithoutDeletions.append(shortSampleName)

print("patients with deletions: ", len(patientsWithDeletions))
print("patients without deletions: ", len(patientsWithoutDeletions))

		
#patientsWithMutations = np.unique(patientsWithDeletions + patientsWithSNVMutations)
#patientsWithoutMutations = np.unique(patientsWithoutDeletions + patientsWithoutSNVMutations)

patientsWithMutations = patientsWithSNVMutations
patientsWithoutMutations = patientsWithoutSNVMutations

#merge patients with SNV mutation patients
print("number of patients with both SNVs and deletions: ", len(patientsWithMutations))
print("number of patients without both SNVs and deletions: ", len(patientsWithoutMutations))


#T-test for the overall case
from statsmodels.sandbox.stats.multicomp import multipletests
from math import sqrt
from scipy.stats import sem
from scipy.stats import t
from scipy import stats
def tTest(data1,data2):
	alpha = 0.05
	# calculate means
	mean1, mean2 = np.mean(data1), np.mean(data2)
	# calculate standard errors
	se1, se2 = sem(data1), sem(data2)
	# standard error on the difference between the samples
	sed = sqrt(se1**2.0 + se2**2.0)
	# calculate the t statistic
	t_stat = (mean1 - mean2) / sed
	# degrees of freedom
	df = len(data1) + len(data2) - 2
	# calculate the critical value
	cv = t.ppf(1.0 - alpha, df)
	# calculate the p-value
	p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
	# return everything
	return p

def tTestOneSided(data1,data2):
	
	t_stat, p = stats.ttest_ind(data1,data2, equal_var = False)
	print "t_stat: ", t_stat
	return p/2, t_stat
	
	alpha = 0.05
	# calculate means
	mean1, mean2 = np.mean(data1), np.mean(data2)
	# calculate standard errors
	se1, se2 = sem(data1), sem(data2)
	# standard error on the difference between the samples
	sed = sqrt(se1**2.0 + se2**2.0)
	# calculate the t statistic
	t_stat = (mean1 - mean2) / sed
	# degrees of freedom
	df = len(data1) + len(data2) - 2
	# calculate the critical value
	cv = t.ppf(1.0 - alpha, df)
	# calculate the p-value
	p = (1.0 - t.cdf(abs(t_stat), df))
	# return everything
	return p, t_stat

### gene-sample approach

cosmicGeneTStat = []
cosmicGenePValuesOneSided = []
geneInd = 0
expression = dict()
for gene in cosmicGenes:
	print "gene: ", gene, geneInd
	geneInd += 1
	if gene not in expressionData[:,0]:
		continue
	
	geneExpression = expressionData[expressionData[:,0] == gene][0]
	
	negativeExpr = []
	
	#First collect the negative set, all patients that do not have a mutation in this cosmic gene
	for patient in patientsWithoutMutations:
		splitPatientId = patient.split("-")
		patientId = splitPatientId[2]
		#Get the index of this sample, match the samplenames to each other that are vastly different
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if sample == "Hybridization REF":
				continue
			splitSampleName = sample.split("-")
			sampleId = splitSampleName[2]
			
			if patientId == sampleId:
				#if the samples are the same, use the samples with the right code only to make the negative set
				#only code < 10 are tumor samples, the rest is normal
				code = int("".join(list(splitSampleName[3])[0:2]))
				
				if code < 10:
					negativeExpr.append(float(geneExpression[sampleInd]))
			
			
		
	#repeat for the patients with the mutations
	for patient in patientsWithMutations:
		splitPatientId = patient.split("-")
		patientId = splitPatientId[2]
		#Get the index of this sample, match the samplenames to each other that are vastly different
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if sample == "Hybridization REF":
				continue
			splitSampleName = sample.split("-")
			sampleId = splitSampleName[2]
			
			if patientId == sampleId:
				#if the samples are the same, use the samples with the right code only to make the negative set
				#only code < 10 are tumor samples, the rest is normal
				code = int("".join(list(splitSampleName[3])[0:2]))
				
				
				if code < 10:
					
					#Do t-test for this patient against all negatives
					z = (float(geneExpression[sampleInd]) - np.mean(negativeExpr)) / float(np.std(negativeExpr))
					pValue = stats.norm.sf(abs(z))
					
					cosmicGenePValuesOneSided.append([gene + "_" + sample, pValue])
					cosmicGeneTStat.append([gene + "_" + sample, z])
					
					#expression[gene + "_" + sample] = [float(geneExpression[sampleInd]), np.mean(negativeExpr)]
					expression[gene + "_" + sample] = [float(geneExpression[sampleInd]), negativeExpr]			


cosmicGenePValuesOneSided = np.array(cosmicGenePValuesOneSided, dtype="object")
sortedInd = cosmicGenePValuesOneSided[:,1].argsort()
cosmicGenePValuesOneSided = cosmicGenePValuesOneSided[sortedInd]	
cosmicGeneTStat = np.array(cosmicGeneTStat, dtype="object")	
cosmicGeneTStat = cosmicGeneTStat[sortedInd]

print cosmicGenePValuesOneSided
print cosmicGeneTStat

reject, pAdjusted, _, _ = multipletests(cosmicGenePValuesOneSided[:,1], method='bonferroni') #fdr_bh or bonferroni

import matplotlib.pyplot as plt

print "Significant COSMIC genes after bonferroni and 1-sided: "
filteredPValues = []
signGenes = []
for pValueInd in range(0, len(cosmicGenePValuesOneSided)):
	
	if reject[pValueInd] == True and np.sign(cosmicGeneTStat[pValueInd,1]) == 1:
		
		pValue = pAdjusted[pValueInd]
		filteredPValues.append([cosmicGenePValuesOneSided[pValueInd,0], cosmicGenePValuesOneSided[pValueInd,1]])
		pair = cosmicGenePValuesOneSided[pValueInd,0]
		splitPair = pair.split("_")
		
		#print splitPair[0] + "\t" + splitPair[1] + "\t" + str(expression[pair][0]) + "\t" + str(expression[pair][1]) + "\t" + str(pAdjusted[pValueInd]) + "\t" + str(cosmicGeneTStat[pValueInd,1])
		#plt.clf()
		#plt.hist(expression[pair][1])
		#plt.axvline(expression[pair][0], color='k', linestyle='dashed', linewidth=1)
		#plt.savefig('expressionPlots/' + pair + '.svg')
		
		# print cosmicGeneTStat[pValueInd,1]
		signGenes.append(cosmicGenePValuesOneSided[pValueInd,0])
filteredPValues = np.array(filteredPValues, dtype="object")

np.savetxt("UCEC_significantAllGenes_bonferroni_oneS_snvs.txt", filteredPValues, fmt="%s", delimiter="\t")


exit()
#### Multiple sample approach


#For each cosmic gene, collect the expression values of the patients in the one group and the other group.
cosmicGenePValues = []
cosmicGeneTStat = []
cosmicGenePValuesOneSided = []
geneInd = 0
expressionMatrix = []
for gene in cosmicGenes:
	print "gene: ", gene, geneInd
	geneInd += 1
	if gene not in expressionData[:,0]:
		continue
	
	geneExpression = expressionData[expressionData[:,0] == gene][0]
	
	negativeExpr = []
	
	#First collect the negative set, all patients that do not have a mutation in this cosmic gene
	for patient in patientsWithoutMutations:
		splitPatientId = patient.split("-")
		patientId = splitPatientId[2]
		#Get the index of this sample, match the samplenames to each other that are vastly different
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if sample == "Hybridization REF":
				continue
			splitSampleName = sample.split("-")
			sampleId = splitSampleName[2]
			
			if patientId == sampleId:
				#if the samples are the same, use the samples with the right code only to make the negative set
				#only code < 10 are tumor samples, the rest is normal
				code = int(splitSampleName[len(splitSampleName)-1])
				 
				if code < 10:
					negativeExpr.append(float(geneExpression[sampleInd]))
			
		
	#repeat for the patients with the mutations
	positiveExpr = []
	positiveSamples = []
	for patient in patientsWithMutations:
		splitPatientId = patient.split("-")
		patientId = splitPatientId[2]
		#Get the index of this sample, match the samplenames to each other that are vastly different
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if sample == "Hybridization REF":
				continue
			splitSampleName = sample.split("-")
			sampleId = splitSampleName[2]
			
			if patientId == sampleId:
				#if the samples are the same, use the samples with the right code only to make the negative set
				#only code < 10 are tumor samples, the rest is normal
				code = int(splitSampleName[len(splitSampleName)-1])
				 
				if code < 10:
					positiveExpr.append(float(geneExpression[sampleInd]))
					positiveSamples.append(sample)
		
	#Do a t-test to see if the expression is different for any of the COSMIC genes
	pValue = tTest(positiveExpr, negativeExpr)
	pValueOneSided, tStat = tTestOneSided(positiveExpr, negativeExpr)
	cosmicGenePValues.append([gene, pValue])
	cosmicGenePValuesOneSided.append([gene, pValueOneSided])
	cosmicGeneTStat.append([gene, tStat])	

print positiveSamples
expressionMatrix = np.array(expressionMatrix, dtype="object")
np.savetxt('expressionMatrix.txt', expressionMatrix, fmt='%s', delimiter='\t')
	
cosmicGenePValues = np.array(cosmicGenePValues, dtype="object")	
cosmicGenePValues = cosmicGenePValues[cosmicGenePValues[:,1].argsort()]
cosmicGeneTStat = np.array(cosmicGeneTStat, dtype="object")	
cosmicGeneTStat = cosmicGeneTStat[cosmicGeneTStat[:,1].argsort()]
# 
# reject, pAdjusted, _, _ = multipletests(cosmicGenePValues[:,1], method='bonferroni') #fdr_bh or bonferroni
# 
# filteredPValues = cosmicGenePValues[reject]
# print "Significant COSMIC genes after bonferroni: "
# for pValue in filteredPValues:
# 	print pValue[0], ": ", pValue[1]
# 
# np.savetxt("significantCosmicGenes_bonferroni.txt", filteredPValues, fmt="%s", delimiter="\t")
# 
# reject, pAdjusted, _, _ = multipletests(cosmicGenePValues[:,1], method='fdr_bh') #fdr_bh or bonferroni
# 
# filteredPValues = cosmicGenePValues[reject]
# print "Significant COSMIC genes after BH: "
# for pValue in filteredPValues:
# 	print pValue[0], ": ", pValue[1]
# 
# np.savetxt("significantCosmicGenes_bh.txt", filteredPValues, fmt="%s", delimiter="\t")

#Repeat for one-sided t-test
	
cosmicGenePValuesOneSided = np.array(cosmicGenePValuesOneSided, dtype="object")	
cosmicGenePValuesOneSided = cosmicGenePValues[cosmicGenePValuesOneSided[:,1].argsort()]

reject, pAdjusted, _, _ = multipletests(cosmicGenePValuesOneSided[:,1], method='bonferroni') #fdr_bh or bonferroni

print "Significant COSMIC genes after bonferroni and 1-sided: "
filteredPValues = []
signGenes = []
for pValueInd in range(0, len(cosmicGenePValuesOneSided)):
	if reject[pValueInd] == True and np.sign(cosmicGeneTStat[pValueInd,1]) == 1:
		
		pValue = pAdjusted[pValueInd]
		filteredPValues.append([cosmicGenePValuesOneSided[pValueInd,0], cosmicGenePValuesOneSided[pValueInd,1]])
		print cosmicGenePValuesOneSided[pValueInd,0], ": ", cosmicGenePValuesOneSided[pValueInd,1]
		print cosmicGeneTStat[pValueInd,1]
		signGenes.append(cosmicGenePValuesOneSided[pValueInd,0])
filteredPValues = np.array(filteredPValues, dtype="object")
print signGenes
np.savetxt("significantCosmicGenes_bonferroni_oneS.txt", filteredPValues, fmt="%s", delimiter="\t")
# 
# reject, pAdjusted, _, _ = multipletests(cosmicGenePValuesOneSided[:,1], method='fdr_bh') #fdr_bh or bonferroni
# 
# print "Significant COSMIC genes after BH and 1-sided: "
# filteredPValues = []
# for pValueInd in range(0, len(cosmicGenePValuesOneSided)):
# 	if reject[pValueInd] == True and np.sign(cosmicGeneTStat[pValueInd,1]) == 1:
# 		pValue = pAdjusted[pValueInd]
# 		filteredPValues.append([cosmicGenePValuesOneSided[pValueInd,0], cosmicGenePValuesOneSided[pValueInd,1]])
# 		print cosmicGenePValuesOneSided[pValueInd,0], ": ", cosmicGenePValuesOneSided[pValueInd,1]
# filteredPValues = np.array(filteredPValues, dtype="object")
# 
# np.savetxt("significantCosmicGenes_bh_oneS.txt", filteredPValues, fmt="%s", delimiter="\t")
# 
# 
