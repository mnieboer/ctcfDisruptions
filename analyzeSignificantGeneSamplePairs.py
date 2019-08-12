"""
	For the significant gene-sample pairs, check how often each gene is found. Are some genes found across multiple samples?
	
	Use this script to output a table showing for each pair if:
	- The pair is significant
	- The p-value
	- If the gene is in COSMIC
	- In how many samples the gene is found in total
	- If there is an SNV or SV affecting this gene in the sample
	- If the expression change can be explained by the SNV/SV yes/no
	- The significance of the expression change. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import re
from scipy import stats

pairList = np.loadtxt(sys.argv[1], dtype='object')

geneOccurrenceCounts = dict()
shortSampleNamePairs = [] #make a list with short sample names to check for SNVs
shortSampleNameLookup = dict()
signPatientIds = []
for pair in pairList[:,0]:
	
	splitPair = pair.split("_")
	signPatientIds.append(splitPair[1])
	
	if splitPair[0] not in geneOccurrenceCounts:
		geneOccurrenceCounts[splitPair[0]] = 0
	geneOccurrenceCounts[splitPair[0]] += 1
	
	splitSampleName = splitPair[1].split("-")
	firstSamplePart = "-".join(splitSampleName[0:3])
	#shortSampleName = firstSamplePart + '-' + splitSampleName[6]
	shortSampleName = firstSamplePart + '-' + '01' #the other data ends with 01 for some reason, but that is another portion of the sample, so the final number is not relevant
	shortSampleNamePairs.append(splitPair[0] + "_" + shortSampleName)
	shortSampleNameLookup[splitPair[1]] = shortSampleName

multipleTimesCounter = 0
timesDistribution = dict()
for gene in geneOccurrenceCounts:
	
	if geneOccurrenceCounts[gene] not in timesDistribution:
		timesDistribution[geneOccurrenceCounts[gene]] = 0
	timesDistribution[geneOccurrenceCounts[gene]] += 1
	if geneOccurrenceCounts[gene] > 1:
		#print("gene: ", gene, " found ", geneOccurrenceCounts[gene], " times")
		multipleTimesCounter += 1
		
print(multipleTimesCounter)
		
plt.bar(timesDistribution.keys(), timesDistribution.values())
#plt.show()


#Check: which of the significant genes also have SNVs or coding SVs in the sample?
#What we can also check right away: can we also see a difference in expression between this sample and the expression in all samples that do not have an SV/SNV?

snvDir = sys.argv[2]
allFiles = [f for f in listdir(snvDir) if isfile(join(snvDir, f))]
print len(allFiles)

degsWithSNVEvidence = []
snvPatientIds = []
for currentFile in allFiles:
	
	if currentFile == "MANIFEST.txt":
		continue
	splitFileName = currentFile.split(".")
	patientID = splitFileName[0]

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
			
			pair = geneName + "_" + patientID
			snvPatientIds.append(pair)
		
			if pair in shortSampleNamePairs:
				
				#if this is a gene that is DEG, and there is an SNV in the gene, this gene may be DEG due to coding effects. 
				degsWithSNVEvidence.append(pair)
			
print("number of genes with SNVs:", len(degsWithSNVEvidence))

#Determine if the expression change for this gene is significant as a result of the SNV.
#1. Get all patients that do not have an SNV for this gene
#For every gene, get a list of all samples in which this gene is affected to exclude these and make a null distribution

expressionFile = sys.argv[3]

expressionData = []
samples = []
with open(expressionFile, 'r') as inF:
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
snvDEGPairs = []
if len(degsWithSNVEvidence) > 0:
	geneSampleRef = dict()
		
	for pair in degsWithSNVEvidence:
		splitPair = pair.split("_")
		gene = splitPair[0]
		sample = splitPair[len(splitPair)-1]
		
		if gene not in geneSampleRef:
			geneSampleRef[gene] = []
		geneSampleRef[gene].append(sample)
	#Get the epxression of all possible genes
	geneSampleExpr = dict()
	for gene in geneSampleRef:
		
		if gene not in expressionData[:,0]:
			continue
		
		geneSamples = geneSampleRef[gene]
		geneSampleExpr[gene] = dict()
		geneExpression = expressionData[expressionData[:,0] == gene][0]
		for geneSample in geneSamples:
			
			shortSampleName = geneSample
			#shortSampleName = '-'.join(geneSample.split("-")[0:3]) + '-01'
			
			#match the sample name with the expression sample name
			for sampleInd in range(0, len(samples)):
				sample = samples[sampleInd]
				if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
					
					splitSampleName = sample.split("-")
					code = int("".join(list(splitSampleName[3])[0:2]))
					
					if code < 10: #above 9 are the normal samples, which we do not want to include here
						sampleInd = samples.index(sample)
						
						geneSampleExpr[gene][geneSample] = float(geneExpression[sampleInd])
	print("done getting expr for samples")
		
	#2. Make a negative expression set
	negativeExpr = dict()
	for gene in geneSampleExpr:
		matchedFullSampleNames = list(geneSampleExpr[gene].keys())
		
		#Get all the samples without an SV for this gene
		unmatchedSamples = np.setdiff1d(samples[1:len(samples)-1], matchedFullSampleNames) #exclude hybrid ref
		negativeSamples = []
		for sample in unmatchedSamples: #sample tumor samples, exclude normals
			splitSampleName = sample.split("-")
			code = int("".join(list(splitSampleName[3])[0:2]))
			
			if code < 10: 
				negativeSamples.append(sample)
			
		#Get the expression of these samples
		negativeSampleExpressionValues = []
		for sample in negativeSamples:
			sampleInd = samples.index(sample)				
			negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
		
		negativeExpr[gene] = negativeSampleExpressionValues
	print("negative expr done")
	
	#3. Compare for each gene with the SNV if the expression is significant
	def getDEPairs(pairs, geneSampleRef, epressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr):
		zScores = dict()								
		for pair in pairs:
			splitPair = pair.split("_")
			gene = splitPair[0]
			pairSample = splitPair[len(splitPair)-1]
			#shortPairSampleName = pairSample
			#shortSampleName = '-'.join(geneSample.split("-")[0:3]) + '-01'
			
			sv = "_".join(splitPair[1:])
			if gene not in expressionData[:,0]:
				continue
			
			
			sampleExpressionValue = geneSampleExpr[gene][pairSample] #expression values of this gene in all samples
			matchedFullSampleNames = list(geneSampleExpr[gene].keys())
						
			
			negativeSampleExpressionValues = negativeExpr[gene]
			
			#Get the expression z-score for this pair
			if np.std(negativeSampleExpressionValues) == 0:
				continue
		
			z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
			pValue = stats.norm.sf(abs(z))
		
			perPairDifferentialExpression[pair] = pValue
			zScores[pair] = z
			
		return perPairDifferentialExpression, zScores
	
	#Get the p-value for each pair in coding & non-coding
	perPairDifferentialExpression, zScoresSnv = getDEPairs(degsWithSNVEvidence, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
	print("done")
	
	#4. multiple testing correction
	perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
	perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
	perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())
	
	from statsmodels.sandbox.stats.multicomp import multipletests
	reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')
	
	perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]
	
	#5. Report. 
	snvDEGPairs = perPairDifferentialExpressionArrayFiltered
	print snvDEGPairs

#Repeat but then for coding SVs.
#These are all SVs that overlap at least 1 gene, output from the rule-based prioritizer. 
codingPairs = np.loadtxt(sys.argv[4], dtype='object')
degsWithSVEvidence = []

codingPairsShortSampleNames = []

for pair in codingPairs:
	
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	codingPairsShortSampleNames.append(gene + "_" + sample)
	

for degPair in pairList[:,0]:
	
	splitDegPair = degPair.split("_")
	degGene = splitDegPair[0]
	#convert the sample name to the format in the SV data for BRCA
	splitSampleName = splitDegPair[1].split("-")
	sampleCode = splitSampleName[2]
	
	convertedSampleName = 'brca' + sampleCode
	
	#Convert for UCEC
	#firstSamplePart = "-".join(splitSampleName[0:3])
	#convertedSampleName = firstSamplePart

	if degGene + "_" + convertedSampleName in codingPairsShortSampleNames:
		
		if degPair not in degsWithSVEvidence:
			#re-do this here for the brca case, where the identifiers are weird in the sv input file
			firstSamplePart = "-".join(splitSampleName[0:3])
			convertedSampleName = firstSamplePart + '-01'
			degsWithSVEvidence.append(degGene + "_" + convertedSampleName)
	
print("number of genes with SVs: ", len(degsWithSVEvidence))
svDEGPairs = []
if len(degsWithSVEvidence) > 0:
	#Also check which genes could be DEG due to the SV. 
	geneSampleRef = dict()
	for pair in degsWithSVEvidence:
		splitPair = pair.split("_")
		gene = splitPair[0]
		sample = splitPair[len(splitPair)-1]
		
		if gene not in geneSampleRef:
			geneSampleRef[gene] = []
		geneSampleRef[gene].append(sample)
	#Get the epxression of all possible genes
	geneSampleExpr = dict()
	for gene in geneSampleRef:
		
		if gene not in expressionData[:,0]:
			continue
		
		geneSamples = geneSampleRef[gene]
		geneSampleExpr[gene] = dict()
		geneExpression = expressionData[expressionData[:,0] == gene][0]
		for geneSample in geneSamples:
			
			
			shortSampleName = geneSample
			#shortSampleName = '-'.join(geneSample.split("-")[0:3]) + '-01'
			
			#match the sample name with the expression sample name
			for sampleInd in range(0, len(samples)):
				sample = samples[sampleInd]
				if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
					
					splitSampleName = sample.split("-")
					code = int("".join(list(splitSampleName[3])[0:2]))
					
					if code < 10: #above 9 are the normal samples, which we do not want to include here
						sampleInd = samples.index(sample)
						
						geneSampleExpr[gene][geneSample] = float(geneExpression[sampleInd])
	print("done getting expr for samples")
		
	#2. Make a negative expression set
	negativeExpr = dict()
	for gene in geneSampleExpr:
		matchedFullSampleNames = list(geneSampleExpr[gene].keys())
		
		#Get all the samples without an SV for this gene
		unmatchedSamples = np.setdiff1d(samples[1:len(samples)-1], matchedFullSampleNames) #exclude hybrid ref
		negativeSamples = []
		for sample in unmatchedSamples: #sample tumor samples, exclude normals
			splitSampleName = sample.split("-")
			code = int("".join(list(splitSampleName[3])[0:2]))
			
			if code < 10: 
				negativeSamples.append(sample)
			
		#Get the expression of these samples
		negativeSampleExpressionValues = []
		for sample in negativeSamples:
			sampleInd = samples.index(sample)				
			negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
		
		negativeExpr[gene] = negativeSampleExpressionValues
	print("negative expr done")
	
	#3. Compare for each gene with the SNV if the expression is significant
	def getDEPairs(pairs, geneSampleRef, epressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr):
		zScores = dict()								
		for pair in pairs:
			splitPair = pair.split("_")
			gene = splitPair[0]
			pairSample = splitPair[len(splitPair)-1]
			shortPairSampleName = pairSample
			shortSampleName = '-'.join(pairSample.split("-")[0:3]) + '-01'
			
			sv = "_".join(splitPair[1:])
			if gene not in expressionData[:,0]:
				continue
			
			
			sampleExpressionValue = geneSampleExpr[gene][pairSample] #expression values of this gene in all samples
			matchedFullSampleNames = list(geneSampleExpr[gene].keys())
						
			
			negativeSampleExpressionValues = negativeExpr[gene]
			
			#Get the expression z-score for this pair
			if np.std(negativeSampleExpressionValues) == 0:
				continue
		
			z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
			pValue = stats.norm.sf(abs(z))
		
			perPairDifferentialExpression[gene + "_" + shortSampleName] = pValue
			zScores[gene + "_" + shortSampleName] = z
			
		return perPairDifferentialExpression, zScores
	
	#Get the p-value for each pair in coding & non-coding
	perPairDifferentialExpression, zScoresSv = getDEPairs(degsWithSVEvidence, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
	print("done")
	
	#4. multiple testing correction
	perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
	perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
	perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())
	
	from statsmodels.sandbox.stats.multicomp import multipletests
	reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')
	
	perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]
	
	#5. Report. 
	svDEGPairs = perPairDifferentialExpressionArrayFiltered
	print svDEGPairs


##Use all this information to generate an output table.

#get cosmic genes

cosmicGeneFile = sys.argv[5]
cosmicGenes = []

with open(cosmicGeneFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			lineCount += 1
		cosmicGenes.append(line.split("\t")[0])
		lineCount += 1
	

header = ''
pairInfo = []
for pair in pairList:
	
	#Get the gene, in how many samples it is DEG, and if it is in COSMIC.
	splitPair = pair[0].split("_")
	gene = splitPair[0]
	cosmic = 'False'
	if gene in cosmicGenes:
		cosmic = 'True'
	noOfSamples = geneOccurrenceCounts[gene]
	
	#Check if there is an SNV for this gene in any sample
	#This needs the original identifier
	shortSampleName = shortSampleNameLookup[splitPair[1]]
	snv = 'False'
	snvSignificant = 'False'
	snvSignificance = '-'
	if gene + '_' + shortSampleName in degsWithSNVEvidence:
		snv = 'True'
		if gene + '_' + shortSampleName in snvDEGPairs[:,0]:
			if np.sign(zScoresSnv[gene + '_' + shortSampleName]) == 1:
				snvSignificant = 'True'
				snvSignificance = snvDEGPairs[snvDEGPairs[:,0] == gene + '_' + shortSampleName,1][0]
				
	
	#Repeat for SVs
	
	shortSampleName = shortSampleNameLookup[splitPair[1]]
	
	sv = 'False'
	svSignificant = 'False'
	svSignificance = '-'
	if gene + '_' + shortSampleName in degsWithSVEvidence:
		sv = 'True'
		if gene + '_' + shortSampleName in svDEGPairs[:,0]:
			if np.sign(zScoresSv[gene + '_' + shortSampleName]) == 1:
				svSignificant = 'True'
				svSignificance = svDEGPairs[svDEGPairs[:,0] == gene + '_' + shortSampleName,1][0]
	
	#Collect all information about this pair
	pairInfo.append([splitPair[0], splitPair[1], pair[1], cosmic, noOfSamples, snv, snvSignificant, snvSignificance, sv, svSignificant, svSignificance])	

pairInfo = np.array(pairInfo, dtype='object')
header = 'gene\tsample\tp-value_CTCF_vs_non-CTCF\tCOSMIC\tnumber_of_samples_in_which_deg\tsnv\tsnv_significant\tsnv_p-value\tsv\tsv_significant\tsv_significance'
np.savetxt(sys.argv[6], pairInfo, header=header, fmt='%s', delimiter='\t')

