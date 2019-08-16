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
mutationInfo = dict() #for each gene-sample combination with an SNV, also keep the relevant information so that we can write it to the output table. 
for currentFile in allFiles:
	
	if currentFile == "MANIFEST.txt":
		continue
	splitFileName = currentFile.split(".")
	patientID = splitFileName[0]

	#Load the contents of the file
	with open(snvDir + "/" + currentFile, 'r') as inF:
		lineCount = 0
		header = dict()
		for line in inF:
			line = line.strip() #remove newlines
			splitLine = line.split("\t")
			if lineCount < 1: #only read the line if it is not a header line
				for fieldInd in range(0, len(splitLine)):
					header[splitLine[fieldInd]] = fieldInd
				lineCount += 1
				continue

			
			geneName = splitLine[0]
			
			pair = geneName + "_" + patientID
			snvPatientIds.append(pair)
			
			#position is usually denoted as 'chr:start-end'
			chromosome = splitLine[header['Chromosome']]
			if 'Start_Position' in header: #thanks TCGA for being consistent
				start = splitLine[header['Start_Position']]
				end = splitLine[header['End_Position']]
			else:
				start = splitLine[header['Start_position']]
				end = splitLine[header['End_position']]
			varType = splitLine[header['Variant_Classification']]
			ref = splitLine[header['Reference_Allele']]
			var = splitLine[header['Tumor_Seq_Allele2']]
			
			mutationInfo[pair] = [chromosome + ":" + start + "-" + end, varType, ref, var]
		
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
degsWithSVEvidenceSVInfo = [] #also keep more information about the SV for later reference. 

codingPairsShortSampleNames = dict()
for pair in codingPairs:
	
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	codingPairsShortSampleNames[gene + "_" + sample] = pair

for degPair in pairList[:,0]:
	
	splitDegPair = degPair.split("_")
	degGene = splitDegPair[0]
	#convert the sample name to the format in the SV data for BRCA
	splitSampleName = splitDegPair[1].split("-")
	sampleCode = splitSampleName[2]
	
	#convertedSampleName = 'brca' + sampleCode
	
	#Convert for UCEC
	firstSamplePart = "-".join(splitSampleName[0:3])
	convertedSampleName = firstSamplePart

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

svFile = sys.argv[6] #provide the orignal SV file to get SV info

svInfo = dict() #keep relevant SV info for each sample. The sample names will not match here on the pair info, this is missing from the sv file
with open(svFile, 'r') as f:
	
	lineCount = 0
	header = []
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
		if lineCount < 1:

			header = splitLine
			lineCount += 1
			continue
		
		#Now extract the chromosome, start and end (there are multiple, 2 for an SV)
		chr1Index = header.index("chr1")
		s1Index = header.index("s1")
		e1Index = header.index("e1")
		o1Index = header.index("o1")

		chr2Index = header.index("chr2")
		s2Index = header.index("s2")
		e2Index = header.index("e2")
		o2Index = header.index("o2")

		
		svTypeIndex = header.index("sv_type")
		svType = splitLine[svTypeIndex]
		
		#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
		if splitLine[chr1Index] == splitLine[chr2Index]:
			if splitLine[s2Index] == 'NaN':
				splitLine[s2Index] = int(splitLine[s1Index])
				
			if splitLine[e2Index] == 'NaN':
				splitLine[e2Index] = int(splitLine[e1Index])
		else:
			if splitLine[chr2Index] == 'NaN':
				continue # This line does not have correct chromosome 2 information (should we be skipping it?)
		
		s1 = int(splitLine[s1Index])
		e1 = int(float(splitLine[e1Index]))
		s2 = int(splitLine[s2Index])
		e2 = int(float(splitLine[e2Index]))
		chr2 = splitLine[chr2Index]
		
		chr1 = splitLine[chr1Index]
		o1 = splitLine[o1Index]
		o2 = splitLine[o2Index]
		
		#Make sure to switch the positions here as well
		#Some positions are swapped
		if int(e2) < int(e1):
			tmpE1 = e1
			e1 = e2
			e2 = tmpE1
			tmpS1 = s1
			s1 = s2
			s2 = tmpS1
		
		#Sometimes only the end is swapped.
		if int(e2) < int(s2):
			tmpS2 = s2
			s2 = e2
			e2 = tmpS2
			
		if int(e1) < int(s1):
			tmpS1 = s1
			s1 = e1
			e1 = tmpS1
		sampleNameIndex = header.index("sample_name")
		sampleName = splitLine[sampleNameIndex]

		if list(chr1)[0] == "c": #add the chr notation only when it is not already there
			
			#Check if this SV needs to be exluded
			svStr = chr1 + "_" + str(s1) + "_" + str(e1) + "_" + chr2 + "_" + str(s2) + "_" + str(e2) + "_" + sampleName
			
			svInfo[svStr] = [chr1, s1, e1, o1, chr2, s2, e2, o2, svType, sampleName]
			
		else:
			
			#Check if this SV needs to be exluded
			svStr = 'chr' + chr1 + "_" + str(s1) + "_" + str(e1) + "_" + 'chr' + chr2 + "_" + str(s2) + "_" + str(e2) + "_" + sampleName
			
			svInfo[svStr] = ['chr' + chr1, s1, e1, o1, 'chr' + chr2, s2, e2, o2, svType, sampleName]

#Get the SVs in CTCF as well and report these
cnvFile = sys.argv[7]
ctcfDels = dict()
geneLocations = [['chr16', 67596310, 67673088]] #location of CTCF
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
		
		#Check if there is a change at the location containing CTCF
		if 'chr' + splitLine[1] == geneLocations[0][0]:
			
			if int(splitLine[2]) <= geneLocations[0][2] and int(splitLine[3]) >= geneLocations[0][1]: #CTCF overlaps this segment
				if float(splitLine[5]) < 0: #assume that this is a deletion
					ctcfDels['CTCF_' + shortSampleName] = ['chr' + splitLine[1], splitLine[2], splitLine[3], splitLine[5]]

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
	
	#get info about the mutation in CTCF, this depends on if it has an SNV or deletion
	ctcfPair = 'CTCF_' + shortSampleName
	ctcfSNVPos = ''
	ctcfSNVType = ''
	ctcfSNVRef = ''
	ctcfSNVAlt = ''
	ctcfSVPos = ''
	ctcfSVSegmentVal = ''
	if ctcfPair in mutationInfo:
		ctcfSNVPos = 'chr' + mutationInfo[ctcfPair][0]
		ctcfSNVType = mutationInfo[ctcfPair][1]
		ctcfSNVRef = mutationInfo[ctcfPair][2]
		ctcfSNVAlt = mutationInfo[ctcfPair][3]
		
	if ctcfPair in ctcfDels:
		ctcfSVPos = ctcfDels[ctcfPair][0] + ':' + ctcfDels[ctcfPair][1] + '-' + ctcfDels[ctcfPair][2]
		ctcfSVSegmentVal = ctcfDels[ctcfPair][3]
	
	#Check if there is an SNV for this gene in any sample
	#This needs the original identifier
	shortSampleName = shortSampleNameLookup[splitPair[1]]
	
	shortPair = gene + '_' + shortSampleName
	snv = 'False'
	snvSignificant = 'False'
	snvSignificance = ''
	snvPos = ''
	snvClassification = ''
	snvRef = ''
	snvAlt = ''
	if gene + '_' + shortSampleName in degsWithSNVEvidence:
		snv = 'True'
		snvPos = 'chr' + mutationInfo[shortPair][0]
		snvClassification = mutationInfo[shortPair][1]
		snvRef = mutationInfo[shortPair][2]
		snvAlt = mutationInfo[shortPair][3]
		if gene + '_' + shortSampleName in snvDEGPairs[:,0]:
			if np.sign(zScoresSnv[gene + '_' + shortSampleName]) == 1:
				snvSignificant = 'True'
				snvSignificance = snvDEGPairs[snvDEGPairs[:,0] == gene + '_' + shortSampleName,1][0]
				
	
	#Repeat for SVs
	sv = 'False'
	svSignificant = 'False'
	svSignificance = ''
	svPosChr1 = ''
	svPosChr2 = ''
	svType = ''
	svO1 = ''
	svO2 = ''
	if gene + '_' + shortSampleName in degsWithSVEvidence:
		sv = 'True'
		
		#depending on the data type, the format of the SVs will be different.
		
		if re.search("BRCA", sys.argv[1], re.IGNORECASE):
			splitSampleName = splitPair[1].split("-")
			sampleCode = splitSampleName[2]
			
			convertedSampleName = 'brca' + sampleCode
			svGenePair = codingPairsShortSampleNames[gene + '_' + convertedSampleName]
			svData = "_".join(svGenePair.split("_")[1:])
			svPosChr1 = svInfo[svData][0] + ':' + str(svInfo[svData][1]) + '-' + str(svInfo[svData][2])
			svPosChr2 = svInfo[svData][4] + ':' + str(svInfo[svData][5]) + '-' + str(svInfo[svData][6])
			svO1 = svInfo[svData][3]
			svO2 = svInfo[svData][7]
			svType = svInfo[svData][8]
			
		if gene + '_' + shortSampleName in svDEGPairs[:,0]:
			if np.sign(zScoresSv[gene + '_' + shortSampleName]) == 1:
				svSignificant = 'True'
				svSignificance = svDEGPairs[svDEGPairs[:,0] == gene + '_' + shortSampleName,1][0]
	
	#Collect all information about this pair
	pairInfo.append([splitPair[0], splitPair[1], pair[1], cosmic, noOfSamples, ctcfSNVPos, ctcfSNVType, ctcfSNVRef, ctcfSNVAlt,
					 ctcfSVPos, ctcfSVSegmentVal, snv, 
					snvPos, snvClassification, snvRef, snvAlt, snvSignificant, snvSignificance,
					sv,svPosChr1, svPosChr2, svType, svO1, svO2, svSignificant, svSignificance])	

pairInfo = np.array(pairInfo, dtype='object')
header = 'Gene\tSample\tP-value_mutated_CTCF_vs_non-mutated_CTCF\tCOSMIC\tNumber_of_samples_in_which_DEG\tCTCF_SNV_pos\tCTCF_SNV_type\tCTCF_SNV_ref\tCTCF_SNV_alt\tCTCF_del_pos\tCTCF_del_segment_value\tSNV?\tSNV_pos\tSNV_classification\tSNV_ref\tSNV_alt\tSNV_significant?\tSNV_p-value\tSV?\tSV_pos_chr1\tSV_pos_chr2\tSV_type\tSV_orientation_chr1\tSV_orientation_chr2\tSV_significant?\tSV_significance'
np.savetxt(sys.argv[8], pairInfo, header=header, fmt='%s', delimiter='\t')

