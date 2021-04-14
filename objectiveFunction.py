#Ashleigh Thomas, contributor
#Julie Chow, contributor
#MS Thesis project for Ashleigh Thomas: Identifying Critical Cell Types for Neurodevelopmental Disorders (NDDs)

#imports
import numpy as np
import pandas as pd
import sys
import random
from random import randint
import scipy.stats as stats
#function that accepts:
#1: a pandas DataFrame (df) that represents UMI counts per barcoded cell per gene from a single cell RNA-seq dataset.
#2: geneSum: an float representing the coexpression of genes, currently not being used.
#3: moduleGenes: a list of strings contaning the names of the genes in the module of interest.
#4: comparison: a flag that if == 1 returns data that will be used to make a heatmap of Z-scored per cell type, averaged per cell type UMI data.
def calcNormSum(df, geneSum, moduleGenes,comparison):
	m = df.shape[0]#df dimensions
	n = df.shape[1]
	if geneSum == 0:#if not using coexpression data, let it be 1 (since it will be multiplied in, so this will have no effect)
		geneSum = 1
	rankDF = df.rank(axis=0,method='dense')#ranking column
	#filter rankDF to only have moduleGenes
	rankDF = rankDF[rankDF.index.isin(moduleGenes)]
	if comparison == 1:
		return rankDF
	rankDF['sum'] = rankDF.sum(axis=1)#sum of ranks of human genes in genome
	rankDF['normSum'] = rankDF['sum'] / (n-1)#summed ranks of UMI counts of all human genes divided by number of cells - 1
	rankDF['normLenModGenesSum'] = rankDF['normSum'].pow(1/len(moduleGenes))#put the summed ranks of UMI counts/#cells-1 to the power of (1/number of module genes)
	normSumDF = rankDF['normLenModGenesSum']#grabbing dataframe from only the column that has been exponentiated
	resultDict = normSumDF.to_dict()#turning said dataframe into a dictionary to return
	return resultDict#dictionary mapping from gene to value

def leftHandSide(genes, umiDF, co_exp, comparison):
	co_exp_dictionary = {}
	for gene in reversed(genes):
		idx = genes.index(gene)
		spec_co_exp = co_exp[(co_exp['geneA'] == gene) | (co_exp['geneB'] == gene)]
		growing_co_exp = 0  # growing summation
		for chain_idx in np.arange(0, idx):
			gene_a = genes[idx]
			gene_b = genes[chain_idx]
			try_c_1 = spec_co_exp[(spec_co_exp['geneA'] == gene_a) & (spec_co_exp['geneB'] == gene_b)]
			try_c_2 = spec_co_exp[(spec_co_exp['geneB'] == gene_a) & (spec_co_exp['geneA'] == gene_b)]
			if try_c_1.empty and try_c_2.empty:  # Then the co-expression doesn't exist, set to 0
				c_value = 0
			elif try_c_2.empty:
				c_value = try_c_1.iloc[0]['coexp']
			else:
				c_value = try_c_2.iloc[0]['coexp']

			growing_co_exp += c_value
		if growing_co_exp == 0:  # if there are no co-expression values associated, just set to 1. will do nothing
			growing_co_exp = 1
		co_exp_dictionary[gene] = growing_co_exp
	geneSumDF = pd.DataFrame.from_dict(co_exp_dictionary, orient='index')
	geneSumDF = geneSumDF.dropna()
	geneSum = geneSumDF[0].sum(axis=0)
	geneSumDF.columns=['value']
	#now to use coexpression and take into account all ~24k genes into the ranking
	lhsDict = calcNormSum(umiDF, geneSum, genes,comparison)#this is multiplying coexpression sum of module genes by umi counts for all genes
	return lhsDict, geneSum

def rightHandSide(distMatrix, cellTypeInterest):
	cellMeta = open('/share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/cell_metadata.csv')
	typeBarcodeDict = {}
	cellMeta.readline()#header
	for line in cellMeta:
		line = line.strip()
		line = line.split(',')
		barcode = line[0]
		cellType = line[1]
		if cellType != cellTypeInterest:
		#	continue
			if line[1] not in typeBarcodeDict:
				typeBarcodeDict[cellType] = [barcode]
			else:
				typeBarcodeDict[cellType].append(barcode)
		else:
			continue
#		if line[1] not in typeBarcodeDict:
#			typeBarcodeDict[cellType] = [barcode]
#		else:
#			typeBarcodeDict[cellType].append(barcode)
#	print(distMatrix)
#	exit()
	dfArr = distMatrix.to_numpy()
	print(dfArr)
	n = dfArr.shape[0]
#	print(n)#should be 9995 for ExN
#	n = dfArr.shape[0]-1
#	triArr = dfArr[np.triu_indices(n, k=1)]
	r,c = np.triu_indices(n,k=1)#TODO why is r or c 9995 when size of array is 9995 (i.e. last index should be 9994)? it's because now the size has changed
	print(r,c)
	triArr = dfArr[r,c]
	np.savetxt(cellTypeInterest + 'UpperTriangle.txt', triArr, delimiter='\t')
#	avg = triArr[np.nonzero(triArr)].mean()
	avg = np.median(triArr[np.nonzero(triArr)])
#	avg = np.mean(triArr)#we did this before, trying median instead to see if it'll work better for skewed 0 matrices
#	avg = np.median(triArr)#tried this, end up with -inf because too many lines have so many 0s
	return avg
def isUnique(series):
	array = series.to_numpy()
	return (array[0] == array).all()

#def multiplication(lhsDict, rhsVal, geneSum, comparison, fileModifier, cellType):#lhsDict: gene to lhs value dictionary; rhsVal: rhs value for this cell type
#TODO put back in geneSum arg if want coexpression in here
def multiplication(lhsDict, rhsVal, comparison, fileModifier, cellType):#lhsDict: gene to lhs value dictionary; rhsVal: rhs value for this cell type
	result = rhsVal
	print('rhs val',rhsVal)
	valsFile = open(fileModifier + cellType + 'LHSAndRHSValues.txt', 'w')
	valsFile.write('RightHandSide' + '\t' + str(rhsVal) + '\n')
	result = float(result)
#	outfile = open('error.out', 'w')
#	result = result * geneSum#TODO uncomment this after mini testing
	if comparison == 1:#lhsDict is actually a genes by cells DF
		result = result * lhsDict
		dfStack =result.stack()
	#	print(result.std())
		#result.std().to_csv('error.out',sep='\t')
		stdDevIsZero = isUnique(result)
		print(stdDevIsZero)
	#	print(result)
	#	print(dfStack)
	#	exit()
#		if result.std() == 0:
#			df = (result - dfStack.mean()) / 1
#		else:
		dfZ = (result - dfStack.mean())/dfStack.std()
#		dfZ = (result - dfStack.mean())/result.std()#dividing by zero somewhere
		dfZ = dfZ.replace([np.inf, -np.inf], 0)
		print(dfZ)
		print('hi')
		return dfZ,1#do more here, does mul/multiply work? also do z score
	lhsVal = 0
	first = True
	for key in lhsDict:
		if first:
			lhsVal = lhsDict[key]
			first = False
		result = result * float(lhsDict[key])
		lhsVal = lhsVal * lhsDict[key]
	valsFile.write('LeftHandSide' + '\t' + str(lhsVal))
	valsFile.close()
	logRes = np.log10(result)
	return logRes, result


def getModuleGenes():
	moduleGenes = []
	moduleGenesFile = open(sys.argv[1], 'r')
	for gene in moduleGenesFile:
		gene = gene.strip()
		moduleGenes.append(gene)
	return moduleGenes

def getCoexpressionDataChunks(moduleGenes, fileModifier):
	#this next segment should be uncommented if you want to use actual coexpression data
#	chunkSize = 10 ** 6
#	fileName = 'coexpressionDFNoMaxRank' + fileModifier + '.txt'
#	for chunk in pd.read_csv(sys.argv[2], sep='\t', chunksize = chunkSize, header=None):
#		chunk.columns = ['geneA', 'geneB', 'coexp']
#		chunk = chunk[(chunk['geneA'].isin(moduleGenes)) & (chunk['geneB'].isin(moduleGenes))]
#		chunk.to_csv(fileName, mode='a',index=False,header=False,sep='\t')
#	coexpDF = pd.read_csv(fileName, sep='\t', names=['geneA', 'geneB', 'coexp'])
	#end of actual coexpression data section
	coexpDF = pd.DataFrame([['one', 'two', 5]], columns=['geneA', 'geneB', 'coexp'])
	return coexpDF

def getUMIFiles():
	umiFile = open(sys.argv[3])
	location = umiFile.readline()
	location = location.strip()
	umiFiles = []
	for line in umiFile:
		line = line.strip()
		umiFiles.append(location + line)
	return umiFiles
		

def getUMICount(umiInfile):
#	umiInfile = sys.argv[3]
	umiDF = pd.read_csv(umiInfile, sep='\t')
#	umiDF = pd.read_csv(umiInfile, sep=' ')#TODO after random barcodes experiment change this back to sep = '\t'
	umiDF.drop(umiDF.columns[len(umiDF.columns)-1], axis=1, inplace=True)
	umiDF = umiDF.dropna()
	return umiDF

def getDistanceMatrix(barcodes):
	wholeDF = pd.read_csv(sys.argv[4], usecols=barcodes,sep=',')#'/share/hormozdiarilab/Data/scRNASeq_AllenBrain/RNN_snn.csv')#before wasn't sep=',' didn't specify
#	print(wholeDF.shape)
#	exit()
	subDF = wholeDF[~wholeDF.index.isin(barcodes)]#get only rows that are barcodes we want#TODO delete the ~ later, this is grabbing all rows except the barcodes for this cell type
	subDFWithin = wholeDF[wholeDF.index.isin(barcodes)]#get only rows that are barcodes we want#TODO delete the ~ later, this is grabbing all rows except the barcodes for this cell type
	return subDF, subDFWithin

def getPrimeUniverse(umiFilesList):
	cellTypeToPrimeUniverseBarcodesDict = {}
	cellTypes = []
	cellTypesToBarcodesDict = {}
	for entry in umiFilesList:
		infile = open(entry)
		fileTerms = entry.split('/')
		cellType = fileTerms[-1]
		cellType = cellType[0:-13]#cuts off 'UMICounts.txt' from the name to get ExN and Mic as cell types
		cellTypes.append(cellType)
		barcodes = infile.readline()
		barcodes = barcodes.strip()
		barcodes = barcodes.split('\t')
		cellTypesToBarcodesDict[cellType] = barcodes
	print('length of ExN barcodes:',len(cellTypesToBarcodesDict['ExN']))#the prime universe of ExN should have length of Mic barcodes
	print('length of Mic barcodes:',len(cellTypesToBarcodesDict['Mic']))#the prime universe of Mic should have length of ExN barcodes
	for ct in cellTypes:
		for entry in cellTypes:
			if entry == ct:
				continue
			if ct in cellTypeToPrimeUniverseBarcodesDict:
				for barcode in cellTypesToBarcodesDict[entry]:
					cellTypeToPrimeUniverseBarcodesDict.append(barcode)
			else:
				cellTypeToPrimeUniverseBarcodesDict[ct] = cellTypesToBarcodesDict[entry]#if ct == ExN, and entry == Mic, make the prime universe of all but ExN have the contents of Mic's barcodes (and vice versa)	
	print('length of ExN prime universe barcodes:', len(cellTypeToPrimeUniverseBarcodesDict['ExN']))
	print('length of Mic prime universe barcodes:', len(cellTypeToPrimeUniverseBarcodesDict['Mic']))
	#exit()
	return cellTypeToPrimeUniverseBarcodesDict
def chooseRandomBarcodes(n, infilePath):
	umiFile = open(infilePath)
	header = umiFile.readline()
	header = header.strip()
	header = header.split(' ')
	umiFile.close()
	maxVal = len(header)-1
	print(maxVal)
	count = 0
	barcodesDict = {}
	while count < n:
		randNum = random.randint(0,maxVal)
		barcodeChosen = header[randNum][1:-1]
	#	if randNum in barcodesDict:
		if barcodeChosen in barcodesDict:
			continue
		else:
	#		barcodesDict[randNum] = 1
			barcodesDict[barcodeChosen] = 1
			count += 1
	return list(barcodesDict.keys())

def rightHandSidePValue(btwnDF, withinDF):
	np.fill_diagonal(withinDF.values, 0)#NOTE this is for when the matrix is square (aka when we are getting only barcodes for this cell type, uncomment this later)
	btwnMaxDF = pd.DataFrame()
	btwnMaxDF["max similarity"] = btwnDF.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
	withinMaxDF = pd.DataFrame()
	withinMaxDF["max similarity"] = withinDF.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
	val = stats.ttest_rel(withinMaxDF['max similarity'], btwnMaxDF['max similarity'], alternative='greater')
	result = (1-val.pvalue)
	return result

#function calls:
moduleGenes = getModuleGenes()
fileModifier = sys.argv[6]
comparison = float(sys.argv[7])
coexpDF = getCoexpressionDataChunks(moduleGenes, fileModifier)#TODO uncomment this for later
umiFilesList = getUMIFiles()#TODO uncomment when done with random barcodes
outfile = open(sys.argv[5], 'w')
outfile.write('Gene\tLog\tRaw\n')
for cellTypeFile in umiFilesList:
	line = cellTypeFile.split('/')
	cellType = line[-1]
	cellType = cellType[:-13]#remove the UMICounts.txt from the end (i.e. EndUMICounts.txt becomes cell type of End)
	print('cell type:', cellType)
	umiDF = getUMICount(cellTypeFile)
	barcodes = umiDF.columns
	btwnMatrix, withinMatrix = getDistanceMatrix(barcodes)#note withinMatrix is the prime distance matrix, btwnMatrix is the matrix of comparisons amongst the same cell type
	rhsVal = rightHandSidePValue(btwnMatrix, withinMatrix)
	lhsDict, geneSum = leftHandSide(moduleGenes, umiDF, coexpDF,comparison)
	objectiveFuncValueLog, rawObjectiveFunctionValue = multiplication(lhsDict, rhsVal, comparison, fileModifier, cellType)#without coexpression
	if comparison != 1:
		outfile.write(cellType + '\t' + str(objectiveFuncValueLog) + '\t' + str(rawObjectiveFunctionValue) + '\n')
	else:
		objectiveFuncValueLog.to_csv(fileModifier +cellType+ '.txt', sep='\t')
outfile.close()
#example of how to run for End cell type
#python3 objectiveFunctionChainRuleAllCellTypesNonMaxRank_max_sim.py m2ExtendedASDwithID.txt polioEFCoexpressionDataWithoutNans.txt /share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/umiCounts/EndUMICounts.txt /share/hormozdiarilab/Data/scRNASeq_AllenBrain/RNN_snn.csv End
#python3 objectiveFunctionChainRuleAllCellTypesNonMaxRank_max_sim.py m2ExtendedASDwithID.txt polioEFCoexpressionDataNoNanNoZero.txt /share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/umiCounts/EndUMICounts.txt /share/hormozdiarilab/Data/scRNASeq_AllenBrain/RNN_snn.csv End
#parameters:
#1: file containing 1 gene per line, represents all of the genes in the module
#2: coexpression file without -nan lines. This contains more than the info we need as it contains the coexpression data for all genes in the dataset, not just the module genes.
#3: file containing: 1st line: location of the UMI count files. All following lines are the names of the UMICount files.
#4: distance matrix outputted from Seurat
#5 output file name
#6 file modifier so temporary files don't get contaminated
#python3 objectiveFunctionChainRuleAllCellTypesNonMaxRank_max_sim.py m2ExtendedASDwithID.txt polioEFCoexpressionDataNoNanNoZero.txt umiInputFile.txt /share/hormozdiarilab/Data/scRNASeq_AllenBrain/RNN_snn.csv randomM2ExtNoMaxRank.txt
#python3 objectiveFunctionChainRuleAllCellTypesNonMaxRank_max_sim.py randomUnconstrainedGenesM2.txt polioEFCoexpressionDataNoNanNoZero.txt umiInputFile.txt /share/hormozdiarilab/Data/scRNASeq_AllenBrain/RNN_snn.csv randomM2ExtNoMaxRank.txt/
