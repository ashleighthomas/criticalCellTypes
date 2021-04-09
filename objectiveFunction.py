# python3 umiCountsRankCoexp.py  \
# /share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/umiCounts/EndUMICounts.txt \
# compactPolioEFCoexp.txt \
# miniFolder/ \
# EndCoexp.txt
#this one at first will use original objective function, later will use chain rule
import numpy as np
import pandas as pd
import sys
import random
from random import randint
import scipy.stats as stats
def calcNormSum(df, geneSum, moduleGenes,comparison):#,barcode):
	m = df.shape[0]
	n = df.shape[1]
	if geneSum == 0:
		geneSum = 1
	rankDF = df.rank(axis=0,method='dense')#ranking column
	#filter rankDF to only have moduleGenes
	rankDF = rankDF[rankDF.index.isin(moduleGenes)]
#	print(rankDF)
#	rankDF.to_csv('M2MicRank.txt', sep='\t')
	if comparison == 1:
#		if barcode in rankDF.columns:
#			rankDF[barcode].to_csv('barcode.txt',sep='\t')
		return rankDF
#		resultDict = rankDF.to_dict()
#		return resultDict
	rankDF['sum'] = rankDF.sum(axis=1)#sum of ranks
	rankDF['normSum'] = rankDF['sum'] / (n-1)#*(maxRank-1))#removing *(maxRank-1) removes the problem of small cell types getting larger OF values
	#normSum being rankDF['sum']/(n-1) was what the successful run was that was reported in file M1NoCoexpressionOutput021821.txt and M2NoCoexpressionOutput021821.txt. versions with coexpression (single coexpression value not made for coexpression of only module genes) can be found at: M1CoexpressionOutput021921.txt and M2CoexpressionOutput021921.txt (also in powerpoint LabUpdate022221.ppt)
	rankDF['normLenModGenesSum'] = rankDF['normSum'].pow(1/len(moduleGenes))#chain rule!
	#rankDF['normSum'] = rankDF['normSum'] * geneSum#this used to be nested in if useAllGenes#commenting out 02/16/21
	normSumDF = rankDF['normLenModGenesSum']#used to be 'normSum'
	resultDict = normSumDF.to_dict()
	return resultDict

def calcNormSumNoCoexpression(df, moduleGenes, comparison):#,barcode):
	m = df.shape[0]
	n = df.shape[1]
	rankDF = df.rank(axis=0,method='dense')#ranking column
	#filter rankDF to only have moduleGenes
	rankDF = rankDF[rankDF.index.isin(moduleGenes)]
	if comparison == 1:
#		if barcode in rankDF.columns:
#			rankDF[barcode].to_csv('barcode.txt',sep='\t')
		return rankDF
#		resultDict = rankDF.to_dict()
#		return resultDict
	rankDF['sum'] = rankDF.sum(axis=1)#sum of ranks
	rankDF['normSum'] = rankDF['sum'] / (n-1)#*(maxRank-1))#removing *(maxRank-1) removes the problem of small cell types getting larger OF values
	#normSum being rankDF['sum']/(n-1) was what the successful run was that was reported in file M1NoCoexpressionOutput021821.txt and M2NoCoexpressionOutput021821.txt. versions with coexpression (single coexpression value not made for coexpression of only module genes) can be found at: M1CoexpressionOutput021921.txt and M2CoexpressionOutput021921.txt (also in powerpoint LabUpdate022221.ppt)
	rankDF['normLenModGenesSum'] = rankDF['normSum'].pow(1/len(moduleGenes))#chain rule!
	#rankDF['normSum'] = rankDF['normSum'] * geneSum#this used to be nested in if useAllGenes#commenting out 02/16/21
	normSumDF = rankDF['normLenModGenesSum']#used to be 'normSum'
	resultDict = normSumDF.to_dict()
	return resultDict

def leftHandSide(genes, umiDF, co_exp, comparison):#,barcode):	
#def leftHandSide(genes, umiDF, comparison):#,barcode):	
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
	lhsDict = calcNormSum(umiDF, geneSum, genes,comparison)#,barcode)#this is multiplying coexpression sum of module genes by umi counts for all genes
	return lhsDict, geneSum

def leftHandSideWithoutCoexpression(genes, umiDF, comparison):#,barcode):	
#def leftHandSide(genes, umiDF, comparison):#,barcode):	
	'''
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
	'''
	#now to use coexpression and take into account all ~24k genes into the ranking
#	lhsDict = calcNormSum(umiDF, geneSum, genes,comparison)#,barcode)#this is multiplying coexpression sum of module genes by umi counts for all genes
	lhsDict = calcNormSumNoCoexpression(umiDF, genes, comparison)
	return lhsDict#, geneSum
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
print('got module genes, about to get coexpression data for module genes')
fileModifier = sys.argv[6]
comparison = float(sys.argv[7])
coexpDF = getCoexpressionDataChunks(moduleGenes, fileModifier)#TODO uncomment this for later
print('got coexpression data')
umiFilesList = getUMIFiles()#TODO uncomment when done with random barcodes
#umiFilePath = sys.argv[3]
#barcodes = chooseRandomBarcodes(9995, 'allBarcodes.txt')
#print(len(barcodes))
#umiFile = open(sys.argv[3])
#umiFilesList = []
#umiFilesList.append(umiFile)
#print(barcodes)
#print(umiFilesList)#['/share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/umiCounts/ExNUMICounts.txt', '/share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/umiCounts/MicUMICounts.txt']
#exit()
outfile = open(sys.argv[5], 'w')
outfile.write('Gene\tLog\tRaw\n')
#outfile = open('outputRandomM2Ext_021221.txt', 'w')
barcode = 'ACTGGCACCGAG'
#primeUniversePerCellTypeDict = getPrimeUniverse(umiFilesList)
#print(primeUniversePerCellTypeDict['ExN'])
#exit()
cellType = 'allCellsRandomSubset'
#print(umiFilePath)
#exit()
#umiDF = getUMICount(umiFilePath)
#print(umiDF)
'''
print('about to call getDistMatrix,***')
btwnMatrix, withinMatrix = getDistanceMatrix(barcodes)#note, this is currently the prime distance matrix
#exit()
print(btwnMatrix.shape)
print(withinMatrix.shape)
#	np.fill_diagonal(distMatrix.values, 0)#NOTE this is for when the matrix is square (aka when we are getting only barcodes for this cell type, uncomment this later)
np.fill_diagonal(withinMatrix.values, 0)#NOTE this is for when the matrix is square (aka when we are getting only barcodes for this cell type, uncomment this later)
df = pd.DataFrame()
df["max_sim_per_cell"] = btwnMatrix.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
df['max_sim_per_cell'].to_csv("%s_max_sim_per_cell_with_cells_not_in_same_cell_type.txt" % cellType)#val per cell
summaryVal = df['max_sim_per_cell'].sum() / df.shape[0]#this now represents the average similarity with cells outside of the cell type per cell type
sumFile = open(cellType + 'BetweenSummaryValue.txt', 'w')
sumFile.write(str(summaryVal))
sumFile.close()
df = pd.DataFrame()
df["max_sim_per_cell"] = withinMatrix.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
df['max_sim_per_cell'].to_csv("%s_max_sim_per_cell_within_barcode_group.txt" % cellType)#val per cell
summaryVal = df['max_sim_per_cell'].sum() / df.shape[0]#this now represents the average similarity with cells outside of the cell type per cell type
sumFile = open(cellType + 'WithinSummaryValue.txt', 'w')
sumFile.write(str(summaryVal))
sumFile.close()
exit()
'''
#rhsOutfile = open('rightHandSidePValuesResults033021.txt', 'w')
for cellTypeFile in umiFilesList:
	line = cellTypeFile.split('/')#TODO uncomment these next 3 lines after done with random barcodese xperiment
	cellType = line[-1]
	cellType = cellType[:-13]#remove the UMICounts.txt from the end (i.e. EndUMICounts.txt becomes cell type of End)
#	cellType = 'allCellsRandomSubset'#TODO remove this after the random barcodes experiment is done
	print('cell type:', cellType)
	print('about to open file', cellTypeFile)
	umiDF = getUMICount(cellTypeFile)
	print(umiDF)
	barcodes = umiDF.columns#TODO uncomment this after random barcodes experiment done
	print('got UMI count')
	#distMatrix = getDistanceMatrix(barcodes)#note, this is currently the prime distance matrix
	btwnMatrix, withinMatrix = getDistanceMatrix(barcodes)#note, this is currently the prime distance matrix
	rhsVal = rightHandSidePValue(btwnMatrix, withinMatrix)
#	rhsOutfile.write(cellType + '\t' + str(rhsVal) + '\n')
	#rhsOutfile.close()
	#exit()
#	np.fill_diagonal(distMatrix.values, 0)#NOTE this is for when the matrix is square (aka when we are getting only barcodes for this cell type, uncomment this later)
	#this next section of code was implemented, now commented out to see what rightHandSidePValue returns (does same thing except 1-pval)
#	np.fill_diagonal(withinMatrix.values, 0)#NOTE this is for when the matrix is square (aka when we are getting only barcodes for this cell type, uncomment this later)
#	df = pd.DataFrame()
#	df["max_sim_per_cell"] = btwnMatrix.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
#	df['max_sim_per_cell'].to_csv("%s_max_sim_per_cell_with_cells_not_in_same_cell_type032921.txt" % cellType)#val per cell
#	summaryVal = df['max_sim_per_cell'].sum() / df.shape[0]#this now represents the average similarity with cells outside of the cell type per cell type
#	sumFile = open(cellType + 'BetweenSummaryValue032921.txt', 'w')
#	sumFile.write(str(summaryVal))
#	sumFile.close()

#	df = pd.DataFrame()
#	df["max_sim_per_cell"] = withinMatrix.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col#TODO change to distMatrix.max
#	df['max_sim_per_cell'].to_csv("%s_max_sim_per_cell_within_barcode_group032921.txt" % cellType)#val per cell
#	summaryVal = df['max_sim_per_cell'].sum() / df.shape[0]#this now represents the average similarity with cells outside of the cell type per cell #type
#	sumFile = open(cellType + 'WithinSummaryValue032921.txt', 'w')
#	sumFile.write(str(summaryVal))
#	sumFile.close()
#	exit()
	#distMatrix['max_sim_per_cell'].to_csv("%s_max_sim_per_cell.txt" % fileModifier)
#	exit()
	print('got distance matrix')
	print('about to call left hand side func')
	lhsDict, geneSum = leftHandSide(moduleGenes, umiDF, coexpDF,comparison)#,barcode)
#	lhsDict = leftHandSideWithoutCoexpression(moduleGenes, umiDF, comparison)# coexpDF,comparison)#,barcode)
	print('got lhs values')
	#rhsVal = rightHandSide(distMatrix, cellType)#TODO uncomment once I'm done running the prime universe calculations
	print('got rhs values')
	print('---------------------------------')
	#objectiveFuncValueLog, rawObjectiveFunctionValue = multiplication(lhsDict, rhsVal, geneSum,comparison, fileModifier, cellType)#with coexpression
	objectiveFuncValueLog, rawObjectiveFunctionValue = multiplication(lhsDict, rhsVal, comparison, fileModifier, cellType)#without coexpression
#	print('finished multiplication, now writing to output file, final value was', objectiveFuncValueLog, rawObjectiveFunctionValue)
	if comparison != 1:
		outfile.write(cellType + '\t' + str(objectiveFuncValueLog) + '\t' + str(rawObjectiveFunctionValue) + '\n')
	else:
#		dfZScored = objectiveFuncValueLog.mean(axis=1)
		#dfZScored.to_csv(fileModifier +cellType+ '.txt', sep='\t', header=False)
		objectiveFuncValueLog.to_csv(fileModifier +cellType+ '.txt', sep='\t')
		
		
outfile.close()
#rhsOutfile.close()
#outfile.close()
#rhsOutfile.close()
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
