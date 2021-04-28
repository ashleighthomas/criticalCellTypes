#Ashleigh Thomas, contributor
#Julie Chow, contributor
#MS Thesis project for Ashleigh Thomas: Identifying Critical Cell Types for Neurodevelopmental Disorders (NDDs)

#imports
import numpy as np
import pandas as pd
import sys
import argparse
import os
import random
from random import randint
import scipy.stats as stats
import argparse
#function that accepts:
#1: a pandas DataFrame (df) that represents UMI counts per barcoded cell per gene from a single cell RNA-seq dataset.
#2: geneSum: an float representing the coexpression of genes, currently not being used.
#3: moduleGenes: a list of strings contaning the names of the genes in the module of interest.
#4: comparison: a flag that if == 1 returns data that will be used to make a heatmap of Z-scored per cell type, averaged per cell type UMI data.
#this function returns a dictionary mapping from each module gene to the exponentiated value
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

#This function calculates the left hand side of the objective function (the chain rule implementation as described in the thesis), including calling calcNormSum above.
#The arguments are as follows:
#genes: a list of module genes
#umiDF: a pandas dataframe containing the UMI counts for the current cell type
#co_exp: a pandas dataframe with the gene coexpression data
#comparison: a flag that (if it equals 1) will be utilized for making a heatmap. This is passed to calcNormSum above.
#This function returns a dictionary from module genes to left hand side values of the objective function, as well as the coexpression value for the current cell type.
def leftHandSide(genes, umiDF, co_exp, comparison):
	co_exp_dictionary = {}
	#implements chain rule for Bayes' Theorem for three events. Forms geneSumDF which is used for coexpression values
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
		co_exp_dictionary[gene] = growing_co_exp#for each gene in module genes, add an entry for that gene of the growing coexpression score
	geneSumDF = pd.DataFrame.from_dict(co_exp_dictionary, orient='index')
	geneSumDF = geneSumDF.dropna()
	geneSum = geneSumDF[0].sum(axis=0)
	geneSumDF.columns=['value']
	#now to use coexpression and take into account all ~24k genes into the ranking
	lhsDict = calcNormSum(umiDF, geneSum, genes,comparison)#this is multiplying coexpression sum of module genes by umi counts for all genes
	return lhsDict, geneSum

#This function calculates the average of the right hand side upper triangle of the cell to cell distance matrix generated by Seurat.
#The arguments are:
#distMatrix: Seurat distance matrix
#cellTypeInterest: current cell type
#This function returns the average of the upper right hand side triangle (cell to cell distance matrix) from Seurat.
#This is not currently being used, as we have transitioned to rightHandSidePValue.
def rightHandSide(distMatrix, cellTypeInterest):
	#opening the Polioudakis metadata file
	cellMeta = open('/share/hormozdiarilab/Experiments/Ashleigh_workspace/nddTissue/forJulie/cell_metadata.csv')
	typeBarcodeDict = {}
	cellMeta.readline()#header
	#going through cellMeta and gathering barcodes into cell types.
	for line in cellMeta:
		line = line.strip()
		line = line.split(',')
		barcode = line[0]
		cellType = line[1]
		if cellType != cellTypeInterest:
			if line[1] not in typeBarcodeDict:
				typeBarcodeDict[cellType] = [barcode]
			else:
				typeBarcodeDict[cellType].append(barcode)
		else:
			continue
	dfArr = distMatrix.to_numpy()
	n = dfArr.shape[0]
	r,c = np.triu_indices(n,k=1)#gets upper right hand triangle without the diagonal values
	triArr = dfArr[r,c]
	avg = np.median(triArr[np.nonzero(triArr)])#getting average of the upper right hand triangle
	return avg

#This function takes in a series that the right hand side of the objective function multiplied by all of the left hand side values. This is involved in the multiplication function. This will return true or false
def isUnique(series):
	array = series.to_numpy()
	return (array[0] == array).all()

#This function multiplies the right and left hand sides of the objective function values to get the final value for the current cell type.
#lhsDict: in non comparison case, this is a gene to left hand value dictionary; in comparison case, this is a genes by cells DF
#rhsVal: this is a float representing the right hand side of the objective function
#comparison: same thing as the previous functions
#fileModifier: This is used for a temporary output file
#geneSum: this is for coexpression, this is not currently in use
#celltype: current cell type of interest
#This function returns 2 values: for comparison, it returns the z scored dataframe and 1 (unused); for non comparison: it returns the log transformed value and the raw value
#This next line is an alternative function header that takes coexpression values into the function
#def multiplication(lhsDict, rhsVal, geneSum, comparison, fileModifier, cellType):#lhsDict: gene to lhs value dictionary; rhsVal: rhs value for this cell type
def multiplication(lhsDict, rhsVal, comparison, fileModifier, cellType):#lhsDict: gene to lhs value dictionary; rhsVal: rhs value for this cell type
	result = rhsVal#taking the right hand side value of the objective function and making result equal to this
	valsFile = open(fileModifier + cellType + 'LHSAndRHSValues.txt', 'w')
	valsFile.write('RightHandSide' + '\t' + str(rhsVal) + '\n')
	result = float(result)
#	result = result * geneSum#uncomment this if you want to use geneSum which is coexpression value
	if comparison == 1:#lhsDict is actually a genes by cells DF in this case
		result = result * lhsDict#multiplying the right hand side by left hand side dataframe
		dfStack =result.stack()#converting to a stack
		stdDevIsZero = isUnique(result)
		dfZ = (result - dfStack.mean())/dfStack.std()
		dfZ = dfZ.replace([np.inf, -np.inf], 0)
		return dfZ,1#returns the z scored dataframe of right and left hand sides multiplied together, along with 1 as a result which does not get used in this comparison==1 (heatmap) case
	lhsVal = 0
	first = True
	#going through left hand side dictionary and multiplying the values together
	for key in lhsDict:
		if first:
			lhsVal = lhsDict[key]
			first = False
		result = result * float(lhsDict[key])
		lhsVal = lhsVal * lhsDict[key]
	valsFile.write('LeftHandSide' + '\t' + str(lhsVal))
	valsFile.close()
	logRes = np.log10(result)#log base 10 transforming result
	return logRes, result#return the log transformed result and the raw result

#This obtains module genes from a file, returns a list of module genes
def getModuleGenes(filepath):
	moduleGenes = []
	moduleGenesFile = open(filepath, 'r')
	for gene in moduleGenesFile:
		gene = gene.strip()
		moduleGenes.append(gene)
	return moduleGenes

#This reads in the coexpression for module genes, also takes in the file modifier for an output file name
def getCoexpressionDataChunks(moduleGenes, fileModifier,coexpEnabled):
	#this next segment should be uncommented if you want to use actual coexpression data
	if coexpEnabled:
		chunkSize = 10 ** 6
		fileName = 'coexpressionDFNoMaxRank' + fileModifier + '.txt'
		for chunk in pd.read_csv(sys.argv[2], sep='\t', chunksize = chunkSize, header=None):#TODO when coexpression is used again, use argparse to enforce this as a positional argument rather than an optional argument as it currently is, then pass in the coexp name instead of calling sys.argv[2]
			chunk.columns = ['geneA', 'geneB', 'coexp']
			chunk = chunk[(chunk['geneA'].isin(moduleGenes)) & (chunk['geneB'].isin(moduleGenes))]
			chunk.to_csv(fileName, mode='a',index=False,header=False,sep='\t')
		coexpDF = pd.read_csv(fileName, sep='\t', names=['geneA', 'geneB', 'coexp'])
		return coexpDF
	#end of actual coexpression data section
	else:
	#this is a section of fake mini coexpression dataframe so that errors don't happen. If you want to use real coexpression data, uncomment the previous section
		coexpDF = pd.DataFrame([['one', 'two', 5]], columns=['geneA', 'geneB', 'coexp'])
		return coexpDF

#This reads in a file where the first line is a header that is the filepath location to the files, and the subsequent lines are the names of the cell type UMI count files. It returns a list of the files
def getUMIFiles(umiInputFile):
	umiFile = open(umiInputFile)
	location = umiFile.readline()
	location = location.strip()
	umiFiles = []
	for line in umiFile:
		line = line.strip()
		umiFiles.append(location + line)
	return umiFiles
		
#This reads in a UMI count input file for a cell type and gets the UMI counts, and returns a pandas DataFrame with NAs dropped
def getUMICount(umiInfile):
	umiDF = pd.read_csv(umiInfile, sep='\t')
	umiDF.drop(umiDF.columns[len(umiDF.columns)-1], axis=1, inplace=True)
	umiDF = umiDF.dropna()
	return umiDF

#This function reads in the output of Seurat (distance matrix, the second input) and the list of barcodes for the current cell type
#This returns 2 dataframes. 
#The first is a dataframe that is grabbing every barcode that is not within the same cell type as the current cell type.
#The second is a dataframe that is only grabbing barcodes within the same cell type. 
def getDistanceMatrix(barcodes, distMatrixInfile):
	wholeDF = pd.read_csv(distMatrixInfile, usecols=barcodes,sep=',')
	subDF = wholeDF[~wholeDF.index.isin(barcodes)]#get only rows that are barcodes we want#TODO delete the ~ later, this is grabbing all rows except the barcodes for this cell type
	subDFWithin = wholeDF[wholeDF.index.isin(barcodes)]#get only rows that are barcodes we want#TODO delete the ~ later, this is grabbing all rows except the barcodes for this cell type
	return subDF, subDFWithin

#This is a function that chooses n random barcodes for an experiment.
def chooseRandomBarcodes(n, infilePath):
	umiFile = open(infilePath)
	header = umiFile.readline()
	header = header.strip()
	header = header.split(' ')
	umiFile.close()
	maxVal = len(header)-1
	count = 0
	barcodesDict = {}
	while count < n:
		randNum = random.randint(0,maxVal)
		barcodeChosen = header[randNum][1:-1]
		if barcodeChosen in barcodesDict:
			continue
		else:
			barcodesDict[barcodeChosen] = 1
			count += 1
	return list(barcodesDict.keys())

#This function is the final right hand side implementation involving p values. This returns for each cell type 1 - the likelihood of the within cell type average max similarity being greater than the between cell type average max similarity. For the 16 cell types from Polioudakis, this ends up being almost 1, as the p value is very close to 0.
def rightHandSidePValue(btwnDF, withinDF):
	np.fill_diagonal(withinDF.values, 0)#this is for when the matrix is square (aka when we are getting only barcodes for this cell type)
	btwnMaxDF = pd.DataFrame()
	btwnMaxDF["max similarity"] = btwnDF.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col
	withinMaxDF = pd.DataFrame()
	withinMaxDF["max similarity"] = withinDF.max(axis=0)#axis=1 is avg of row (previously), now doing avg of col
	val = stats.ttest_rel(withinMaxDF['max similarity'], btwnMaxDF['max similarity'], alternative='greater')
	result = (1-val.pvalue)
	return result

#function calls:
parser = argparse.ArgumentParser()
parser.add_argument('mg',metavar='moduleGenes', type=str,help='The file name and the path to the file of module genes. 1st argument')
parser.add_argument('umi',metavar='umiCountFile', type=str, help="The path to and file that is comprised of the 1st line with a path to the UMI count input files, and every subsequent line is the name of each cell type's UMI count file. 2nd argument")
parser.add_argument('dist', metavar='distanceMatrix',type=str,help='The path to and file containing the cell by cell distance matrix (can be generated in Seurat). 3rd argument')
parser.add_argument('o',metavar='outputFile',type=str,help='The output file name. 4th argument')
parser.add_argument('f',metavar='fileModifer',type=str,help='File modifier to save coexpression input in chunks; it is helpful to set this to a descriptor that indicates what run or set of module genes you are using. 5th argument')
parser.add_argument('h',metavar='comparison',type=str,help='Heatmap comparison argument, where 1 means output a csv with ranked UMI counts z-scored per cell type, make this a float value. 6th argument')
parser.add_argument('-c', metavar='coexpression', type=str,help='the path to the file plus the file name of the coexpression file, not currently in use, so optional')
args = parser.parse_args()
moduleGenesInputPath = args.mg
umiInputFile=args.umi
distMatrixInfile =args.dist
outfileName=args.o
heatmapComparison=args.h
fileModifier = args.f
comparison = float(heatmapComparison)
moduleGenes = getModuleGenes(moduleGenesInputPath)
coexpEnabled = False#set this to True if you want to take into account real coexpression information, then getCoexpressionDataChunks will read it in properly
coexpDF = getCoexpressionDataChunks(moduleGenes, fileModifier, coexpEnabled)
umiFilesList = getUMIFiles(umiInputFile)
outfile = open(outfileName, 'w')
outfile.write('CellType\tLog\tRaw\n')
for cellTypeFile in umiFilesList:
	line = cellTypeFile.split('/')
	cellType = line[-1]
	cellType = cellType[:-13]#remove the UMICounts.txt from the end (i.e. EndUMICounts.txt becomes cell type of End)
	print('cell type:', cellType)
	umiDF = getUMICount(cellTypeFile)
	barcodes = umiDF.columns
	btwnMatrix, withinMatrix = getDistanceMatrix(barcodes, distMatrixInfile)#note withinMatrix is the prime distance matrix, btwnMatrix is the matrix of comparisons amongst the same cell type
	rhsVal = rightHandSidePValue(btwnMatrix, withinMatrix)
	lhsDict, geneSum = leftHandSide(moduleGenes, umiDF, coexpDF,comparison)
	objectiveFuncValueLog, rawObjectiveFunctionValue = multiplication(lhsDict, rhsVal, comparison, fileModifier, cellType)#without coexpression
	if comparison != 1:
		outfile.write(cellType + '\t' + str(objectiveFuncValueLog) + '\t' + str(rawObjectiveFunctionValue) + '\n')
	else:
		objectiveFuncValueLog.to_csv(fileModifier +cellType+ '.txt', sep='\t')
outfile.close()
