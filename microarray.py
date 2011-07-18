import math, csv, array, operator, json, os
import mds


class MicroArray:
	def __init__(self, fileObj, **kwargs):
		'fileObj - a filename or file object referencing to the microarray csv file'
		if isinstance(fileObj, str):
			self.fileObj = open(fileObj, "rb")
		elif isinstance(fileObj, file):
			self.fileObj = fileObj
		else:
			raise Exception('input to Microarray is not of type string or file')
		
		self.genes = self.fileReader(self.fileObj, **kwargs)
		self.numGenes = len(self.genes)
		
		uidCount = 1
		for geneId in self.genes.iterkeys():
			self.genes[geneId].uid = uidCount
			uidCount += 1
		
	def gene(self, geneId):
		if geneId in self.genes:
			return self.genes[geneId]
		else:
			return None
	
	def fileReader(self, fileObj, idCol = 1, ignoreRows = [1], ignoreCols = [2]):
		'''
			values of idCol, ignoreRows, and ignoreCols begin from 1 (instead of 0).... so 1 is the first column
		'''
		
		if isinstance(ignoreRows, int):
			ignoreRows = [ignoreRows]
		if isinstance(ignoreCols, int):
			ignoreCols = [ignoreCols]
		
		ignoreCols = set(ignoreCols)
		ignoreRows = set(ignoreRows)
		
		
		reader = csv.reader(fileObj)
		
		genes = {}
		for i,row in enumerate(reader):
			lineNum = i+1
			
			if lineNum in ignoreRows:
				continue
			
			if row[idCol - 1] == "":
				continue
			
			values = []
			tempId = ''
			for j, cell in enumerate(row):
				colNum = j + 1
				
				if colNum in ignoreCols:
					continue
				
				if colNum == idCol:
					tempId = cell
				else:
					if cell == '':
						cell = '0'
					values.append(int(cell))
			
			#Checking if the variance of the gene expression profile is zero, if it is so, then it is left out
			if self.variance(values) == 0:
				continue
			
			genes[tempId] = Gene(tempId, values)
						
		fileObj.close()
		return genes
		
	def variance(self, values):
		mean = sum(values)/float(len(values))
		return sum( pow((mean-i), 2) for i in values)/float(len(values))
	
	
class Gene:
	def __init__(self, id, profile):
		'profile refers to the row of values related to the gene in the microarray	'
		self.id = id
		self.profile = array.array('f', profile)
		if not isinstance(profile, list):
			print profile
			raise Exception(" gene profile is not a list")
		
		self.numCols = len(self.profile)
		self.mean = sum(self.profile)/float(len(self.profile))
		self.meanDiff = array.array('f', [self.profile[k] - self.mean for k in xrange(self.numCols)])
		self.varianceLike = pow( sum( pow((self.profile[k] - self.mean), 2) for k in xrange(self.numCols)), 0.5)


	
