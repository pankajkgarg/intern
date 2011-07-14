import math, csv, array, operator
import goTree
import gene_ontology as go

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


	
class GoTerm():
	"This class takes care of t-test for a Go Tree node"
	
	def __init__(self, termObj, genes, totalGenes):
		'termObj - instance of go.Term'
		'genes - A list or set of instances of Gene'
		'totalGenes - total no. of genes in the annotation file (to calcuate the Information content)'
		self.id = termObj.id
		self.name = termObj.name
		self.tags = termObj.tags
		
		self.genes = genes
		
		self.numGenes = len(self.genes)
		
		try:
			self.ic = -1 * math.log(float(self.numGenes)/totalGenes)	# Information content
		except:
			print float(self.numGenes), totalGenes
			raise Exception
		
		if self.numGenes:
			self.numCols = len(self.genes[0].profile)
		else:
			self.numCols = 0
		
		self.computeCorrelation()	
		
		#Take care of the t-test here
		
	def computeCorrelation(self):
		
		if self.numGenes == 1:
			self.correlationList = [1.0]
			self.meanCorrelation = 1.0
			return
			#return [[1]]
		
		correlationList = array.array('f', [0] * int(((self.numGenes * self.numGenes) - self.numGenes)/2) )
		
		k = 0
		for i in xrange(self.numGenes):			
			#correlationMatrix[i][i] = 1.0
			
			if i%100 == 0:
				print i, '\tTotal: ', self.numGenes
			
			for j in xrange(i):	
				numerator = sum(map(operator.mul, self.genes[i].meanDiff, self.genes[j].meanDiff ))
				denominator = self.genes[i].varianceLike * self.genes[j].varianceLike 
				
				correlation = numerator/denominator
				modCorrelation  = math.fabs(correlation)
				
				#iMPORTANT: substitute append by something else, because appends are slow
				correlationList[k] = (modCorrelation)
				k += 1
		
		self.meanCorrelation = float(sum(correlationList))/len(correlationList)
		self.correlationList = correlationList
		
class GoTree:
	def __init__(self, 
			ontologyFile,  
			annotationFile, 
			microarray,
			evidenceCodes = []
		):
		'''
		Super class, It combines all other classes to reach the final point.
		
		microarray - An instance of MicroArray class
		'''
		
		self.significanceLevel = 0.05
		
		self.ontologyFile = ontologyFile
		self.annotationFile = annotationFile
		
		self.totalGenes = microarray.numGenes
		
		self.associations = goTree.annotations(self.annotationFile, evidenceCodes)
		self.tree = go.Tree.from_obo(ontologyFile)
		
		numProcessed = 0
		self.terms = {}
		for goId, geneSet in self.associations.iteritems():
			numProcessed += 1
			if numProcessed % 100 == 0:
				print 'Num processed: ', numProcessed
			
			tempTerm = self.tree.ensure_term(goId)
			
			modGenes = []
			genes = list(geneSet)
			for geneId in genes:
				tempGene = microarray.gene(geneId)
				if tempGene is not None:
					modGenes.append(tempGene)
			
			if len(modGenes) < 3:
				#Probably no data found in the microarray for the corresponding set of genes
				continue
			
			tempTerm = GoTerm(tempTerm, modGenes, self.totalGenes)
			self.terms[goId] = tempTerm
		
		print 'No. of terms before filtering:', len(self.terms)
		#Terms after filtering
		modTerms = self.filtering()
		print 'Number of terms after filtering: ', len(modTerms)
		
		simMatrix = self.resnick(modTerms)
		
		#Converting into dissimilarity matrix
		tempList = []
		[tempList.extend(i) for i in simMatrix]
		maxValue = float(max(tempList))
		disMatrix = [[0 for _ in xrange(len(modTerms))] for _ in xrange(len(modTerms))]
		for i in xrange(len(modTerms)):
			for j in xrange(i):
				disMatrix[i][j] = disMatrix[j][i] = 1 - (simMatrix[i][j]/maxValue)
		
		
		#Doing mds
		mds = MDS(disMatrix)
		mds.process()
		
		#self.terms = {}
		#for term in modTerms:
		#	self.terms[term.id] = term
		
	
	def filtering(self):
		modTerms = set()
		for goId, term in self.terms.iteritems():
			try:
				parents = term.tags['is_a']
			except:
				print term.tags
				print "\n\n"
				continue
				
			for parent in parents:
				parentId = parent.value
				if parentId in self.terms:
					tValue = self.tTest(term.correlationList, self.terms[parentId].correlationList)
					
					if tValue < 0:
						modT = -1 * tValue
					else:
						modT = tValue
					
					#Using one tailed pValue
					pValue = self.pCalculator(modT, len(term.correlationList) + len(self.terms[parentId].correlationList) - 2)
					
					print tValue, pValue
					if math.isnan(pValue):
						pValue = 0.001
					
					if pValue < self.significanceLevel and tValue > 0:
						modTerms.add(term)
					else:
						print pValue
					
					
		return modTerms
						
	def resnick(self, terms):
		'Terms - List of filtered go terms ( Instance of GoTerm )'
		lenTerms = len(terms)
		simMatrix = [[0 for _ in xrange(lenTerms)] for _ in xrange(lenTerms)]
		
		terms = list(terms)
		
		ancestors = []
		for i in xrange(lenTerms):
			ancestors.append(self.tree.ancestors(terms[i]))
			
		for i in xrange(lenTerms):
			for j in xrange(i+1):
				common = ancestors[i].intersection(ancestors[j])
				scores = []
				for term in common:
					scores.append( len(self.associations[term.id]) ) 
				
				score = min(scores) # minimum because, lesser the number of genes, higher the Information Content
				simMatri[i][j] = simMatrix[j][i] = -1 * math.log(float(score)/self.totalGenes)
				
		return simMatrix
		
					
		
	def tTest(self, series1, series2):
		num1 = len(series1)
		num2 = len(series2)
		
		mean1 = sum(series1)/float(num1)
		mean2 = sum(series2)/float(num2)
		numerator = mean1 - mean2
		
		variance1 = sum( pow((mean1 - i), 2) for i in series1) / float(num1)
		variance2 = sum( pow((mean2 - i), 2) for i in series2) / float(num2)
		denominator = pow( (pow(variance1, 2)/num1) + (pow(variance2, 2)/num2), 0.5)
		
		try:
			result = (numerator/denominator)
		except:
			print variance1
			print series1
			print variance2
			print series2
			
			
		return result
			
		
	memoize = {}
	def dfParam(self, df):
		'Gives the part dependent upon degree of freedom in the probability distribution function'
		df = int(df)
		if df not in self.memoize:
			if (df%2) == 0:
				#if even
				factor = 0.5
			else:
				# if odd
				factor = 1.0 / math.pi
				
			numerator = 1.0
			denominator = math.sqrt(df)
			for i in xrange(df - 1, 1, -2):
				numerator = numerator * i
			for i in xrange(df - 2, 1, -2):
				denominator = denominator * i
			
			result = (float(numerator)/denominator) * factor
			self.memoize[df] = result
		return self.memoize[df]
		
	def tDist(self, x, df):
		'Gives the probability density function of t-distribution'	
		return  self.dfParam(df) * pow( (1 + (x*x/float(df))), -1 * 0.5 * (df + 1) )
	
	def pCalculator(self, tValue,  df):
		'''
		The p value is calculated by measuring the area under the curve from the tValue to a certain high value.
		So, always give a positive tValue
		'''
		#Integration through Simpson rule
		initialX = tValue
		finalX = tValue + 10.0
		numIntervals = 500
		h = (finalX - initialX)/float(numIntervals)
		
		pValue = self.tDist(initialX, df) + self.tDist(finalX, df)
		tempSum = 0
		for i in xrange(1, numIntervals/2):
			tempSum += self.tDist(initialX + (h*2*i), df)
		pValue += (2*tempSum)
		tempSum = 0
		for i in xrange(1, (numIntervals/2) + 1):
			tempSum += self.tDist(initialX + (h*((2*i) - 1)), df)
		pValue += 4*tempSum
		
		pValue = (h/3.0) * pValue
		return pValue	#only gives the one sided p value
		
		#Integration by a method which assumes the function to be constant
		# in small intervals (Also gives good results but simpsons is established)
		
		# pValue = 0
		# numIntervals = 2000
		# initialX = tValue
		# finalX = tValue + 10.0
		
		# interval = (finalX - initialX)/ numIntervals
		# i = tValue
		# while i < (tValue + 10.0):
		
			# temp = self.tDist(i, df)
			# if temp < 0.000001:
				# break
			# pValue += (temp * interval)
			
			# i += interval
			
		# return pValue
		
if __name__ == "__main__":
	import sys, time
	
	microarray = MicroArray('data/Birnbaum.csv')
	
	temp = GoTree('data/gene_ontology.obo', 'data/gene_association.tair', microarray)
	print sys.argv
	t = time.time()
	print temp.pCalculator(float(sys.argv[1]), int(sys.argv[2]))
	print time.time() - t
		
			
		
		
	
