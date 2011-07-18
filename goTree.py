import math, csv, array, operator, json, os
import gene_ontology as go
import mds
from microarray import MicroArray
from gzip import GzipFile
from collections import defaultdict

def annotations(fileObj, evidenceCodes = []):
	'''
	Takes only 3 things from the annotation file, db_object_id and go_id, while evidence_code are used for filterDataing
	Returns a dict mapping go_id to a set of db_object_id
	
	fileObj - file containing annotations, can be filename, or file object
	evidenceCodes - the evidence codes to select, if the list is left empty, then all are allowed
	'''
	
	ann = go.AnnotationFile(fileObj)
	
	
	if len(evidenceCodes) == 0:
		filterData = False
	else:
		filterData = True
		evidenceCodes = set(evidenceCodes)
		
	associations = defaultdict(set)	# Maps go_id to a set of genes id
	for temp in ann:
		#filterDataing based on evidence codes
		if filterData  and temp.evidence_code not in evidenceCodes:
			continue
		
		associations[temp.go_id].add(temp.db_object_symbol)		#IMPORTANT: change to db_object_id
	
	return associations

		

	
class GoTerm():
	"This class takes care of t-test for a Go Tree node"
	
	def __init__(self, termObj, genes, totalGenes):
		''''
		termObj - instance of go.Term
		genes - A list or set of instances of Gene
		totalGenes - total no. of genes in the annotation file 
			(to calcuate the Information content)
		'''
		self.id = termObj.id
		self.name = termObj.name
		self.tags = termObj.tags
		
		self.genes = genes
		
		self.numGenes = len(self.genes)
		
		# Information content
		try:
			self.ic = -1 * math.log(float(self.numGenes)/totalGenes)	
		except Exception, info:
			print info
			raise Exception('Trouble with information content  ' + str(self.numGenes) + '  ' +  str(totalGenes))
		
		if self.numGenes:
			self.numCols = len(self.genes[0].profile)
		else:
			self.numCols = 0
		
		self.computeCorrelation()	
		
		
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
			
			#if i%100 == 0:
			#	print i, '\tTotal: ', self.numGenes
			
			for j in xrange(i):	
				numerator = sum(map(operator.mul, self.genes[i].meanDiff, self.genes[j].meanDiff ))
				denominator = self.genes[i].varianceLike * self.genes[j].varianceLike 
				
				correlation = numerator/denominator
				modCorrelation  = math.fabs(correlation)
				
				correlationList[k] = modCorrelation
				k += 1
		
		self.meanCorrelation = float(sum(correlationList))/len(correlationList)
		self.correlationList = correlationList
		

		
class GoTree:
	def __init__(self, 
			ontologyFile,  
			annotationFile, 
			microarray,
			name = 'results',
			evidenceCodes = []
		):
		'''
		Super class, It combines all other classes to reach the final point.
		
		microarray - An instance of MicroArray class
		name - a unique name for the resultant file
		'''
		
		self.significanceLevel = 0.05
		
		self.name = name
		
		self.ontologyFile = ontologyFile
		self.annotationFile = annotationFile
		
		self.totalGenes = microarray.numGenes
		
		self.associations = annotations(self.annotationFile, evidenceCodes)
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
		disMatrix = self.sim2dis(simMatrix)
		
		#Doing mds
		mdsInstance = mds.MDS(disMatrix)
		coordinates = mdsInstance.process()
		
		self.save(modTerms, coordinates)
		
		#self.terms = {}
		#for term in modTerms:
		#	self.terms[term.id] = term
	
	def sim2dis(self, matrix):
		"Given a similarity matrix, converts it to dissimilarity matrix"
		size = len(matrix)
		disMatrix = [[0 for _ in xrange(size)] for _ in xrange(size)]
		
		tempList = []
		map(tempList.extend, matrix)
		maxValue = float(max(tempList))
		
		for i in xrange(size):
			for j in xrange(i):
				disMatrix[i][j] = disMatrix[j][i] = 1 - (matrix[i][j]/maxValue)
		return disMatrix
		
		
	def save(self, terms, coordinates):
		data = []
		for i,term in enumerate(terms):	
			data.append({
				"i": i,
				"x": coordinates[i][0],
				"y": coordinates[i][1],
				"label": term.name,
				"score": term.meanCorrelation,
			})
		
		f = open( os.path.join('results', self.name + '.js') , "wb")
		text = "var data = " + json.dumps(data)
		f.write(text)
		f.close()
		
	
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
						#TODO: WHY is this happening.. why am I getting nan values
						pValue = 0.001
					
					if pValue < self.significanceLevel and tValue > 0:
						modTerms.add(term)
					else:
						print 'Not selected: ',pValue
					
		modTerms = list(modTerms)		
		return modTerms
						
	def resnick(self, terms):
		'Terms - List of filtered go terms ( Instance of GoTerm )'
		
		#Ensure that the terms is a list, and not a set
		#	This is necessary, so that later labels doesn't mix up
		if not isinstance(terms, (list, tuple)):
			raise Exception("given parameter is not a list")
		
		lenTerms = len(terms)
		simMatrix = [[0 for _ in xrange(lenTerms)] for _ in xrange(lenTerms)]
		
		ancestors = []
		for i in xrange(lenTerms):
			ancestors.append(self.tree.ancestors(terms[i].id))
			
		for i in xrange(lenTerms):
			for j in xrange(i+1):
				common = ancestors[i].intersection(ancestors[j])
				scores = []
				for term in common:
					scores.append( len(self.associations[term.id]) ) 
				
				if len(scores) != 0:
					score = min(scores) # minimum because, lesser the number of genes, higher the Information Content
				else:
					score = 0
					
				if score == 0:
					#TODO: Analyze this situation, how come the minimum score is 0
					simMatrix[i][j] = simMatrix[j][i] = 0
				else:
					simMatrix[i][j] = simMatrix[j][i] = -1 * math.log(float(score)/self.totalGenes)
				
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
	
		
