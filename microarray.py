import math
import goTree
import gene_ontology as go

class MicroArray:
	def __init__(self, matrix, genesId, conditions):
		self.matrix = matrix
		self.conditions = conditions
		
		self.numGenes = len(genesId)
		self.numCols = len(self.matrix[0])
		
		self.genes = {}
		for i in xrange(self.numGenes):
			self.genes[genesId[i]] = Gene(genesId[i], matrix[i])
		
	def geneProfile(self, geneId):
		if geneId in self.genes:
			return self.genes[geneId]
		else:
			return None
	
	
	
		
class Gene:
	def __init__(self, id, profile):
		'profile refers to the row of values related to the gene in the microarray	'
		self.id = id
		self.profile = profile
	
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
		
		self.ic = -1 * math.log(float(self.numGenes)/totalGenes)	# Information content
		
		if self.numGenes:
			self.numCols = len(self.genes[0].profile)
		else:
			self.numCols = 0
			
		self.matrix = [tempGene.profile for tempGene in genes]
		
		self.correlationMatrix = self.correlation()	
		
		#Take care of the t-test here
		
	def correlation(self):
		correlationMatrix = [[0 for _ in xrange(self.numGenes)] for _ in xrange(self.numGenes)]
		
		meanArray = []
		for i in xrange(self.numGenes):
			meanArray.append(sum(self.matrix[i])/ float(len(self.matrix[i])))
			
		
		correlationList = []
		
		for i in xrange(self.numGenes):			
			correlationMatrix[i][i] = 1.0
			for j in xrange(i):	
				numerator = sum( (self.matrix[i][k] - meanArray[i]) * (self.matrix[j][k] - meanArray[j]) for k in xrange(self.numCols))
				denominator = (sum( pow((self.matrix[i][k] - meanArray[i]), 2) for k in xrange(self.numCols)) * sum( pow((self.matrix[j][k] - meanArray[j]), 2) for k in xrange(self.numCols)) ) ** 0.5
				correlationMatrix[i][j] = correlationMatrix[j][i] = numerator / denominator
				
				correlationList.append(numerator/denominator)
			
		self.correlationList = correlationList
		self.meanCorrelation = float(sum(correlationList))/len(correlationList)
		return correlationMatrix
	
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
		
		self.significanceLevel = 0.01
		
		self.ontologyFile = ontologyFile
		self.annotationFile = annotationFile
		
		self.totalGenes = microarray.numGenes
		
		self.associations = goTree.annotations(self.annotationFile, evidenceCodes)
		self.tree = go.Tree.from_obo(ontologyFile)
		
		self.terms = {}
		for goId, geneSet in self.associations.iteritems():
			tempTerm = self.tree.ensure_term(goId)
			
			modGenes = []
			genes = list(geneSet)
			for gene in genes:
				geneProfile = microArray.geneProfile(gene)
				if geneProfile is not None:
					modGenes.append(Gene(gene, geneProfile))
			
			self.terms[goId] = GoTerm(tempTerm, modGenes, self.totalGenes)
		
		#Terms after filtering
		modTerms = self.filtering()
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
			parents = term.tags['is_a']
			for parentId in parents:
				if parentId in self.terms:
					tValue = self.tTest(term.correlationList, self.terms[parentId].correlationList)
					
					if tValue < 0:
						modT = -1 * tValue
					else:
						modT = tValue
					
					#Using one tailed pValue
					pValue = self.pCalculator(modT, len(term.correlationList) + len(self.terms[parentId].correlationList) - 2)
					
					if pValue < significanceLevel and tValue > 0:
						modTerms.add(term)
					
		return modTerms
						
	def resnick(self, terms):
		'Terms - List of filtered go terms ( Instance of GoTerm )'
		lenTerms = len(terms)
		simMatrix = [[0 for _ in xrange(lenTerms)] for _ in xrange(lenTerms)]
		
		ancestors = []
		for i in xrange(lenTerms):
			ancestors.append(self.tree.ancestors(terms[i])
			
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
		
		return (numerator/denominator)
			
		
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
		finalX = tValue + 20.0
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
			
		# return (pValue * 2)
		
if __name__ == "__main__":
	import sys, time
	
	
	
	GoTree(
	
	temp = GoTree()
	print sys.argv
	t = time.time()
	print temp.pCalculator(float(sys.argv[1]), int(sys.argv[2]))
	print time.time() - t
		
			
		
		
	
