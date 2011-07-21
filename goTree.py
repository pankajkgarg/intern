import math, csv, array, operator, json, os, random
import gene_ontology as go
import mds
from utils import *
from microarray import MicroArray
from gzip import GzipFile
from collections import defaultdict

#if __name__ != "__main__":
#	import sys
#	sys.stdout = open("data/debug.txt", "ab")

def annotations(fileObj, evidenceCodes = []):
	'''
	Takes only 3 things from the annotation file, db_object_id and go_id, while evidence_code are used for filterDataing
	Returns a dict mapping go_id to a set of db_object_id
	
	fileObj - file containing annotations, can be filename, or file object
	evidenceCodes - the evidence codes to select, IF the list is left empty, then all are allowed
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
	def __init__(self, termObj, genes = [], totalGenes = 0):
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
		# try:
			# self.ic = -1 * math.log(float(self.numGenes)/totalGenes)	
		# except Exception, info:
			# print info
			# raise Exception('Trouble with information content  ' + str(self.numGenes) + '  ' +  str(totalGenes))
		
		if self.numGenes:
			self.numCols = len(self.genes[0].profile)
		else:
			self.numCols = 0
		
		if "is_a" in self.tags:
			#Only compute the correlation for non root nodes
			# Since, they are anyway left out
			self.computeCorrelation()	
		else:
			self.meanCorrelation = 0.0
			self.correlationList = None
		
		
	def computeCorrelation(self):
		if self.numGenes == 0:
			return
		
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
			uid,
			name = 'results',
			evidenceCodes = []
		):
		'''
		Super class, It combines all other classes to reach the final point.
		
		microarray - An instance of MicroArray class
		name - a unique name for the resultant file
		'''
		
		self.significanceLevel = 0.05
		
		self.uid = uid
		self.name = name
		
		self.ontologyFile = ontologyFile
		self.annotationFile = annotationFile
		self.microarray = microarray
		
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
				#Probably no data found in the microarray for the 
				#	corresponding set of genes
				continue
			
			tempTerm = GoTerm(tempTerm, modGenes, self.totalGenes)
			self.terms[goId] = tempTerm
		
		print 'No. of terms before filtering:', len(self.terms)
		#Terms after filtering
		self.modTerms = modTerms = self.pairFilter()
		print 'No. of terms after filtering: ', len(modTerms)
		
		simMatrix = self.resnick(modTerms)
		
		#Converting into dissimilarity matrix
		disMatrix = self.sim2dis(simMatrix)
		
		#Doing mds
		mdsInstance = mds.MDS(disMatrix)
		coordinates = mdsInstance.process()
		
		for i,row in enumerate(coordinates): 
			#TODO: Look at the random noise added below			
			modTerms[i].x = row[0] + (random.random() * row[0] * 0.1)
			modTerms[i].y = row[1] + (random.random() * row[1] * 0.1)
		
		zoomLabels = self.subSelect(modTerms)
		
		self.save(modTerms, zoomLabels)
		
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
		
		
	def save(self, terms, zoomLabels):
		data = []
		for i,term in enumerate(terms):	
			data.append({
				"i": i,
				"x": term.x,
				"y": term.y,
				"goId": term.id,
				"label": term.name,
				"score": term.meanCorrelation,
			})
		
		f = open( os.path.join('results', self.uid + '_data.js') , "wb")
		text = "var data = " + json.dumps(data)
		f.write(text)
		f.close()
		
		
		f = open( os.path.join('results', self.uid + '_zoomLabels.js') , "wb")
		text = "var zoomLabels = " + json.dumps(zoomLabels)
		f.write(text)
		f.close()
		
		info = {}
		info['uid'] = self.uid
		info['name'] = self.name
		info['microArrayFile'] = self.microarray.fileObj.name
		if isinstance(self.ontologyFile, file):
			info['ontologyFile'] = self.ontologyFile.name
		else:
			info['ontologyFile'] = self.ontologyFile
		
		if isinstance(self.annotationFile, file):
			info['annotationFile'] = self.annotationFile.name
		else:
			info['annotationFile'] = self.annotationFile
		
		f = open(os.path.join('results', self.uid + '_info.js') , "wb")
		text = "var info = " + json.dumps(info)
		f.write(text)
		f.close()
		
	
	def pairFilter(self):
		modTerms = set()
		
		#Finding out the root nodes
		rootNodes = set()
		for goId, term in self.terms.iteritems():
			if "is_a" not in term.tags:
				rootNodes.add(term.id)
		
		#Finding out the direct descendants of root node
		underRoot = set()		
		for goId, term in self.terms.iteritems():
			if "is_a" not in term.tags:
				# Must be a root node, and hence no parents
				continue
			
			parents = term.tags['is_a']
				
			for parent in parents:
				parentId = parent.value
				if parentId in rootNodes:
					underRoot.add(term.id)
					
		for goId, term in self.terms.iteritems():
			if goId in rootNodes or goId in underRoot:
				#Direct descendants are being neglected now, 
				# but an exception is made for their children
				# i.e. If their children doesn't pass the t-test
				# then the direct descendant can be selected
				continue
			
			parents = term.tags['is_a']
				
			for parent in parents:
				parentId = parent.value
				if parentId in self.terms:	

					tValue = tTest(term.correlationList, self.terms[parentId].correlationList)
					
					if tValue < 0 and parentId not in underRoot:
						# The mean of the parent is greater, so neglecting the term
						continue
					
					modT = math.fabs(tValue)
					
					#Using one tailed pValue
					pValue = pCalculator(modT, len(term.correlationList) + len(self.terms[parentId].correlationList) - 2)
					
					print tValue, pValue
					
					if math.isnan(pValue):
						#TODO: WHY is this happening.. why am I getting nan values
						pValue = 0.001
					
					if pValue < self.significanceLevel and tValue > 0:
						modTerms.add(term)
					elif parentId in underRoot and pValue  < self.significanceLevel and tValue < 0:
						modTerms.add(self.terms[parentId])
					else:
						print 'Not selected: ',pValue
					
		modTerms = list(modTerms)		
		return modTerms
		
	def greedyFilter(self):
		termsList = []
		for term in self.terms.itervalues():
			termsList.append((term, term.meanCorrelation))
			
		termsList = sorted(termsList, key = operator.itemgetter(1), reverse = True)
		
		modTerms = []
		for term, meanCorrelation in termsList[:200]:
			modTerms.append(term)
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
		
	def subSelect(self, terms):
		'''
		Implements the algorithm for selecting summarizing labels.
		Takes the terms and return the labels
		'''
		if not isinstance(terms, list):
			raise Exception("Input to subSelect, must be a list")
		
		children = defaultdict(set)
		for i, term in enumerate(terms):
			tempAncestors = self.tree.ancestors(term.id)
			for ancestor in tempAncestors:
				 children[ancestor.id].add(term.id)
		
		ancestorScores = [(termId, len(temp)) for termId, temp in children.iteritems() ]
		
		#Finding the immediate children of all terms in ancestor
		directChildren = defaultdict(set)
		ancestorsId = children.keys()
		for termId in ancestorsId:
			term = self.tree.ensure_term(termId)
			if "is_a" in term.tags:
				parents = term.tags["is_a"]
				for parent in parents:
					parentId = parent.value
					directChildren[parentId].add(termId)
					
		
		#Find minimum and maximum score
		minScore = min(ancestorScores, key = operator.itemgetter(1))[1]
		maxScore = max(ancestorScores, key = operator.itemgetter(1))[1]
		
		labels = []
		lastLabelSet = set()
		#Dividing data into 10 bins
		binSize = (maxScore - minScore)/15.0
		i = lastI = maxScore
		while i >=  (minScore-binSize):
			selectedIds = [ancestorId for ancestorId, score in ancestorScores if score > i]
			
			selectedSet = set(selectedIds)
			
			leafNodesId = []
			for termId in selectedIds:
				isLeaf = True
				tempChildren = list(directChildren[termId])
				for childId in tempChildren:
					if childId in selectedSet:
						isLeaf = False
						break
				if isLeaf:
					leafNodesId.append(termId)
			
			if len(leafNodesId) == 0:
				lastI = i
				i -= binSize
				continue
			
			finalSelection = []
			alreadyIn = set()
			remaining = leafNodesId[:]
			target = 0.9 * len(terms)
			while True:
				if len(alreadyIn) > target:
					break
				maxLen = 0
				maxId = None
				maxSet = None
				
				for nodeId in remaining:
					uniqueNodes = children[nodeId].difference(alreadyIn)
					if len(uniqueNodes) > maxLen:
						maxLen = len(uniqueNodes)
						maxSet = uniqueNodes
						maxId = nodeId
				
				if maxSet is None:
					break
				
				finalSelection.append(maxId)
				alreadyIn.update(maxSet)
			
			
			if len(finalSelection) >  2* len(lastLabelSet) and len(lastLabelSet) > 3 and lastI != i and binSize/(lastI - i) < 10:
				i = i + ((lastI - i)/2.0)
				continue
			elif set(finalSelection) != lastLabelSet:
				if len(lastLabelSet) == len(finalSelection):
					labels.pop()
				labels.append(finalSelection)
				lastLabelSet = set(finalSelection)
				
				# print "Initially ", len(selectedIds), "Leaf nodes", len(leafNodesId), "Final Selection", len(finalSelection)
				# print 'leafNodes', leafNodesId
				# print 'Final selection',  finalSelection
				
			lastI = i
			i -= binSize
				
			
		
		termDict = {}
		for term in terms:
			termDict[term.id] = term
		
		# Key will represent zoom level, and values will be a set of terms 
		#	with properties x coordinate, y coordinate and name, and go id
		finalLabels = {}
		for zoomLevel, labelList in enumerate(labels):
			tempLabels = [] 
			for labelId in labelList:
				meanX = 0.0
				meanY = 0.0
				for childId in children[labelId]:
					meanX += termDict[childId].x
					meanY += termDict[childId].y
				
				tempTerm = {}	
				tempTerm["x"] = meanX/float(len(children[labelId]))
				tempTerm["y"] = meanY/float(len(children[labelId]))
				tempTerm["score"] = len(children[labelId])
				tempTerm["label"] = self.tree.ensure_term(labelId).name
				tempTerm["goId"] = labelId
				tempLabels.append(tempTerm)
			
			finalLabels[zoomLevel] = tempLabels
			
		return finalLabels
		
if __name__ == "__main__":
	import sys, time
	t = time.time()
	microarray = MicroArray('data/Birnbaum.csv')
	
	temp = GoTree('data/gene_ontology.obo', 'data/gene_association.tair', microarray, "43jo")
	
	print 'Total time taken: ', time.time() - t
		
