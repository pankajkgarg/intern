'''
LICENSE:
	Copyright (C) 2011, Pankaj Kumar Garg
	This program is distributed under GNU General Public License
'''

__author__  = "Pankaj Kumar Garg"
__email__   = "pankajn17@gmail.com"
__copyright__ = "Copyright (c) 2011, Pankaj Kumar Garg"
__license__ = "GPLv3"

#First implementation of mds
#Non-metric mds using SMACOF algorithm

import os
os.chdir(os.path.dirname(__file__))

import random, time, os, operator, array
from pprint import pprint



class MDS:
	def __init__(self, matrix, numDims = 2):
		'matrix - Dissimalarity matrix'
		'numDims - number of Dimensions to create the mapping into'
		self.disMatrix = matrix		#disMatrix - dissimalarity matrix
		
		#Checking dimensions of matrix
		if len(matrix) != len(matrix[0]):
			raise Exception('Not a square matrix')
		self.numPoints = len(matrix)
		
		self.numDims = numDims
		
	def initialization(self):
		'''
		Initialize the distance matrix. The matrix is random, 
			except that the closest pairs are kept close.
		Returns a matrix of dimension 
			-> N(number of elements) * M(number of dimensions)
		'''
		#The maximum length of the graph is decided as the 
		#	square root of the number of points in the distance matrix
		maxCoord = self.numPoints 	
		
		
		#Finding closest pairs
		self.closestPairs = {}
		for i in xrange(self.numPoints):
			temp = self.disMatrix[i][:]
			temp[i] = 10000000
			self.closestPairs[i] = min(enumerate(temp), key = operator.itemgetter(1))[0]
				
		
		matrix = []
		for i in xrange(self.numPoints):
			matrix.append(array.array('f'))
			for j in xrange(self.numDims):
				if self.closestPairs[i] < i :
					matrix[i].append(matrix[self.closestPairs[i]][j] + (random.random() * maxCoord * 0.1 * random.choice([-1, 1])) )
				else:
					matrix[i].append( (random.random() * maxCoord ) )
		
		return matrix			
		#return [ [(random.random() *  maxCoord) for i in xrange(self.numDims)] for _ in xrange(self.numPoints)]
		
	def updateDistanceMatrix(self):
		'Shall update the distanceMatrix based on X(vector matrix)'
		matrix = [[0 for _ in xrange(self.numPoints)] for _ in xrange(self.numPoints)]
		for i in xrange(self.numPoints):
			matrix[i][i] = 0
			for j in xrange(i):
# Only for 2 dimension
				matrix[i][j] = matrix[j][i] = pow(pow((self.X[i][0] - self.X[j][0]), 2) + pow((self.X[i][1] - self.X[j][1]), 2) , 0.5)
				
				#matrix[i][j] = matrix[j][i] = sum((self.X[i][k] - self.X[j][k]) ** 2 for k in xrange(self.numDims)) ** 0.5
				
		self.distanceMatrix = matrix
	
	def updateB(self):
		'Updates the B matrix'
		matrix = [array.array('f', [0 for _ in xrange(self.numPoints)]) for _ in xrange(self.numPoints)]
		for i in xrange(self.numPoints):
			matrix[i][i] = 0
			for j in xrange(i):
				matrix[i][j] = matrix[j][i] = -1 * self.disMatrix[i][j] * (1.0/self.distanceMatrix[i][j]) if self.distanceMatrix[i][j] != 0 else 0
		for i in xrange(self.numPoints):
			matrix[i][i] = -1 * sum(matrix[i])
		self.B = matrix
	
	def multiplyMatrix(self, A, B):
		#Checking the dimensions of the matrix
		n1, m1 = len(A), len(A[0])
		n2, m2 = len(B), len(B[0])
		if m1 != n2:
			raise Exception ('matrices not compatible for multiplication')
		
		matrix = []
		for i in xrange(n1):
			matrix.append(array.array('f', ([0]*m2)))
			for j in xrange(m2):
				matrix[i][j] = sum(A[i][k] * B[k][j] for k in xrange(m1))
		return matrix
		
	def stress(self):
		#Returns the stress 
		total = 0.0
		for i in xrange(self.numPoints):
			for j in xrange(i):
				total += (self.distanceMatrix[i][j] - self.disMatrix[i][j]) ** 2
		return total
		
	def iteration(self):
		self.updateDistanceMatrix()
		self.updateB()
		
		nInverse = 1.0 / self.numPoints
		#print 'before multiplication'; pprint(self.distanceMatrix); print; pprint(self.B); print ; pprint(self.X)
		self.X = self.multiplyMatrix(self.B, self.X)
		#print 'after multiplication';		pprint (self.X)
		self.X = [ [(nInverse * self.X[i][j]) for j in xrange(self.numDims)] for i in xrange(self.numPoints)]
		
	def process(self, plot = False):
		#Retruns a n*2 matrix
		
		self.X = self.initialization()

		#Plotting the graph
#Assuming the data is 2-Dimensional
		if plot:
			import matplotlib.pyplot as plt
			plt.ion()
		
		self.iteration()
		xData = [i[0] for i in self.X]
		yData = [i[1] for i in self.X]
		
		if plot:
			p, = plt.plot(xData, yData, 'o')
		
		stress, oldStress = 0, 0
		for j in xrange(400):	#Putting maximum 400 iterations
			self.iteration()

			if plot and j%20 == 0 and j!= 0:
				xData = [i[0] for i in self.X]
				yData = [i[1]  for i in self.X]
				p.set_xdata(xData)
				p.set_ydata(yData)
		
				plt.draw()
			
			if j == 0:
				stress = self.stress()
			elif j == 1:
				firstDiff = oldStress - stress
			if j%4 == 0 and j!= 0:
				oldStress = stress
				stress = self.stress()
				if (oldStress - stress) < (0.005 * firstDiff):
					break
			
			if j%20 == 0:
				print 'Iteration: ', j, '\tStress: ', self.stress()
			#time.sleep(1)
		
		# Publishing coordinates to a file
		coordString = ""
		for i in self.X:
			coordString += "%f,%f\n" % (i[0], i[1])
		coordString = coordString[:-1]
		
		f = open("geneCoord.txt", "wb")
		f.write(coordString)
		f.close()
			
		if plot:
			plt.ioff()
			plt.show()
		
		return self.X

	
def test():
	numPoints = 100
	disMatrix = [[0 for _ in xrange(numPoints)] for _ in xrange(numPoints)]
	
	
	for i in xrange(numPoints):	
		disMatrix[i][i] = 0
		for j in xrange(i):
			disMatrix[i][j] = disMatrix[j][i] = int(random.random() * numPoints)
	
	random.seed(random.random() * 10000)
	#disMatrix = [[0,2,3], [2,0,1], [3,1,0]]
	#disMatrix = [[0,5], [5,0]]
	
	disMatrix = irisData()
	#pprint(disMatrix)
	
	mds	= MDS(disMatrix)
	mds.process(plot=True)

def irisData():
	'Shall return the iris data as dissimalarity matrix'
	f = open('mds_sampledata.txt', 'rb')
	data = f.read()
	f.close()
	data = data.split('\n')[3:-2]
	modData = [i.split('\t')[:4] for i in data]
	modData = [[float(j) for j in i] for i in modData]
	numDims = 4
	numPoints = len(modData)
	
	matrix = [ [0 for _ in xrange(numPoints)] for _ in xrange(numPoints)]	#Initializing matrix
	
	
	for i in xrange(numPoints):	
		for j in xrange(i):
			matrix[i][j] = matrix[j][i] = sum( (modData[i][k] - modData[j][k]) ** 2 for k in xrange(numDims) ) ** 0.5
			
	return matrix

def geneData():
	f = open("tempNodeSimMat.txt", "rb")
	rows = f.read().split("\n")
	f.close()
	
	allData = []
	simMatrix = []
	for row in rows:
		if len(row.strip()) == 0:
			continue
		simMatrix.append([float(node.strip()) for node in row.split(",")])
		
	for row in simMatrix:
		allData.extend(row)
	maxValue = max(allData)
	
	disMatrix = []
	matrixLen = len(simMatrix)
	for i in xrange(matrixLen):
		disMatrix.append([])
		for j in xrange(matrixLen):
			simMatrix[i][j] = 1.0/simMatrix[i][j]
		simMatrix[i][i] = 0
		disMatrix[i] = simMatrix[i]
	return disMatrix
	
if __name__ == '__main__':
	test()		

