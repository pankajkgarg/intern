'Utility functions'

'''
LICENSE:
	Copyright (C) 2011, Pankaj Kumar Garg
	This program is distributed under GNU General Public License
'''

__author__  = "Pankaj Kumar Garg"
__email__   = "pankajn17@gmail.com"
__copyright__ = "Copyright (C) 2011, Pankaj Kumar Garg"
__license__ = "GPLv3"

import math
from pprint import pprint

def tTest(series1, series2):
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
def dfParam(df):
	'Gives the part dependent upon degree of freedom in the probability distribution function'
	global memoize
	
	df = int(df)
	
	#If df > 295, then it results in "nan" values, so putting a upper limit
	if df > 295:
		df = 295
	
	if df not in memoize:
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
		memoize[df] = result
	return memoize[df]
		
def tDist(x, df):
	'Gives the probability density function of t-distribution'	
	return  dfParam(df) * pow( (1 + (x*x/float(df))), -1 * 0.5 * (df + 1) )

dfReason = 0
totalNan = 0
def pCalculator(tValue,  df):
	'''
	The p value is calculated by measuring the area under the curve from the tValue to a certain high value.
	So, always give a positive tValue
	'''
	global dfReason, totalNan
	#Integration through Simpson rule
	initialX = tValue
	finalX = tValue + 10.0
	numIntervals = 500
	h = (finalX - initialX)/float(numIntervals)
	
	pValue = tDist(initialX, df) + tDist(finalX, df)
	tempSum = 0
	for i in xrange(1, numIntervals/2):
		tempSum += tDist(initialX + (h*2*i), df)
	pValue += (2*tempSum)
	tempSum = 0
	for i in xrange(1, (numIntervals/2) + 1):
		tempSum += tDist(initialX + (h*((2*i) - 1)), df)
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
	
		# temp = tDist(i, df)
		# if temp < 0.000001:
			# break
		# pValue += (temp * interval)
		
		# i += interval
		
	# return pValue
	
	
def benchmark(func):
	"""
	A decorator that print the time of function take
	to execute.
	"""
	import time
	def wrapper(*args, **kwargs):
		t = time.time()
		res = func(*args, **kwargs)
		print 'Executed the function ',func.__name__, 'in ', time.time()-t, 'seconds'
		return res
	return wrapper
