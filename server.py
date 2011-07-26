# -*- coding: utf-8 -*-
'''
LICENSE:
	Copyright (C) 2011, Pankaj Kumar Garg
	This program is distributed under GNU General Public License
'''

__author__  = "Pankaj Kumar Garg"
__email__   = "pankajn17@gmail.com"
__copyright__ = "Copyright (c) 2011, Pankaj Kumar Garg"
__license__ = "GPLv3"

import os,sys
os.chdir(os.path.abspath(os.path.dirname(__file__)))
curdir = os.path.dirname(__file__)
if curdir not in sys.path:
    sys.path.append(curdir)

import web
from web import form
import goTree
from jinja2 import Environment, PackageLoader
from microarray import MicroArray
from fetch import fetch
import re, Queue, threading, string, random, os, sqlite3, json, time, base64
from pprint import pprint

#render = web.template.render('templates/')
env = Environment(loader=PackageLoader('__main__', 'templates'))

allowed = (
    ('admin','biology'),
    ('a','a'),
)


urls = (
	'/', 'home',
	'/login', 'Login',
	'/view/(?P<jobId>[^/]+)', 'downloadView',
	#'/downloadView/(?P<jobId>[^/]+)', 'downloadView',
	'/results/(?P<fileName>[^/]+)', 'results',
	'/history', 'history',
	'/deleteAll', "clearHistory",
	'/deleteRow/(?P<id>[^/]+)', 'deleteRow',
)

web.config.debug = False
app = web.application(urls, globals())

curdir = os.path.dirname(__file__)
session = web.session.Session(app, web.session.DiskStore(os.path.join(curdir, 'data', 'sessions')),)

def dbConn():
	conn = sqlite3.connect(os.path.join('results', 'results.db'))
	conn.row_factory = sqlite3.Row
	return conn

def db(func):
	'A decorator to manage sqlite connections'
	'Only works upon methods'
	
	def wrapper(*args, **kwargs):
		conn = dbConn()
		cursor = conn.cursor()
		kwargs["cursor"] = cursor
		result = func(*args, **kwargs)
		conn.commit()
		cursor.close()
		return result
	return wrapper
		

def dbSetup():
	'Function to setup required tables etc.'
	conn = dbConn()
	cursor = conn.cursor()
	cursor.execute('create table if not exists results (id text, name text, description text, timestamp int)')
	conn.commit()
	cursor.close()

dbSetup()

file2Url = {}
def annotationChoices():
	'If convert is true, data shall be returned in pickled format'
	global file2Url
	url = 'http://www.geneontology.org/GO.downloads.annotations.shtml'
	filePath = getFile('annotationList.html', url, cacheTime = 8*3600)
	f = open(filePath)
	text = f.read()
	f.close()
	
	table = re.search(r'(?is)<table[^<]+?id\s*=\s*"annot".+?</table>', text)
	if table is not None:
		table = table.group()
		annotations =  re.findall('(?is)<tr>.+?</tr>', table)
		for annotation in annotations:
			name = re.search(r'(?i)<span\s+class\s*=\s*[\'\"]spp[\'\"]\s*>(.+?)</span>', annotation)
			if name is None:
				continue
			name = name.group(1)
			name = re.sub(r'\s+', '_', name)
			
			link = re.search(r'(?i)<a\s+href\s*=\s*"(.+?)"\s*>\s*annotations\s*</a>', annotation)
			if link is None:
				continue
			link = link.group(1)
			
			file2Url[name] = link
			
	return file2Url	
	
	
def getFile(fileName, url, cacheTime = 7*24*3600):
	'''Returns the path of filename on the disk, if the file 
		doesn't exist, then it fetches them, and if the file is too
		old, then it replace them with new one.
	'''
	filePath = os.path.join('data', fileName)
	if os.path.isfile(filePath):
		info = os.stat(filePath)
		if (time.time() - info.st_mtime) < cacheTime:
			return filePath
	
	print 'Fetching url: ', url
	urlData = fetch(url)
	if urlData['responseCode'] != 200:
		raise Exception("Couldn't fetch url: ", url)
	
	f = open(filePath, 'wb')
	f.write(urlData['html'])
	f.close()
	
	return filePath
	

vemail = form.regexp(r"^$|.*@.*", "must be a valid email address")
vnums = form.regexp(r"[0-9]*\s*,?\s*", 'must be a list of comma separatad numbers e.g "2,7,12" ')
vchars = form.regexp(r"[^,]*\s*,?\s*", 'must be a list of comma separatad evidence codes e.g "exp,abc" ')
vnum = form.regexp(r"[0-9]+", "must be a column number")

vfile = form.regexp(r".+", "Required!")


inputForm = form.Form(
	#form.File("annotationFile", vfile, description = 'Annotation File'),
	#form.File('ontologyFile', vfile, description = 'Ontology File'),
	form.Dropdown("annotation", args = annotationChoices().keys(), description="Annotation", class_ ="chzn-select"),
	
	form.File('microarrayFile', vfile, description = 'Microarray file (Only csv files)'),
	form.Textbox("ignoreRows", vnums, description = "Rows in microarray to ignore"),
	form.Textbox("ignoreCols", vnums, description = "Cols in microarray to ignore"),
	form.Textbox("idCol", vnum, description = "Column containing the Gene Id"),
	
	# form.Checkbox("treeTypes", value="biological_process", description="biological process", checked=True),
	# form.Checkbox("treeTypes", description="molecular function", value="molecular_function"),
	# form.Checkbox("treeTypes", description="cellular component", value="cellular_component"),
	form.Dropdown("treeTypes", args=["biological_process", "molecular_function", "cellular_component"], value="biological_process", description="Tree type", class_ ="chzn-select"),
	
	#form.Dropdown("evidenceCodes", args=["all", "exp"], value="all", description="Evidence Codes"),
	form.Textbox("evidenceCodes", vchars, description= "Evidence codes to ignore"),
	
	form.Textbox("jobName", description = "Job Name (Optional)"),
	form.Textarea("description", description="Description"),
	form.Textbox('email', vemail, description = 'Email (optional)'),
	form.Button("submit", type="submit", description="Submit"),
)

def queueBased(jobQueue):
	while True:
		kwargs = jobQueue.get()
		
		if kwargs == None:
			print 'stopping the process'
			break	
		
		return process(kwargs)
		
def process(kwargs):

	if 'ignoreRows' in kwargs and kwargs['ignoreRows'] != '':
		kwargs['ignoreRows'] = re.sub(r'\s+', '', kwargs['ignoreRows']).split(',')
		kwargs['ignoreRows'] = [int(temp) for temp in kwargs['ignoreRows']]
	else:
		kwargs['ignoreRows'] = []
	
	if 'ignoreCols' in kwargs and kwargs['ignoreCols'] != '':
		kwargs['ignoreCols'] = re.sub(r'\s+', '', kwargs['ignoreCols']).split(',')
		kwargs['ignoreCols'] = [int(temp) for temp in kwargs['ignoreCols']]
	else:
		kwargs['ignoreCols'] = []
	

	if 'idCol' in kwargs:
		kwargs['idCol'] = int(kwargs['idCol'])
		
	if 'evidenceCodes' in kwargs and kwargs['evidenceCodes'] != '':
		kwargs['evidenceCodes'] = re.sub(r'\s+', '', kwargs['evidenceCodes']).split(',')
	else:
		kwargs['evidenceCodes'] = []
	
	microarray = MicroArray(kwargs['microarrayFile'], idCol = kwargs['idCol'], ignoreRows = kwargs['ignoreRows'], ignoreCols = kwargs['ignoreCols'] )

	#Shall do the processing and shall save the results in 'results' directory 
	#	with the filename  "kwargs['name']" + ".js"

	temp = goTree.GoTree(kwargs['ontologyFile'], kwargs['annotationFile'], microarray, treeTypes = [kwargs["treeTypes"]], uid = kwargs['jobId'], name = kwargs['jobName'], evidenceCodes = kwargs['evidenceCodes'])
	
	
	#os.remove(kwargs['ontologyFile'])
	#os.remove(kwargs['annotationFile'])
	os.remove(kwargs['microarrayFile'])
	
	#Add to the database
	conn = dbConn()
	cursor = conn.cursor()
	cursor.execute('insert into results values (?, ?, ?, ?)', (kwargs['jobId'], kwargs['jobName'], kwargs['description'], int(time.time())) )
	conn.commit()
	cursor.close()
	
	
	if kwargs["email"] is not None and len(kwargs['email']) > 0:
		email(**kwargs)


def giveHistory(cursor, limit = 0):
	if limit != 0:
		limitString = " limit " + str(limit)
	else:
		limitString = ""
	results = cursor.execute('select * from results order by timestamp desc '+ limitString)
	results = cursor.fetchall()
	
	modResult = []
	for row in results:
		modResult.append(dict(row))
	
	#pprint(modResult)
	
	return modResult
		
			
class home:
	@db
	def GET(self, cursor):
		#if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
		#	raise web.seeother('/login')
		
		annotationChoices()
		f = inputForm()
		
		formText = f.render()
		return env.get_template("home.html").render(form = formText, history = giveHistory(cursor, 5))
		#return render.register(f)
		
	@db	
	def POST(self, cursor):
		#print type(web.data())
		#pprint(web.data())
		inputData = web.input(annotationFile = {}, ontologyFile = {})
		
		if inputData.annotationFile != {} and inputData.annotationFile.filename[-3:] == '.gz':
			annotationAdd = '.gz'
		else:
			annotationAdd = ''
		
		if inputData.ontologyFile != {} and inputData.ontologyFile.filename[-3:] == '.gz':
			ontologyAdd = '.gz'
		else:
			ontologyAdd = ''
		
		
		f = inputForm()
		if not f.validates():
			formText = f.render()
			return env.get_template("home.html").render(form = formText, history=giveHistory(cursor, 5))
		else:
			dataDict = dict(f.value)
			dataDict["jobId"] = randomId()
			
			
			#Post data contains file content
			# The following code writes that data into a local file and replaces 
			# the dictionary value by filename
			
			# ontologyPath = os.path.join("data", dataDict["jobId"] + "_ontology.obo")
			# f = open(ontologyPath, "wb")
			# f.write(dataDict["ontologyFile"])
			# f.close()
			
			#dataDict["ontologyFile"] = os.path.join("data", "gene_ontology.obo")
			
#			annotationPath = os.path.join("data", dataDict["jobId"] + "_annotations.txt" + annotationAdd)
#			f = open(annotationPath, "wb")
#			f.write(dataDict["annotationFile"])
#			f.close()
#			dataDict["annotationFile"] = annotationPath
			
			dataDict["ontologyFile"] = getFile("gene_ontology.obo", "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo")
			
			if dataDict["annotation"] not in file2Url:
				raise Exception("Annotation file name not in file2Url")
			dataDict["annotationFile"] = getFile(dataDict["annotation"] + '.gz', file2Url[dataDict["annotation"]])
			
			microarrayPath = os.path.join("data", dataDict["jobId"] + "_microarray.csv")
			f = open(microarrayPath, "wb")
			f.write(dataDict["microarrayFile"])
			f.close()
			dataDict["microarrayFile"] = microarrayPath
			
			print dataDict
			 
			#jobQueue.put(dataDict)
			process(dataDict)
			raise web.seeother('/view/' + dataDict["jobId"])


class view:
	@db
	def GET(self, jobId, cursor, withHistory = True):
		#if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
		#	raise web.seeother('/login')
		cursor.execute("select * from results where id = ?", (jobId, ))
		info = dict(cursor.fetchone())
		
		info['dataFile'] = '/results/' +  info['id'] + '_data.js'
		info['zoomLabelsFile'] = '/results/' + info['id'] + '_zoomLabels.js'
		
		if withHistory:
			history = giveHistory(cursor, 5)
		else:
			history = []
		
		return env.get_template("view.html").render(info = info, history = history)

class downloadView:
	@db
	def GET(self, jobId, cursor):
		temp = view()
		
		text = temp.GET(jobId, withHistory = False)
		
		def substitute(text):
			cssLinks = re.findall(r"(?i)<link.+?rel\s*=\s*[\"\']stylesheet[\"\'].+?href\s*=\s*[\"\'].+?[\"\'].*?>", text)
			
			for link in cssLinks:
				tempText = "<style type='text/css'>"
				
				cssFile = re.search(r"href\s*=\s*[\"\']/?(.+?)[\"\'].*?>", link)
				if cssFile is not None:
					cssFile = cssFile.group(1)
					#print cssFile
					f = open(cssFile, "rb")
					cssText = f.read()
					f.close()
				else:
					cssText = ""
				
				tempText += cssText + "</style>"
				tempText = tempText.decode("utf-8")
				#print link
				text = text.replace(link, tempText)
			
			jsLinks = re.findall(r"(?i)<script.+?src\s*=\s*[\"\'].+?[\"\'].*?>\s*</script>", text)
			
			for link in jsLinks:
				tempText = "<script type='text/javascript'>"
				
				jsFile = re.search(r"src\s*=\s*[\"\']/?(.+?)[\"\']", link)
				if jsFile is not None:
					jsFile = jsFile.group(1)
					#print jsFile
					f = open(jsFile, "rb")
					jsText = f.read()
					f.close()
				else:
					jsText = ""
				
				tempText += jsText + "</script>"
				tempText = tempText.decode("utf-8")
				#print link
				text = text.replace(link, tempText)
				
			return text
			
		return substitute(text)
		
		

class history:
	@db
	def GET(self, cursor):
		#if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
		#	raise web.seeother('/login')
		return env.get_template("base.html").render(history = giveHistory(cursor), allHistory = True)
		
class results:
	def GET(self, fileName):
		f = open(os.path.join("results", fileName))
		data = f.read()
		f.close()
		return data				

class clearHistory:
	@db
	def POST(self, cursor):	
		cursor.execute('delete from results where 1')
		return
		
		
class deleteRow:
	@db
	def POST(self, id, cursor):
		data = web.data()
		print 'Data given to delete'
		pprint(data)
		
		cursor.execute('delete from results where id = ?', (id, ))	

class Login:
	def GET(self):
		auth = web.ctx.env.get('HTTP_AUTHORIZATION')
		authreq = False
		if auth is None:
			authreq = True
		else:
			auth = re.sub('^Basic ','',auth)
			username,password = base64.decodestring(auth).split(':')
			if (username,password) in allowed:
				raise web.seeother('/')
			else:
				authreq = True
		if authreq:
			web.header('WWW-Authenticate','Basic realm="Please authorize yourself"')
			web.ctx.status = '401 Unauthorized'
			return		
	
def randomId(idLen = 7):	
	allChars = string.ascii_lowercase + "0123456789"
	return ''.join(random.choice(allChars) for i in xrange(idLen))

import smtplib
from email.mime.text import MIMEText
def email(**data):
	''' 
	Do note that while installing on a new server, 
	sendmail must be installed
	Install it using "sudo apt-get install sendmail-bin"
	'''
	
	msg = '''"Hey,
	 
Gocharts job: "%(jobName)s"  has been completed successfully. Please see the results by visiting the server. 

Your job id is: "%(jobId)s"

Thanks,
Your very own gocharts server ;)"
''' % data
	
	msg = MIMEText(msg)
	msg["Subject"] = "Gocharts job: " + data["jobName"] + " complete."
	msg["From"] = "gocharts@localhost"
	msg["To"] = data["email"]
	
	s = smtplib.SMTP("localhost")
	s.sendmail(msg["From"], [msg["To"]], msg.as_string())
	s.quit()	

application = app.wsgifunc()	
	
if __name__ == "__main__":


	#app.run()
	jobQueue = Queue.Queue()
	import __builtin__
	__builtin__.jobQueue = jobQueue
	
#	thread = threading.Thread(target = process, args = (jobQueue,), )
#	thread.start()
#	print "started the thread"
	
	#For handling Ctrl+C interrupts correctly
	import signal, sys
	def signal_handler(signal, frame):
		global thread, app
		print 'You pressed Ctrl+C!'
		jobQueue.put(None)
		app.stop()
		thread.exit()
		sys.exit(0)
	#signal.signal(signal.SIGINT, signal_handler)
	#signal.pause()
	#print "setup of Ctrl + C complete"
	
	
	app.run()
