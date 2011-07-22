#from pprint import pprint
#pprint(globals())
	
import web
from web import form
import goTree
from jinja2 import Environment, PackageLoader
from microarray import MicroArray
import re, Queue, threading, string, random, os, sqlite3, json, time, base64
from pprint import pprint

#render = web.template.render('templates/')
env = Environment(loader=PackageLoader('__main__', 'templates'))

allowed = (
    ('admin','biology'),
)


urls = (
	'/', 'home',
	'/login', 'Login',
	'/view/(?P<jobId>[^/]+)', 'view',
	'/results/(?P<fileName>[^/]+)', 'results',
	'/history', 'history',
	'/deleteAll', "clearHistory",
	'/deleteRow/(?P<id>[^/]+)', 'deleteRow',
)

web.config.debug = True
app = web.application(urls, globals())
 
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
	

vemail = form.regexp(r"^$|.*@.*", "must be a valid email address")
vnums = form.regexp(r"[0-9]*\s*,?\s*", 'must be a list of comma separatad numbers e.g "2,7,12" ')
vnum = form.regexp(r"[0-9]+", "must be a column number")

vfile = form.regexp(r".+", "Required!")

inputForm = form.Form(
	form.File("annotationFile", vfile, description = 'Annotation File'),
	#form.File('ontologyFile', vfile, description = 'Ontology File'),
	
	form.File('microarrayFile', vfile, description = 'Microarray file (Only csv files)'),
	form.Textbox("ignoreRows", vnums, description = "Rows in microarray to ignore"),
	form.Textbox("ignoreCols", vnums, description = "Cols in microarray to ignore"),
	form.Textbox("idCol", vnum, description = "Column containing the Go term Id"),
	
	# form.Checkbox("treeTypes", value="biological_process", description="biological process", checked=True),
	# form.Checkbox("treeTypes", description="molecular function", value="molecular_function"),
	# form.Checkbox("treeTypes", description="cellular component", value="cellular_component"),
	form.Dropdown("treeTypes", args=["biological_process", "molecular_function", "cellular_component"], value="biological_process", description="Tree type"),
	
	form.Dropdown("evidenceCodes", args=["all", "exp"], value="all", description="Evidence Codes"),
	
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
		
	if 'evidenceCodes' in kwargs:
		if kwargs['evidenceCodes'] == 'all':
			kwargs['evidenceCodes'] = []
	
	microarray = MicroArray(kwargs['microarrayFile'], idCol = kwargs['idCol'], ignoreRows = kwargs['ignoreRows'], ignoreCols = kwargs['ignoreCols'] )

	#Shall do the processing and shall save the results in 'results' directory 
	#	with the filename  "kwargs['name']" + ".js"

	temp = goTree.GoTree(kwargs['ontologyFile'], kwargs['annotationFile'], microarray, treeTypes = [kwargs["treeTypes"]], uid = kwargs['jobId'], name = kwargs['jobName'], evidenceCodes = kwargs['evidenceCodes'])
	
	
	#os.remove(kwargs['ontologyFile'])
	os.remove(kwargs['annotationFile'])
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
		if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
			raise web.seeother('/login')

		f = inputForm()
		
		formText = f.render()
		return env.get_template("home.html").render(form = formText, history = giveHistory(cursor, 5))
		#return render.register(f)
		
	@db	
	def POST(self, cursor):
		#print type(web.data())
		#pprint(web.data())
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
			
			dataDict["ontologyFile"] = os.path.join("data", "gene_ontology.obo")
			
			annotationPath = os.path.join("data", dataDict["jobId"] + "_annotations.txt")
			f = open(annotationPath, "wb")
			f.write(dataDict["annotationFile"])
			f.close()
			dataDict["annotationFile"] = annotationPath
			
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
	def GET(self, jobId, cursor):
		if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
			raise web.seeother('/login')
		cursor.execute("select * from results where id = ?", (jobId, ))
		info = dict(cursor.fetchone())
		
		info['dataFile'] = '/results/' +  info['id'] + '_data.js'
		info['zoomLabelsFile'] = '/results/' + info['id'] + '_zoomLabels.js'
		
		return env.get_template("view.html").render(info = info, history = giveHistory(cursor, 5))

class history:
	@db
	def GET(self, cursor):
		if web.ctx.env.get('HTTP_AUTHORIZATION') is None:
			raise web.seeother('/login')
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
			web.header('WWW-Authenticate','Basic realm="Auth example"')
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
	
	dbSetup()
	app.run()
