#from pprint import pprint
#pprint(globals())
	
import web
from web import form
import goTree
from jinja2 import Environment, PackageLoader
from microarray import MicroArray
import re, Queue, threading, string, random, os
from pprint import pprint

#render = web.template.render('templates/')
env = Environment(loader=PackageLoader('__main__', 'templates'))


print __name__

urls = (
	'/', 'home',
)

app = web.application(urls, globals(), autoreload = None)

vemail = form.regexp(r"^$|.*@.*", "must be a valid email address")
vnums = form.regexp(r"[0-9]+\s*,?\s*", 'must be a list of comma separatad numbers e.g "2,7,12" ')
vnum = form.regexp(r"[0-9]+", "must be a column number")

vfile = form.regexp(r".+", "Required!")

inputForm = form.Form(
	form.File("annotationFile", vfile, description = 'Annotation File'),
	form.File('ontologyFile', vfile, description = 'Ontology File'),
	
	form.File('microarrayFile', vfile, description = 'Microarray file (Only cvs files)'),
	form.Textbox("ignoreRows", vnums, description = "Rows in microarray to ignore"),
	form.Textbox("ignoreCols", vnums, description = "Cols in microarray to ignore"),
	form.Textbox("idCol", vnum, description = "Column containing the Go term Id"),
	
	form.Dropdown("evidenceCodes", args=["all", "exp"], value="all", description="Evidence Codes"),
	
	form.Textbox("jobName", description = "Job Name (Optional)"),
	form.Textbox('email', vemail, description = 'Email (optional)'),
	form.Button("submit", type="submit", description="Submit"),
)

def process(jobQueue):
	while True:
		kwargs = jobQueue.get()
		
		if kwargs == None:
			print 'stopping the process'
			break	
	
		if 'ignoreRows' in kwargs:
			kwargs['ignoreRows'] = re.sub(r'\s+', '', kwargs['ignoreRows']).split(',')
			kwargs['ignoreRows'] = [int(temp) for temp in kwargs['ignoreRows']]
		else:
			kwargs['ignoreRows'] = []
		
		if 'ignoreCols' in kwargs:
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

		temp = goTree.GoTree(kwargs['ontologyFile'], kwargs['annotationFile'], microarray, name = kwargs['jobId'], evidenceCodes = kwargs['evidenceCodes'])
	
		if kwargs["email"] is not None:
			email(**kwargs)
	
			
			
class home:
	def GET(self):
		f = inputForm()
		
		formText = f.render()
		
		return env.get_template("home.html").render(form = formText, )
		#return render.register(f)
	def POST(self):
		f = inputForm()
		if not f.validates():
			formText = f.render()
			return env.get_template("home.html").render(form = formText, )
		else:
			dataDict = dict(f.value)
			dataDict["jobId"] = randomId()
			
			ontologyPath = os.path.join("data", "ontologyFile.obo")
			f = open(ontologyPath, "wb")
			f.write(dataDict["ontologyFile"])
			f.close()
			dataDict["ontologyFile"] = ontologyPath
			
			annotationPath = os.path.join("data", "annotations.txt")
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
			 
			jobQueue.put(dataDict)
			return "Done"
	
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
	
	thread = threading.Thread(target = process, args = (jobQueue,), )
	thread.start()
	print "started the thread"
	
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
	print "setup of Ctrl + C complete"
	
	app.run()
