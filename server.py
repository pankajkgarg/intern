import web
from web import form
import goTree
from microarray import MicroArray
import re

render = web.template.render('templates/')

urls = (
	'/', 'home',
)

app = web.application(urls, globals())

vemail = form.regexp(r".*@.*", "must be a valid email address")
vnums = form.regexp(r"[0-9]\s*,\s*", 'must be a list of comma separatad numbers e.g "2,7,12" ')
vnum = form.regexp(r"[0-9]+", "must be a column number")

inputForm = form.Form(
	form.File(name = "annotationFile", description = 'Annotation File'),
	form.File(name = 'ontologyFile', description = 'Ontology File'),
	
	form.File(name='microarrayFile', description = 'Microarray file (Only cvs files)'),
	form.Textbox("ignoreRows", vnums, description = "Rows in microarray to ignore"),
	form.Textbox("ignoreCols", vnums, description = "Cols in microarray to ignore"),
	form.Textbox("idCol", vnum, description = "Column containing the GoId"),
	
	form.Dropdown("evidenceCodes", args=["all", "exp"], value="all", description="Evidence Codes"),
	
	form.Textbox("jobName", description = "Job Name (Optional)"),
	form.Textbox('email', vemail, description = 'Email (optional)'),
	form.Button("submit", type="submit", description="Submit"),
)

def process(microarrayFile, **kwargs):
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
	
	
	microarray = MicroArray(kwargs['microarrayFile'], idCol = kwargs['idCol'], ignoreRows = kwargs['ignoreRows'], ignoreCols = kwargs['ignoreCols'] )
	
	#Shall do the processing and shall save the results in 'results' directory 
	#	with the filename  "kwargs['name']" + ".js"
	temp = GoTree(kwargs['ontologyFile'], kwargs['annotationFile'], kwargs, microarray, name = kwargs['id'], evidenceCodes = kwargs['evidenceCodes'])
	
	if kwargs["email"] is not None:
		email(to = kwargs["email"], jobId = kwargs["id"], jobName = kwargs["jobName"])
	
			

class home:
	def GET(self):
		f = inputForm()
		
		formText = f.render()
		
		return render.register(f)
	def POST(self):
		f = inputForm()
		if not f.validates():
			formText = f.render()
			return render.register(f)
		else:
			
		
		
if __name__ == "__main__":
	app.run()
