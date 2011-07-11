import web

urls = (
	'/', 'home',
)

app = web.application(urls, globals())

class home:
	def GET(self):
		
		pass
		
if __name__ == "__main__":
	app.run()
