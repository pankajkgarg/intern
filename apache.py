from server.py import app
import web

curdir = os.path.dirname(__file__)
session = web.session.Session(app, web.session.DiskStore(os.path.join(curdir,'sessions')),)

import sys
sys.stdout = sys.stderr



application = app.wsgifunc()
