import os, sys
curdir = os.path.dirname(__file__)

if curdir not in sys.path:
    sys.path.append(curdir)

from server import app
import web




sys.stdout = sys.stderr

curdir = os.path.dirname(__file__)
session = web.session.Session(app, web.session.DiskStore(os.path.join(curdir,'sessions.abrad')),)

application = app.wsgifunc()
