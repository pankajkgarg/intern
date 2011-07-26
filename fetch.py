'Provides routines to fetch urls'
'''
LICENSE:
	Copyright (C) 2011, Pankaj Kumar Garg
	This program is distributed under GNU General Public License
'''


__author__  = "Pankaj Kumar Garg"
__email__   = "pankajn17@gmail.com"
__copyright__ = "Copyright (c) 2011, Pankaj Kumar Garg"
__license__ = "GPLv3"

import socket
import random, subprocess, shlex, cStringIO, gzip
import urllib, urllib2
from pprint import pprint


def fetch(url, postData = None, headers = {}, timeout = 120, debug = False, unarchive = False):
	'''
	Since, this module uses python library to fetch urls, there are few things to notice
	Proxy setup is automatically done by the python
	Behind proxy, "https" links cann't be opened
	To solve this, look at: http://code.activestate.com/recipes/456195/
	
	postData:  pass a dictionary
	decode is True by default
	responseCode shall be int
	headers returned shall be  a dict
	'''
	data = {}
	data["url"] = url
	
	socket.setdefaulttimeout(timeout)
	
	if "User-Agent" not in headers:
		headers["User-Agent"] = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
	#"Referer"
		
	if postData is not None:
		post = urllib.urlencode(postData)
	else:
		post = None
	
	req = urllib2.Request(url, post, headers)
	try:
		response = urllib2.urlopen(req)
		data["html"] = response.read()
	except urllib2.URLError, e:	
		if debug:
			if hasattr(e, 'reason'):
				print 'We failed to reach a server.'
				print 'Reason: ', e.reason
			elif hasattr(e, 'code'):
				print 'The server couldn\'t fulfill the request.'
				print 'Error code: ', e.code
		data["html"] = ""
		data["effectiveUrl"] = ""
		data["headers"] = {}
		if hasattr(e, "code"):
			data["responseCode"] = e.code
		else:
			data["responseCode"] = 444
	except Exception:
		#timeout error
		data["html"] = ""
		data["effectiveUrl"] = ""
		data["headers"] = {}
		data["responseCode"] = 444
	else:
		# everything is fine
		data["effectiveUrl"] = response.geturl()
		data["headers"] = dict(response.info())
		data["responseCode"] = 200
		
		headerStr = "\n".join(key + ": " + value for key,value in data["headers"].iteritems())
		
		pprint(data['headers'])
		#decoding
		if unarchive:
			if "content-encoding" in data["headers"] and data["headers"]["content-encoding"] == "gzip":
				try:
					htmlStream = cStringIO.StringIO(urlList[i]['html'][:])
					gzipper = gzip.GzipFile(fileobj = htmlStream, mode="rb")
					urlList[i]['html'] = gzipper.read()
					gzipper.close()
					htmlStream.close()
				
				except:
					urlList[i]['html'] = urlList[i]['html']
		
	return data
