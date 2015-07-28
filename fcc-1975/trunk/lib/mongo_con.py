'''
File: mongo_con.py
Author: Min Feng
Version: 0.1
Create: 2015-07-06 17:18:40
Description: establish the mongodb connection
'''

class connect:

	def __init__(self):
		import pymongo
		self._con = pymongo.MongoClient('129.2.12.64', 27017)

	def __enter__(self):
		return self._con

	def __exit__(self, type, value, traceback):
		if self._con:
			self._con.close()

