
def main():
	_opts = _init_env()

	import os
	import landsat
	import pymongo
	with pymongo.MongoClient('129.2.12.64', 27017) as _c:
		_db = _c.tasks
		_rs = _db.mss

		for _r in _rs.find():
			print _r
			# _p = landsat.parse(os.path.basename(_r['file']))
			# _rs.update_one({'_id': _r['_id']}, {'$set': {'code': str(_p)}})
			break

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	return _p.parse_args()

def _init_env():
	import os, sys

	_dirs = ['lib', 'libs']
	_d_ins = [os.path.join(sys.path[0], _d) for _d in _dirs if \
			os.path.exists(os.path.join(sys.path[0], _d))]
	sys.path = [sys.path[0]] + _d_ins + sys.path[1:]

	_opts = _usage()

	import logging_util
	logging_util.init(_opts.logging)

	import config
	config.load(_opts.config)

	import file_unzip as fz
	fz.clean(fz.default_dir(_opts.temp))

	return _opts

if __name__ == '__main__':
	main()

