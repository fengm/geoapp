'''
File: mss_fcc_proc.py
Author: Min Feng
Version: 0.1
Create: 2015-07-06 13:17:37
Description:
'''

import logging

def filter_scenes(fs):
	return fs

	# import config
	# import landsat
	# import os

	# with open(config.cfg.get('conf', 'exclude_list')) as _fi:
	# 	_ls = [str(landsat.parse(_l.strip())) for _l in _fi.read().splitlines() if _l.strip()]
	# 	_fs = []
	# 	for _f in fs:
	# 		if str(landsat.parse(os.path.basename(_f))) in _ls:
	# 			logging.info('skip scene %s' % _f)
	# 			continue

	# 		_fs.append(_f)
	# 	return _fs

def load_images(t, add_next=False):
	import mongo_con

	with mongo_con.connect() as _c:
		_db = _c.tasks
		_rs = _db.mss_l1t

		print 'search tile %s' % t

		_p = int(t[1:4])
		_r = int(t[5:8])

		_fs = filter_scenes([_o['file'] for _o in _rs.find({'$and': [{'tile': 'p%03dr%03d' % (_p, _r)}, \
				{'l1_type': 'L1T'}]})])
		if len(_fs) == 0:
			logging.warning('failed to find any images at p%03dr%03d' % (_p, _r))
			return []

		if add_next:
			for _y in xrange(-1, 2):
				for _x in xrange(-1, 2):
					if _y == 0 and _x == 0:
						continue

					print 'search neighor p%03dr%03d' % (_p + y, _r + x)
					_fs.extend(filter_scenes([_o['file'] for _o in _rs.find({'tile': 'p%03dr%03d' % (_p + _y, _r + _x), \
							'l1_type': 'L1T'})]))

		print 'find', len(_fs), 'images'

		return _fs

def main():
	_opts = _init_env()

	with open(_opts.output, 'w') as _fo:
		_fo.write('\n'.join(load_images(_opts.tile, _opts.neighor)))

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	_p.add_argument('-t', '--tile', dest='tile', required=True)
	_p.add_argument('-n', '--neighor', dest='neighor', action='store_true')
	_p.add_argument('-o', '--output', dest='output', required=True)

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

