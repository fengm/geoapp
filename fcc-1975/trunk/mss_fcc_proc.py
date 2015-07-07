'''
File: mss_fcc_proc.py
Author: Min Feng
Version: 0.1
Create: 2015-07-06 13:17:37
Description:
'''

import logging

def process_mss(f_inp, d_out):
	import mod_mss_fcc
	import file_unzip

	with file_unzip.file_unzip() as _zip:
		_d_out = _zip.generate_file('out')
		mod_mss_fcc.mss_fcc(f_inp, _d_out, _zip)

		import config
		if not config.cfg.getboolean('conf', 'debug'):
			# clear the folder
			import os
			for _f in os.listdir(_d_out):
				_p = _f.split('_')[-1]
				if _p not in ['dat.tif', 'err.tif', 'txt.tree', 'txt.names', 'txt.data']:
					logging.info('clean file: %s' % os.path.join(_d_out, _f))
					try:
						os.remove(os.path.join(_d_out, _f))
					except Exception:
						pass

		file_unzip.compress_folder(_d_out, d_out, ['.txt', '.tmp', '.cases', '.data'])

def process_file(f_inp, d_out):
	import config

	if config.cfg.getboolean('conf', 'debug'):
		process_mss(f_inp, d_out)
		return

	try:
		process_mss(f_inp, d_out)
	except KeyboardInterrupt, err:
		print '\n\n* User stopped the program'
		raise err
	except Exception, err:
		import traceback

		logging.error(traceback.format_exc())
		logging.error(str(err))

		print '\n\n* Error:', err

def main():
	_opts = _init_env()
	del _opts

	_f_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_6/p026r039/p026r039_com.img.gz'
	_d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/forest_6/test1'

	process_file(_f_inp, _d_out)

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

