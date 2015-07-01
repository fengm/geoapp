'''
File: classify_fcc_mss.py
Author: Min Feng
Version: 0.1
Create: 2013-11-26 17:31:44
Description: run the model for classifying forest cover change from MSS
'''

def main():
	_opts = init_env()

	import config
	_d_tmp = config.cfg.get('path', 'path_temp')

	import file_unzip
	with file_unzip.file_unzip(_d_tmp) as _zip:

		import os
		os.path.exists(_opts.output) or os.makedirs(_opts.output)

		import classify_fcc

		_d_inp = classify_fcc.unzip_input_scene(_opts.input, _zip)
		if _d_inp == None:
			return

		import classify_fcc
		classify_fcc.classify_forest_mss(_d_inp, _opts.output, _zip)

def usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('-i', '--input', dest='input', required=True)
	_p.add_argument('-o', '--output', dest='output', required=True)

	return _p.parse_args()

def init_env():
	import os, sys
	_d_in = os.path.join(sys.path[0], 'lib')
	if os.path.exists(_d_in):
		sys.path.append(_d_in)

	_opts = usage()

	import logging_util
	logging_util.init(_opts.logging)

	import config
	config.load()

	return _opts

if __name__ == '__main__':
	main()

