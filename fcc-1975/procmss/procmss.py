
'''
File: procmss.py
Author: Min Feng
Version: 0.1
Create: 2015-06-16 16:17:37
Description:
'''

import logging

def process_file(f_inp, f_out):
	import mod_mss
	import file_unzip

	with file_unzip.file_unzip() as _zip:
		import landsat
		import os

		_p = landsat.parse(os.path.basename(f_inp))
		if _p == None:
			raise Exception('failed to parse file name')

		_d_tmp = _zip.generate_file()
		os.makedirs(_d_tmp)

		mod_mss.proc_mss(f_inp, os.path.join(_d_tmp, os.path.basename(f_out)), _zip)
		file_unzip.compress_folder(_d_tmp, os.path.dirname(f_out))

def main():
	_opts = _init_env()

	import os
	import config
	import multi_task

	_f_inp = _opts.input if _opts.input else config.cfg.get('conf', 'input')
	logging.info('input: %s' % _f_inp)

	_fs = [_f_inp] if _f_inp.endswith('.tar.gz') else multi_task.load_from_list(_f_inp, _opts)

	_d_out = _opts.output if _opts.output else config.cfg.get('conf', 'output')
	logging.info('output: %s' % _d_out)
	os.path.exists(_d_out) or os.makedirs(_d_out)

	_ps = []
	for _f in _fs:
		_f_out = os.path.join(_d_out, os.path.basename(_f)[:-6] + 'img')
		if any(map(os.path.exists, [_f_out, _f_out + '.gz'])):
			logging.info('skip processed %s' % _f)
			continue
		_ps.append((_f, _f_out))

	multi_task.Pool(process_file, _ps, _opts.task_num, True).run()

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	_p.add_argument('-i', '--input', dest='input')
	_p.add_argument('-o', '--output', dest='output')

	import multi_task
	multi_task.add_task_opts(_p)

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

