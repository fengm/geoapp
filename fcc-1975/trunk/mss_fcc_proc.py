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
		import config
		if config.cfg.getboolean('conf', 'cache_output'):
			mod_mss_fcc.mss_fcc(f_inp, d_out, _zip)
		else:
			_d_out = _zip.generate_file('out')
			mod_mss_fcc.mss_fcc(f_inp, _d_out, _zip)
			file_unzip.compress_folder(_d_out, d_out, [])

def main():
	_opts = _init_env()

	# _tile = 'p026r039'
	# _tile = 'p240r078'

	# _f_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_6/%s/%s_com.img.gz' % (_tile, _tile)
	# _d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/forest_6/test_%s' % _tile

	import multi_task
	import os

	if _opts.input.endswith('.txt'):
		_fs = multi_task.load_from_list(_opts.input, _opts)
	else:
		_fs = [_opts.input]

	_ps = []
	for _f in _fs:
		_ps.append((_f, os.path.join(_opts.output, os.path.basename(_f).split('.')[0])))

	multi_task.Pool(process_mss, _ps, _opts.task_num).run()

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	_p.add_argument('-i', '--input', dest='input', required=True)
	_p.add_argument('-o', '--output', dest='output', required=True)

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

