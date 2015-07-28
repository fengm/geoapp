
def visualize_file(f, d_out):
	import os

	_f_out = os.path.join(d_out, os.path.basename(f[:-7]) + '_b321.tif')
	if os.path.exists(_f_out):
		return

	_cmd = 'visualize_bands.py -i %s -b 3 2 1 -sr -c -o %s' % (f, _f_out)
	os.system(_cmd)

def main():
	_opts = _init_env()

	# _d_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_3'
	# _d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_3_vis'
	_d_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_7'
	_d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_7_vis'

	_ps = []
	import os
	for _root, _dirs, _files in os.walk(_d_inp):
		for _file in _files:
			# if _file.endswith('_com.img.gz') and ('p026r039' in _file):
			if _file.endswith('_com.img.gz'):
				_ps.append((os.path.join(_root, _file), _d_out))

	os.path.exists(_d_out) or os.makedirs(_d_out)
	# _ps = _ps[:2]

	import multi_task
	multi_task.Pool(visualize_file, _ps, 1).run_single()

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

