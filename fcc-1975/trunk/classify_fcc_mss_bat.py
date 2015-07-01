'''
File: classify_fcc_mss_bat.py
Author: Min Feng
Version: 0.1
Create: 2013-12-15 16:17:14
Description:
'''

def process_landsat_scene(f_inp, f_out):
	import logging

	try:
		import config
		_d_tmp = config.cfg.get('path', 'path_temp')

		import file_unzip
		import os
		with file_unzip.file_unzip(_d_tmp) as _zip:

			_d_tmp = _zip.generate_file()
			os.makedirs(_d_tmp)

			import classify_fcc

			_d_inp = classify_fcc.unzip_input_scene(f_inp, _zip)
			if _d_inp == None:
				return

			classify_fcc.classify_forest_mss(_d_inp, _d_tmp, _zip)

			file_unzip.compress_folder(_d_tmp, f_out)

	except KeyboardInterrupt:
		print '\n\n* User stopped the program'
	except Exception, err:
		import traceback

		logging.error('%s, %s, %s'  % (f_inp, f_out, traceback.format_exc()))
		logging.error('%s, %s, %s'  % (f_inp, f_out, str(err)))

		print '\n\n* Error:', err

def main():
	_opts = _init_env()

	import multi_task
	import config

	_f_inp = config.cfg.get('path', 'input')
	_f_out = config.cfg.get('path', 'output')

	print 'input', _f_inp
	print 'output', _f_out

	_fs = multi_task.load_list(_f_inp, _opts.instance_num, _opts.instance_pos)

	import os
	_ts = []
	for _f in _fs:
		_f_ot = os.path.join(_f_out, os.path.basename(os.path.normpath(_f)))
		if os.path.exists(_f_ot):
			continue
		_ts.append((_f, _f_ot))
		# break

	import config
	_d_tmp = config.cfg.get('path', 'path_temp')

	import os
	if os.path.exists(_d_tmp):
		try:
			import shutil
			shutil.rmtree(_d_tmp)
		except Exception:
			pass

	_pool = multi_task.Pool(process_landsat_scene, _ts, _opts.task_num)
	_pool.run()

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

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

