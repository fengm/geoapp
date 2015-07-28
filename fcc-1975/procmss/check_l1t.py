

def check_image(f):
	import mod_mss
	import file_unzip

	with file_unzip.file_unzip('/export2/data/mfeng/tmp/mss_check') as _zip:
		_fs, _fm, _fo = mod_mss.load_file(f, _zip)
		_ms = mod_mss.parse_metadata(_fm)

		import pymongo
		with pymongo.MongoClient('129.2.12.64', 27017) as _c:
			_db = _c.tasks
			_rs = _db.mss

			import os, landsat
			_p = landsat.parse(os.path.basename(f))

			if _rs.find_one({'code': str(_p)}) != None:
				return

			_rs.insert_one({'file': f, 'type': _ms['DATA_TYPE'], 'code': str(_p)})

def main():
	_opts = _init_env()
	del _opts

	import pymongo
	_c = pymongo.MongoClient('129.2.12.64', 27017)
	# _c = pymongo.MongoClient('129.2.12.64', 18090)
	_db = _c.tasks
	_rs = _db.mss

	_rs.delete_many({})

	# print _rs.count()
	# print _rs.find_one_and_delete({})
	# print _rs
	# print _rs.insert_one({'name': 'good1'})

	# print dir(_rs)

	# _f_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/list/mss_1975_list.txt'
	_f_inp = '/data/glcf-nx-002/data/PALSAR/fcc_1975/list/mss_1975_list_gls.txt'
	_f_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/list/mss_1975_list_check_gls.txt'

	_ps = []
	for _f in open(_f_inp).read().splitlines():
		if _f.strip():
			_ps.append(_f)

	# _ps = _ps[:30]
	import multi_task
	multi_task.Pool(check_image, _ps, 12, True).run()

	print _rs.count()
	_ls = []
	for _r in _rs.find():
		_ls.append(','.join([_r[u'file'], _r[u'type']]))

	with open(_f_out, 'w') as _fo:
		_fo.write('\n'.join(_ls))

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

