
def load_list(fs):
	import landsat, os

	_rs = []

	for _f in fs:
		with open(_f) as _fi:
			for _l in _fi:
				_p = landsat.parse(os.path.basename(_l))
				_rs.append(str(_p))

	return _rs

def main():
	_opts = _init_env()

	_fs = ['/data/glcf-nx-002/data/PALSAR/fcc_1975/list/mss_1975_list_gls.txt',
			'/data/glcf-nx-002/data/PALSAR/fcc_1975/list/mss_1975_list.txt']

	_f_out = '/export2/data/mfeng/mss_1975_l1t_download_non.txt'
	_f_ext = '/export2/data/mfeng/mss_1975_l1t_download_non_existed.txt'

	_ts = load_list(_fs)
	_ts = []

	import psycopg2
	with psycopg2.connect("dbname='glcf_v5' user='glcfapp' host='glcfdb03.umd.edu' password='V!L!V!'") as _con:
		_cur = _con.cursor()
		_cur.execute("SELECT scene_id from data_granule where isl1t=FALSE")

		import landsat
		with open(_f_out, 'w') as _fo, open(_f_ext, 'w') as _fe:
			for _r in _cur:
				_p = landsat.parse(_r[0].strip())
				_l = str(_p)
				if _l not in _ts:
					_fo.write(_r[0].strip() + '\n')
				else:
					_fe.write(_r[0].strip() + '\n')

		# print dir(_cur)
		# _res = _cur.fetchall()
		# print _res


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

