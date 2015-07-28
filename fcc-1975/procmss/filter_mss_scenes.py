
def main():
	_opts = _init_env()

	# _f = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/list/global_mss_1975_sr_ic_0626.txt'
	# _f_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/list/global_mss_1975_sr_ic_0626_l1t.txt'
	_f_out = '/a/glcffs02/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/list/global_mss_1975_sr_gls1975_0628.txt'

	import pymongo
	with pymongo.MongoClient('129.2.12.64', 27017) as _c:
		_db = _c.tasks
		_rs = _db.mss_l1t
		_gs = _db.mss


		# _rs.update_many({}, {'$set': {'l1_type': 'L1T'}})
		print _rs.find_one()

		# _ids = []
		for _o in _rs.find({'note': 'gls1975'}):
			if _gs.find_one({'code': _o['code'], 'type': 'L1G'}) != None:
				# _ids.append(_o['_id'])
				print _o
				break
				# _rs.update_one(_o, {'$set': {'l1_type': 'L1G'}})

		# print len(_ids)

		# _nu = 0
		# import landsat, os
		# for _l in open(_f_out).read().splitlines():
		# 	_r = landsat.parse(os.path.basename(_l))

		# 	if _gs.find_one({'code': str(_r), 'type': 'L1G'}) != None:
		# 		_nu += 1
				# print _l

			# if _o != None:
			# 	_nu += 1
			# 	print _l

			# _o = {'code': str(_r), 'tile': _r.tile, 'date': _r.ac_date, 'sensor': _r.sensor,
			# 		'mission': _r.mission, 'file': _l, 'note': 'gls1975'}

			# _rr = _rs.find_one({'code': str(_r)})
			# if _rr == None:
			# 	_nu += 1

			# 	# print _l
			# 	# print _rr
			# 	# break
			# 	_rs.insert(_o)

		# print _nu

	# 	_rr = []
	# 	import landsat, os

	# 	_ls = open(_f).read().splitlines()
	# 	for _l in _ls:
	# 		_r = landsat.parse(os.path.basename(_l))
	# 		if _rs.find_one({'code': str(_r)})['type'] == 'L1T':
	# 			_rr.append(_l)

	# 	with open(_f_out, 'w') as _fo:
	# 		_fo.write('\n'.join(_rr))

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

