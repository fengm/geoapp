'''
File: mss_composite.py
Author: Min Feng
Version: 0.1
Create: 2015-06-25 14:43:21
Description:
'''

import logging

def save_bands(f_out, bnds):
	import geo_raster_c as ge
	import os

	(lambda x: os.path.exists(x) or os.makedirs(x))(os.path.dirname(f_out))

	_bnd = bnds[0]
	ge.geo_raster.create(f_out, [len(bnds), _bnd.height, _bnd.width], _bnd.geo_transform,
			_bnd.proj, driver='HFA', nodata=_bnd.nodata, pixel_type=_bnd.pixel_type)

	_img = ge.open(f_out, True)
	for _b in xrange(len(bnds)):
		_img.get_band(_b+1).write(bnds[_b].data, 0, 0)
	_img.flush()

def composite(fs, f_out, fzip):
	import geo_raster_c as ge
	import mod_grid

	_imgs = [ge.open(fzip.unzip(_f)) for _f in fs]
	_bnds, _bnd_ndvi, _bnd_ndwi, _bnd_indx, _bnd_lyrs = mod_grid.combine(_imgs)

	save_bands(f_out, _bnds)
	_bnd_ndvi.save(f_out[:-4] + '_ndvi.img')
	_bnd_ndwi.save(f_out[:-4] + '_ndwi.img')
	_bnd_indx.save(f_out[:-4] + '_indx.img')
	_bnd_lyrs.save(f_out[:-4] + '_lyrs.img')

	with open(f_out[:-4] + '_lyrs.txt', 'w') as _fo:
		_fo.write('\n'.join(['%s:%s' % (_b, fs[_b]) for _b in xrange(len(fs))]))

def load_tiles_tile(t):
	import re
	_m = re.match('p(\d+)r(\d+)', t)
	_p = int(_m.group(1))
	_r = int(_m.group(2))

	return load_tiles(_p, _r)

def filter_scenes(fs):
	import config
	import landsat
	import os

	with open(config.cfg.get('conf', 'exclude_list')) as _fi:
		_ls = [str(landsat.parse(_l.strip())) for _l in _fi.read().splitlines() if _l.strip()]
		_fs = []
		for _f in fs:
			if str(landsat.parse(os.path.basename(_f))) in _ls:
				logging.info('skip scene %s' % _f)
				continue

			_fs.append(_f)
		return _fs

def load_tiles(p, r):
	import pymongo
	with pymongo.MongoClient('129.2.12.64', 27017) as _c:
		_db = _c.tasks
		_rs = _db.mss_l1t

		_fs = filter_scenes([_o['file'] for _o in _rs.find({'$and': [{'tile': 'p%03dr%03d' % (p, r)}, \
				{'l1_type': 'L1T'}]})])
		if len(_fs) == 0:
			logging.warning('failed to find any images at p%03dr%03d' % (p, r))
			return []

		for _y in xrange(-1, 2):
			for _x in xrange(-1, 2):
				if _y == 0 and _x == 0:
					continue

				_fs.extend(filter_scenes([_o['file'] for _o in _rs.find({'tile': 'p%03dr%03d' % (p + _y, r + _x), \
						'l1_type': 'L1T'})]))

		return _fs

def process_tile(tile, d_out):
	import file_unzip
	import os
	with file_unzip.file_unzip('/export2/data/mfeng/tmp/mss_ndvi') as _zip:
		_fs = load_tiles_tile(tile)
		if len(_fs) == 0:
			return

		# print 'found', len(_fs)
		_d_tmp = _zip.generate_file()
		os.makedirs(_d_tmp)
		composite(_fs, os.path.join(_d_tmp, '%s_com.img' % tile), _zip)
		file_unzip.compress_folder(_d_tmp, os.path.join(d_out, tile))

def main():
	_opts = _init_env()

	import multi_task

	# _tile = 'p026r039'

	_d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_7'
	_f = '/data/glcf-nx-002/data/PALSAR/fcc_1975/src/combine/conf/test_tiles.txt'

	_ps = []

	if False:
		for _t in multi_task.load_from_list(_f, _opts):
			if not _t:
				continue
			_ps.append((_t, _d_out))

	# _ps = [('p026r039', _d_out)]
	_ps = [('p009r024', _d_out)]
	multi_task.Pool(process_tile, _ps, 8, True).run()

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


