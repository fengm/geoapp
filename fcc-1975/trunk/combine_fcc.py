'''
File: combine_fcc.py
Author: Min Feng
Version: 0.1
Create: 2015-07-23 17:55:03
Description:
'''

def update(bs1, bs2):
	bs1[1].data[bs1[0].data == 2] = 1.0
	bs2[1].data[bs2[0].data == 2] = 1.0

	_dat = (bs2[1].data > bs1[1].data) & (bs1[0].data != bs1[0].nodata)

	import mod_grid
	mod_grid.update(bs1[0].data, _dat, bs2[0].data)
	mod_grid.update(bs1[1].data, _dat, bs2[1].data)

def combine(fs, f_clr, d_out):
	if len(fs) == 0:
		return

	import os

	print 'target', fs[0]
	_bnd, _err = load_img(fs[0])

	import progress_percentage
	_ppp = progress_percentage.progress_percentage(len(fs) - 1)

	for _i in xrange(1, len(fs)):
		_ppp.next(count=True, message=os.path.basename(fs[_i]))
		update((_bnd, _err), load_img(fs[_i], _bnd))

	_ppp.done()

	import landsat
	_inf = landsat.parse(fs[0])

	os.path.exists(d_out) or os.makedirs(d_out)

	_f_out = os.path.join(d_out, '%s_com_dat.tif' % _inf.tile)
	_f_err = os.path.join(d_out, '%s_com_err.tif' % _inf.tile)

	print 'write 1'
	_bnd.save(_f_out, color_table=load_color(f_clr), opts=['compress=lzw'])
	_err.save(_f_err, opts=['compress=lzw'])

def load_color(f):
	import geo_raster_c as ge
	return ge.load_colortable(f)

	# _bnd = ge.open(f).get_band()
	# return _bnd.color_table

def load_img(f, bnd=None):
	import geo_raster_c as ge

	if bnd == None:
		return ge.open(f).get_band().cache(), \
				ge.open(f.replace('_dat.tif', '_err.tif')).get_band().cache()
	else:
		return ge.open(f).get_band().read_block(bnd), \
			ge.open(f.replace('_dat.tif', '_err.tif')).get_band().read_block(bnd)

def load_list(f):
	with open(f) as _fi:
		return filter(lambda x: x, [_f.strip() for _f in _fi.read().splitlines()])

def main():
	_opts = _init_env()

	combine(load_list(_opts.input), _opts.color, _opts.output)

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	_p.add_argument('-i', '--input', dest='input', required=True)
	_p.add_argument('-c', '--color', dest='color', \
			default='/data/glcf-st-004/data/workspace/fengm/prog/fcc_1975/v2/conf/colors/colors_dat.txt')
	_p.add_argument('-o', '--output', dest='output', required=True)

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

