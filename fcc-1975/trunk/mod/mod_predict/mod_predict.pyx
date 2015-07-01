import numpy as np
import logging
import config
import run_commands
import re

cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def predict_block(cmd, f_case, np.ndarray[np.uint8_t, ndim=2] dat_out,
		np.ndarray[np.float32_t, ndim=2] dat_err, locs):
	logging.debug('predict block %s, %s' % (len(locs), f_case))

	_forest = not config.cfg.getboolean('config', 'percent_forest')

	_f_mod = f_case[:len(f_case) - 6]
	_cmd = cmd + ' -f ' + _f_mod + ('' if _forest else ' -e ')

	_ls = run_commands.run_exe(_cmd)[0].splitlines()[3:]
	cdef int _nu = 0

	_reg = '^\s*\d+\s+\S+\s+(\S+)\s+\[(.+)\]\s*$' if _forest else \
			'\S?\s*-?\d*\.?\d*\s+\S+\s+(-?\d*\.?\d*)\s*\+\-\s*(-?\d*\.?\d*)'

	logging.debug('predicted lines %s' % len(_ls))
	for _l in _ls:
		_m = re.match(_reg, _l)
		if _m:
			dat_out[locs[_nu][0], locs[_nu][1]] = int(float(_m.group(1)))
			dat_err[locs[_nu][0], locs[_nu][1]] = float(_m.group(2))
			_nu += 1

	logging.info('predicted %d values (%d)' % (_nu, len(locs)))
	if _nu != len(locs):
		logging.error('number of the predict values (%s) does not match the expected values (%s)\ncmd: %s' % (_nu, len(locs), _cmd))
		raise Exception('number of the predict values (%s) does not match the expected values (%s)' % (_nu, len(locs)))
	assert(_nu == len(locs))

def predict(cmd, bnd_qa, bnds, f_case, f_dat, f_err, fzip):
	logging.debug('start building')
	logging.info('case file: %s' % f_case)

	cdef np.ndarray[np.uint8_t, ndim=2] _dat_out
	cdef np.ndarray[np.float32_t, ndim=2] _dat_err
	cdef int _row, _col
	cdef float _v_slp
	cdef int _v, _nodata

	_forest = not config.cfg.get('config', 'percent_forest')
	_nodata = 0 if _forest else 255

	import numpy
	_dat_out = numpy.zeros((bnd_qa.height, bnd_qa.width), dtype=numpy.uint8)
	_dat_out.fill(_nodata)
	_dat_err = numpy.empty((bnd_qa.height, bnd_qa.width),
			dtype=numpy.float32)
	_dat_err.fill(-9999)

	_bs = bnds.keys()
	_bs.sort()

	import progress_percentage
	_ppp = progress_percentage.progress_percentage(bnd_qa.height)

	import os
	_max_cache_pixels = config.cfg.getint('predict', 'max_cache_pixels')
	_cache_lines = config.cfg.getint('predict', 'cache_lines')

	cdef np.ndarray[np.uint8_t, ndim=2] _dat
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_qa = bnd_qa.data
	cdef int _rrr, _row_r

	_o_data = open(f_case, 'w')
	_loc_s = []

	import geo_raster_c as ge
	for _rrr in xrange(0, bnd_qa.height, _cache_lines if _cache_lines > 0 else bnd_qa.height):
		_ppp.next(_cache_lines)

		_rows = min(bnd_qa.height - _rrr, _cache_lines)
		_dats = []
		for _b in _bs:
			if isinstance(bnds[_b], ge.geo_band):
				_dats.append(bnds[_b].read_rows(_rrr, _rows))
			else:
				_dats.append(bnds[_b].data[_rrr: _rrr + _rows, :])

		for _row in xrange(_rows):
			_row_r = _rrr + _row

			for _col in xrange(bnd_qa.width):
				_v_qa = _dat_qa[_row + _rrr, _col]

				if _v_qa != 1:
					continue

				_vs = ['?', '1']
				for _b in xrange(len(_bs)):
					_dat = _dats[_b]
					_v = _dat[_row, _col]
					_vs.append('%s' % _v)

				# _o_data.write(','.join(['%s' % _v for _v in _vs]))
				_o_data.write(','.join(_vs))
				_o_data.write('\n')

				_loc_s.append((_row_r, _col))

		if len(_loc_s) >= _max_cache_pixels:
			_o_data.flush()
			_o_data.close()
			predict_block(cmd, os.path.abspath(f_case), _dat_out, _dat_err, _loc_s)

			_o_data = open(f_case, 'w')
			_loc_s = []

	if len(_loc_s) > 0:
		_o_data.flush()
		_o_data.close()
		predict_block(cmd, os.path.abspath(f_case), _dat_out, _dat_err, _loc_s)

	# remove the temperary file
	del _loc_s

	_ppp.done()

	# output predicted results
	import geo_raster_c as ge
	logging.info('output image file: %s' % f_dat)

	if not _forest:
		_dat_out[(_dat_out != _nodata) & (_dat_out > 100)] = 100

	logging.debug('output predicted water grids')
	_bnd_dat = ge.geo_band_cache(_dat_out, bnd_qa.geo_transform,
			bnd_qa.proj, _nodata, ge.pixel_type())
	_bnd_dat.save(f_dat)

	logging.info('output error file: %s' % f_err)
	_bnd_err = ge.geo_band_cache(_dat_err, bnd_qa.geo_transform,
			bnd_qa.proj, -9999, ge.pixel_type('float'))
	_bnd_err.save(f_err)

