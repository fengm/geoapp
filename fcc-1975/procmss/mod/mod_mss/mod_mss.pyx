
import logging
import numpy as np

cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def format_val(v):
	import re

	_m = re.match('^\"(.+)\"$', v)
	if _m:
		return _m.group(1)

	if re.match('^\-?\d+$', v):
		return int(v)

	if re.match('^\-?\d*\.\d+$', v):
		return float(v)

	return v

def parse_metadata(f):
	_ms = {}
	with open(f) as _fi:
		for _l in _fi.read().splitlines():
			_vs = _l.split(' = ')
			if len(_vs) == 2:
				_ms[_vs[0].strip()] = format_val(_vs[1].strip())

	return _ms

def load_file(f, fzip):
	'''load files'''
	import os

	_d_tmp = fzip.generate_file()
	os.makedirs(_d_tmp)

	import tarfile
	tarfile.open(f, 'r:gz').extractall(_d_tmp)

	_fs = []
	_fo = []
	_fm = None
	for _root, _dirs, _files in os.walk(_d_tmp):
		for _file in _files:
			_f = os.path.join(_root, _file)

			if _f.endswith('_MTL.txt'):
				_fm = _f
				continue

			if _f.lower().endswith('.tif'):
				_fs.append(_f)
				continue

			_fo.append(_f)

	return sorted(_fs), _fm, _fo

def read_ozone(f):
	from pyhdf.SD import SD, SDC

	_hdf = SD(f, SDC.READ)
	_lyr = _hdf.select(2)

	return _lyr

def read_reanalysis(f):
	from pyhdf.SD import SD, SDC

	_hdf = SD(f, SDC.READ)
	_lyr_slp = _hdf.select(2)
	_lyr_wv = _hdf.select(3)

	return _lyr_slp, _lyr_wv

def read_dem(f):
	from pyhdf.SD import SD, SDC

	_hdf = SD(f, SDC.READ)
	_lyr = _hdf.select(0)

	return _lyr

def read_grid(g, lon, lat):
	_inf = g.info()[2]

	if len(_inf) not in [2, 3]:
		raise Exception('wrong grid size')

	_r = lambda v, x: max(min(v, x), 0)

	_row = _r(int((90.0 - lat) / (180 / _inf[-2])), _inf[-2])
	_col = _r(int((180.0 + lon) / (360.0 / _inf[-1])), _inf[-1])

	if len(_inf) == 2:
		return g[_row, _col]
	else:
		return [g[_b, _row, _col] for _b in xrange(_inf[0])]

def build_band(dat):
	import geo_raster_c as ge

	_dat = np.array(dat[:, :])
	_saz = _dat.shape
	assert(len(_saz) == 2)

	return ge.geo_band_cache(_dat, [-180, 360.0 / _saz[0], 0, 90, 0, -180 / _saz[1]], 
			ge.proj_from_epsg(), nodata=-9999, pixel_type=ge.pixel_type('short'))

def ancdata(t, lon, lat):
	import config
	_d_anc = config.cfg.get('conf', 'anc_path')

	import datetime
	_t = datetime.datetime(max(t.year, 1978), t.month, t.day)

	import os
	_f_ozone = os.path.join(_d_anc, 'EP_TOMS', _t.strftime('ozone_%Y'), _t.strftime('TOMS_%Y%j.hdf'))
	_f_reana = os.path.join(_d_anc, 'REANALYSIS', t.strftime('RE_%Y'), t.strftime('REANALYSIS_%Y%j.hdf'))
	_f_dem = os.path.join(_d_anc, 'CMGDEM.hdf')

	_g = lambda g: read_grid(g, lon, lat)

	_g_ozone = read_ozone(_f_ozone)
	_g_reana = read_reanalysis(_f_reana)
	_g_dem = read_dem(_f_dem)

	_vs = map(_g, [_g_ozone, _g_reana[0], _g_reana[1]])

	return _vs[0] / 1000.0,  \
			map(lambda v: (v*0.0099999998+277.64999)/10.0, _vs[2]), \
			map(lambda v: (v*1.0+119765.)/100.0, _vs[1]), \
			build_band(_g_dem)

def linear_interpolate(vs, t):
	_d = float(24 / len(vs))
	_v = t.split(':')
	_t = float(_v[0]) + float(_v[1]) / 60.0

	_t0 = int(_t / _d)
	_t1 = min(_t0 + 1, len(vs) - 1)

	_td = (vs[_t1] - vs[_t0]) * ((_t - (_d * _t0)) / _d)
	return vs[_t0] + _td

def comptransray(xtaur, xmus):
	import math

	ddiftt=(2.0/3.0+xmus)+(2.0/3.0-xmus)*math.exp(-xtaur/xmus)
	ddiftt=ddiftt/((4.0/3.0)+xtaur)-math.exp(-xtaur/xmus)
	ddirtt=math.exp(-xtaur/xmus)

	return ddirtt+ddiftt

def fintexp1(xtau):
	a = [-0.57721566,0.99999193,-0.24991055, 0.05519968,-0.00976004,0.00107857]
	xx=a[0]
	xftau=1.0

	for i in xrange(5):
		xftau=xftau*xtau
		xx=xx+a[i]*xftau

	import math
	return xx-math.log(xtau)

def fintexp3(xtau):
	import math
	return (math.exp(-1.0 * xtau)*(1.0-xtau)+xtau*xtau*fintexp1(xtau))/2.0

def local_csalbr(xtau):
	import math
	xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*math.exp(-xtau))
	return xalb/(4.+3*xtau)

def local_chand(xphi, xmuv, xmus, xtau):
	# input parameters: xphi,xmus,xmuv,xtau
	# xphi: azimuthal difference between sun and observation (xphi=0,
	# in backscattering) and expressed in degree (0.:360.)
	# xmus: cosine of the sun zenith angle
	# xmuv: cosine of the observation zenith angle
	# xtau: molecular optical depth
	# output parameter: xrray : molecular reflectance (0.:1.)
	# constant : xdep: depolarization factor (0.0279)

	# parameter (fac = 0.017453293)
    # real xdep,pl(10)
    # real fs0,fs1,fs2
    # real as0(10),as1(2),as2(2)
    #                   real xphi,xmus,xmuv,xtau,xrray,pi,phios,xcosf1,xcosf2
    #                   real xcosf3,xbeta2,xfd,xph1,xph2,xph3,xitm, xp1, xp2, xp3
    #                   real cfonc1,cfonc2,cfonc3,xlntau,xitot1,xitot2,xitot3
    #                   integer i
	as0 = [0.33243832,-6.777104e-02,.16285370 \
			 ,1.577425e-03,-0.30924818,-1.240906e-02,-0.10324388 \
			 ,3.241678e-02,0.11493334,-3.503695e-02]
	as1 = [0.19666292, -5.439061e-02]
	as2 = [0.14545937,-2.910845e-02]

	# pi=3.1415927
	# fac=pi/180.
	fac = 0.017453293

	import math

	phios=180.0 - xphi
	xcosf1=1.0
	xcosf2=math.cos(phios*fac)
	xcosf3=math.cos(2*phios*fac)
	xbeta2=0.5
	xdep=0.0279
	xfd=xdep/(2-xdep)
	xfd=(1-xfd)/(1+2*xfd)
	xph1=1+(3*xmus*xmus-1)*(3*xmuv*xmuv-1)*xfd/8.
	xph2=-xmus*xmuv*math.sqrt(1-xmus*xmus)*math.sqrt(1-xmuv*xmuv)
	xph2=xph2*xfd*xbeta2*1.5
	xph3=(1-xmus*xmus)*(1-xmuv*xmuv)
	xph3=xph3*xfd*xbeta2*0.375
	xitm=(1-math.exp(-xtau*(1/xmus+1/xmuv)))*xmus/(4*(xmus+xmuv))
	xp1=xph1*xitm
	xp2=xph2*xitm
	xp3=xph3*xitm
	xitm=(1-math.exp(-xtau/xmus))*(1-math.exp(-xtau/xmuv))
	cfonc1=xph1*xitm
	cfonc2=xph2*xitm
	cfonc3=xph3*xitm
	xlntau=math.log(xtau)

	pl = [0.0] * 10
	pl[0]=1.0
	pl[1]=xlntau
	pl[2]=xmus+xmuv
	pl[3]=xlntau*pl[3]
	pl[4]=xmus*xmuv
	pl[5]=xlntau*pl[5]
	pl[6]=xmus*xmus+xmuv*xmuv
	pl[7]=xlntau*pl[7]
	pl[8]=xmus*xmus*xmuv*xmuv
	pl[9]=xlntau*pl[9]
	fs0=0.0
	for i in xrange(10):
		fs0=fs0+pl[i]*as0[i]
	fs1=pl[0]*as1[0]+pl[1]*as1[1]
	fs2=pl[0]*as2[0]+pl[1]*as2[1]

	xitot1=xp1+cfonc1*fs0*xmus
	xitot2=xp2+cfonc2*fs1*xmus
	xitot3=xp3+cfonc3*fs2*xmus
	xrray=xitot1*xcosf1
	xrray=xrray+xitot2*xcosf2*2
	xrray=xrray+xitot3*xcosf3*2
	xrray=xrray/xmus

	return xrray

def create_output_raster(dat, bnd, f_out, nodata):
	import geo_raster_c as ge

	_cols = bnd.width
	_rows = bnd.height
	_geo = bnd.geo_transform
	_prj = bnd.proj

	ge.geo_raster.create(f_out, [dat.shape[0], _rows, _cols], _geo, _prj, driver='HFA',
			nodata=nodata, pixel_type=ge.pixel_type('short'))

	_img = ge.open(f_out, True)
	for _b in xrange(dat.shape[0]):
		_img.get_band(_b+1).write(dat[_b, :, :], 0, 0)
	_img.flush()

def proc_mss(f, f_out, fzip):
	'''process the MSS'''

	cdef np.ndarray[np.int16_t, ndim=3] _dat_out = None
	cdef np.ndarray[np.uint8_t, ndim=2] _dat
	cdef np.ndarray[np.int16_t, ndim=2] _dat_dem

	cdef int _val, _out, _alti

	_fs, _fm, _fo = load_file(f, fzip)
	if _fm == None:
		raise Exception('failed to find the metadata')

	_bs = [4, 5, 6, 7]

	_ms = parse_metadata(_fm)

	_lmin = [_ms['RADIANCE_MINIMUM_BAND_%d' % _d] for _d in _bs]
	_lmax = [_ms['RADIANCE_MAXIMUM_BAND_%d' % _d] for _d in _bs]

	_qmin = [_ms['QUANTIZE_CAL_MIN_BAND_%d' % _d] for _d in _bs]
	_qmax = [_ms['QUANTIZE_CAL_MAX_BAND_%s' % _d] for _d in _bs]

	# _nr = _ms['REFLECTIVE_LINES']
	# _nc = _ms['REFLECTIVE_SAMPLES']

	_xts = 90.0 - _ms['SUN_ELEVATION']
	# _xfs = _ms['SUN_AZIMUTH']

	_time = _ms['SCENE_CENTER_TIME']

	_lat11 = _ms['CORNER_UL_LAT_PRODUCT']
	_lon11 = _ms['CORNER_UL_LON_PRODUCT']
	_lat12 = _ms['CORNER_UR_LAT_PRODUCT']
	_lon12 = _ms['CORNER_UR_LON_PRODUCT']
	_lat21 = _ms['CORNER_LL_LAT_PRODUCT']
	_lon21 = _ms['CORNER_LL_LON_PRODUCT']
	_lat22 = _ms['CORNER_LR_LAT_PRODUCT']
	_lon22 = _ms['CORNER_LR_LON_PRODUCT']

	import datetime
	import math
	import os

	_date = datetime.datetime.strptime(_ms['DATE_ACQUIRED'], '%Y-%m-%d')
	_f_out = f_out
	if os.path.isdir(_f_out):
		_f_out = os.path.join(_f_out, 'sr_l%s_p%sr%s_%s.img' % ( \
			_ms['SPACECRAFT_ID'].split('_')[1], _ms['WRS_PATH'], \
			_ms['WRS_ROW'], _date.strftime('%Y%m%d')))

	_cpi = math.atan(1.0) * 4.0
	# _pi = _cpi

	_om=(.9856*float(int(_date.strftime('%j'))-4)) * math.pi / 180.0
	_dsol=1.0 / ((1.0 - 0.01673 * math.cos(_om)) ** 2.0)

	_mlat=(_lat11+_lat12+_lat21+_lat22)/4.0
	_mlon=(_lon11+_lon12+_lon21+_lon22)/4.0

	# compute mean ozone, water vapor and pressure (slp) for the whole scene
	_v_ozone, _v_wv, _v_slp, _bnd_dem = ancdata(_date, _mlon, _mlat)

	# _wvint = linear_interpolate(_v_wv, _time)
	_slpint = linear_interpolate(_v_slp, _time)

	_esun = [1805.42,1511.42,1241.4,908.018]
	_ah2o = [0.00251053,0.00768567,0.0249089,0.105129]
	_bh2o = [0.753595,0.687902,0.609682,0.426469]
	_aoz = [-0.084081,-0.0582646,-0.00965367,0.0]
	_a1 = [9.82973e-06,0.00716206,0.0297642,8.56296e-06]
	_b0 = [0.639237,0.614214,0.690885,0.140782]
	_b1 = [-0.240851,0.193631,0.564946,-0.110042]
	_tauray = [0.09817,0.04713,0.02862,0.01428]

	_irad = [_esun[i] * math.cos(_xts * _cpi / 180.0) * _dsol / _cpi for i in xrange(4)]

	_gain = [(_lmax[i] - _lmin[i]) / (_qmax[i] - _qmin[i]) for i in xrange(4)]
	_offset = [_lmin[i] for i in xrange(4)]

	_xtv = 0.0
	_m=1.0/ math.cos(_xts*_cpi/180.0)+1/math.cos(_xtv*_cpi/180.0)

	_uh2o=3.0
	# _uoz = 0.35
	_uoz = _v_ozone
	_pres = _slpint / 1013.0
	_xphi = 0.0

	import geo_raster_c as ge
	import geo_raster_ex_c as gx

	_bnd_tmp = None
	_grid = None

	for _ib in xrange(len(_fs)):
		logging.info('process band %s' % _ib)
		_tgwv = math.exp(-1 * _ah2o[_ib] * ((_m * _uh2o) ** _bh2o[_ib]))
		_tgoz = math.exp(_aoz[_ib] * _m * _uoz)
		_tgog=math.exp(-1 * (_a1[_ib]*_pres)*(_m**(math.exp(-1 * \
				(_b0[_ib]+_b1[_ib]*_pres)))))
		_xtaur=_tauray[_ib]*_pres
		_xmus=math.cos(_xts*_cpi/180.0)
		_xmuv=math.cos(_xtv*_cpi/180.0)

		_xtts = comptransray(_xtaur, _xmus)
		_xttv = comptransray(_xtaur, _xmuv)
		# Compute total transmission (product downward by  upward)
		_ttatm = _xtts * _xttv
		_satm = local_csalbr(_xtaur)
		_roray = local_chand(_xphi, _xmus, _xmuv, _xtaur)

		_bnd = ge.open(_fs[_ib]).get_band()
		if _bnd_tmp == None:
			_bnd_tmp = _bnd
			_grid = gx.projection_transform.from_band(_bnd_tmp, \
					ge.proj_from_epsg())
			_dat_out = np.zeros((len(_fs), _bnd_tmp.height, _bnd_tmp.width), dtype=np.int16)

		_bnd = _bnd.read_block(_bnd_tmp)

		# import progress_percentage
		# _ppp = progress_percentage.progress_percentage(_bnd.height)

		_dat = _bnd.data
		_dat_dem = _bnd_dem.read_block(_bnd_tmp).data

		for _row in xrange(_bnd.height):
			# _ppp.next()
			for _col in xrange(_bnd.width):
				# if not (_col == 2139 - 1 and _row == 1004-1):
				# 	continue

				_val = _dat[_row, _col]
				if _val == 0:
					_dat_out[_ib, _row, _col] = -9999
					continue

				_alti = _dat_dem[_row, _col]
				_alti = max(_alti, 0)

				_pres = _slpint * math.exp(-1 * _alti / 8500.0) / 1013.0
				_xtts = comptransray(_xtaur, _xmus)
				_xttv = comptransray(_xtaur, _xmuv)

				_ttatm = _xtts * _xttv
				_satm = local_csalbr(_xtaur)
				_roray = local_chand(_xphi, _xmus, _xmuv, _xtaur)
				_rotoa = (_val * _gain[_ib] + _offset[_ib]) / _irad[_ib]
				_xx = _rotoa / (_tgog * _tgoz)
				_xx = (_xx - _roray) / (_ttatm * _tgwv)
				_xx = _xx / (1.0 + _satm * _xx)

				_dat_out[_ib, _row, _col] = int(_xx * 10000.0)
				# if _ib == 0 and _row == 1004 - 1 and _col == 2139 - 1:
				# 	print _dat_out[_ib, _row, _col], _tgog, _tgoz, _roray, _ttatm, _tgwv

		# _ppp.done()

	create_output_raster(_dat_out, _bnd, _f_out, -9999)

	# also keep the metadata
	with open(_f_out[:-3] + 'met', 'w') as _fo:
		_fo.write(open(_fm).read())

