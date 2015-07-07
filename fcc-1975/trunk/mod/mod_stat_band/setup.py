from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("mod_stat_band", ["mod_stat_band.pyx"],
		extra_compile_args=["-O3", "-march=native", "-ffast-math","-funroll-loops"],
		define_macros=[("NPY_NO_DEPRECATED_API", None)])]

setup(
	 name = "",
	 cmdclass = {"build_ext": build_ext},
	 ext_modules = ext_modules
)
