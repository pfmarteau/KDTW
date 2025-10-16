from distutils.core import setup, Extension
import numpy
# define the extension module
KDTW = Extension('KDTW', sources=['KDTW.c'])

# run the setup
#setup(ext_modules=[KDTW])
setup(
    name='kdtw',
    version='1.0',
    description='Positive Definite Kernel close to the Dynamic Time Warping measure',
    author='Pierre-François Marteau',
    author_email='pierrefrancois.marteau@gmail.com',
    license='MIT License',
    packages=['kdtw'],
    install_requires=['numpy'],
    ext_modules=[
        Extension("KDTW", ["KDTW.c"],
                  include_dirs=[numpy.get_include()]),
    ],
)

