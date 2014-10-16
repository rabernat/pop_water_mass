from setuptools import setup

from numpy.distutils.core import Extension

ext1 = Extension(name = 'jmd95',
                 sources = ['watermasstools/jmd95.F90'])

def readme():
    with open('README.md') as f:
        return f.read()

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='watermasstools',
      version='0.1',
      description='Tools for calculating water mass transformation from POP model',
      url='https://bitbucket.org/ryanaberanthey/pop_water_mass',
      author='Ryan Abernathey',
      author_email='rpa@ldeo.columbia.edu',
      license='MIT',
      packages=['watermasstools'],
      ext_modules=[ext1],
      install_requires=[
          'numpy','scipy','netCDF4'
      ],
      test_suite = 'nose.collector',
      zip_safe=False)
