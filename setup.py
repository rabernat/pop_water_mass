from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='watermasstools',
      version='0.1',
      description='Tools for analyzing ocean satellite data',
      url='https://bitbucket.org/ryanaberanthey/pop_water_mass',
      author='Ryan Abernathey',
      author_email='rpa@ldeo.columbia.edu',
      license='MIT',
      packages=['watermasstools'],
      install_requires=[
          'numpy','scipy','netCDF4'
      ],
      test_suite = 'nose.collector',
      zip_safe=False)
