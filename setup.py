from setuptools import setup

import kraken_trawl

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='kraken_trawl',
      version=kraken_trawl.__version__,
      description=kraken_trawl.__description__,
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GPLv3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics kraken',
      url=kraken_trawl.__url__,
      author=kraken_trawl.__author__,
      author_email=kraken_trawl.__author_email__,
      license=kraken_trawl.__license__,
      packages=['kraken_trawl'],
      install_requires=[
          'click',
          'pandas',
      ],
      test_suite='nose.collector',
      tests_require=[],
      entry_points={
          'console_scripts': ['kraken-trawl=kraken_trawl.kraken_trawl:kraken_trawl'],
      },
      include_package_data=True,
      zip_safe=False)
