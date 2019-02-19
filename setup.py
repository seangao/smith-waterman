#!/usr/bin/ python

from distutils.core import setup

setup(name='Smith-Waterman',
      version='1.0',
      description=' Smith-Waterman local alignment algorithm for protein sequences',
      author='Sean Gao',
      author_email='sean.gao@yale.edu',
      packages=['smith_waterman'],
      entry_points={
        'console_scripts': [
            'smith_waterman = smith_waterman.__main__:main'
        ]
      }
     )
