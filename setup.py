# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 16:19:07 2015

@author: carolc24
"""

from setuptools import setup

setup(name='growthsim',
      version='1.0',
      description='A tool for simulating populations of engineered cells',
      author='Caroline Cannistra',
      author_email='carolc24@uw.edu',
      packages=['growthsim'],
      zip_safe=False)

f = open('growthsim/__init__.py', 'r');
init = f.read();
f.close();
g = open('__init__.py', 'w');
g.write(init);
g.close();
