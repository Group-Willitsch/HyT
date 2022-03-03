#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 11:16:19 2019

@author: tomislavdamjanovic
"""

from numba import jit


def fun1(x,n):
   B=0.
   for i in range(n):
      B=B+x

@jit
def fun2(x,n):
   B=0.
   for i in range(n):
      B=B+x

import time

start1 = time.time()
fun1(0.5,10000000)
end1 = time.time()
print(end1-start1)

start1 = time.time()
fun2(0.5,10000000)
end1 = time.time()
print(end1-start1)
