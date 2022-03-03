# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 13:32:14 2019

@author: Gruppe Willitsch
"""

import visa

rm = visa.ResourceManager()
print(rm.list_resources())

tm = rm.open_resource("COM24")

a = tm.query('KRDG? 1')
#
#print(a)

