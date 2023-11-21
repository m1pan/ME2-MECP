#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 18:27:18 2022

@author: ncinosi
"""

# Gauss-Legendre quadrature nodes and weights
tg = np.array([ [0],
     [1/np.sqrt(3), -1/np.sqrt(3)], \
     [0, np.sqrt(3/5), -np.sqrt(3/5)], \
     [np.sqrt(3/7-2/7*np.sqrt(6/5)), -np.sqrt(3/7-2/7*np.sqrt(6/5)), \
            np.sqrt(3/7+2/7*np.sqrt(6/5)), -np.sqrt(3/7+2/7*np.sqrt(6/5))], \
     [0, 1/3*np.sqrt(5-2*np.sqrt(10/7)), -1/3*np.sqrt(5-2*np.sqrt(10/7)), \
     1/3*np.sqrt(5+2*np.sqrt(10/7)), -1/3*np.sqrt(5+2*np.sqrt(10/7))] ] )

wg = np.array([ [2], \
     [1, 1], \
     [8/9, 5/9, 5/9], \
     [(18+np.sqrt(30))/36, (18+np.sqrt(30))/36, (18-np.sqrt(30))/36, (18-np.sqrt(30))/36], \
     [128/225, (322+13*np.sqrt(70))/900, (322+13*np.sqrt(70))/900, \
      (322-13*np.sqrt(70))/900, (322-13*np.sqrt(70))/900] ])

