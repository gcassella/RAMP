#%%
import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import simps

import os, csv, re

from time import time

import tkinter as tk
from PIL import Image, ImageTk

from enum import Enum
from random import random

os.environ['PYOPENCL_CTX']='0:0'
os.environ['PYOPENCL_COMPILER_OUTPUT']='1'

class MonitorVars(Enum):
    THETA = 0
    TOF = 1

#%%

def LoadSQW(fn):
    with open(fn, 'r') as fin:
        lines = fin.readlines()
        q = np.fromstring(lines[17], dtype=np.float32, sep=' ')
        w = np.fromstring(lines[20], dtype=np.float32, sep=' ')
        rho = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[6])[0])
        sigma_abs = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[9])[0])
        sigma_coh = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[10])[0])

    sqw = np.loadtxt(fn, skiprows=23).astype(np.float32) + np.finfo(float).eps.astype(np.float32)
    sigma_scat = sigma_coh

    pw = simps(np.multiply(sqw.T, q), axis=1, x=q) / np.linalg.norm(sqw)
    pw_cdf = np.array([simps(pw[:i], x=w[:i]) for i in range(1,len(w)+1)], dtype=np.float32)
    pw_cdf = pw_cdf / np.max(pw_cdf)
    pw_cdf[np.argmax(pw_cdf):] = 1

    pq = np.array([sqw[:,i] / np.linalg.norm(sqw[:,i]) for i in range(0,len(w))])
    pq_cdf = np.array(
        [[simps(pq[j,:i], x=q[:i]) / simps(pq[j,:len(q)+1], x=q[:len(q)+1]) \
        for i in range(1,len(q)+1)] for j in range(0,len(w))], dtype=np.float32)

    return (q, w, rho, sigma_abs, sigma_scat, pw_cdf, pq_cdf, sqw)

q, w, rho, sigma_abs, sigma_scat, pw_cdf, pq_cdf, sqw = LoadSQW('He4_liq_coh.sqw')

#%%
plt.plot(q,pq_cdf[0])
plt.show()

#%%
from random import random

v = random()

N = 100000

pq_cdf = pq_cdf.flatten()
omega = np.zeros((N,))
Q = np.zeros((N,))

for n in range(N):
    u = random()
    i = 0
    mindiff = 1.
    
    windex = 0

    for pw in pw_cdf:
        if abs(pw - u) < mindiff:
            mindiff = abs(pw - u)
            omega[n] = w[i]
            windex = i

        i+=1

    v = random()
    i=0
    mindiff = 1.

    for pq in pq_cdf[windex*len(q):windex*len(q)+len(q)]:
        if abs(pq - v) < mindiff:
            mindiff = abs(pq - v)
            Q[n] = q[i]

        i+=1

H, xe, ye = np.histogram2d(Q, omega, bins=100)
X, Y = np.meshgrid(xe, ye)
H = H.T
plt.pcolormesh(X, Y, H)
plt.show()

#%%

H,_ = np.histogram(omega, bins=len(w))
plt.plot(w, H)
plt.show()

#%%

H,_ = np.histogram(Q, bins=len(q))
plt.plot(q, H)
plt.show()

#%%

#pq_cdf = pq_cdf.flatten()

N = 100000
Q = np.zeros((N,))

for n in range(N):
    v = random()
    i=0
    mindiff = 1.
    for pq in pq_cdf[110*len(q):110*len(q)+len(q)]:
        if abs(pq - v) < mindiff:
            mindiff = abs(pq - v)
            Q[n] = q[i]
        i+=1

H,_ = np.histogram(Q, bins=len(q))
plt.plot(q, H)
plt.show()

#%%

sqw = np.loadtxt('diffuse.sqw', skiprows=22).astype(np.float32) + np.finfo(float).eps.astype(np.float32)
pq = np.array([sqw[:,i] / np.linalg.norm(sqw[:,i]) for i in range(0,len(w))])

plt.plot(q,pq[0])
plt.show()