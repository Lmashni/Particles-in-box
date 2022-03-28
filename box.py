#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 19:13:15 2019

@author: lyth
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.animation as anim


N =3
itr = 1000
dt =0.001
movie = False
traj = True
anim = False
#s = np.array([[1,1],[1,10],[1,1]],dtype ='float')
s =np.random.random((N,3,2))
#s[:,1,:] = (s[:,1,:]- np.ones(s[:,1,:].shape)* (1/np.sqrt(2)))*0.1
s[0,0,0],s[0,0,1],s[1,0,0],s[1,0,1],s[2,0,0],s[2,0,1] = 0.6 , 0.5,0.7 , 0.5,0.8 , 0.5
s[0,1,0],s[0,1,1],s[1,1,0],s[1,1,1],s[2,1,0],s[2,1,1] = 0   , 2,0   , 2,0   , 2
s[:,2,:] = np.zeros(s[:,2,:].shape)
#s[:,0,:] = np.zeros(s[:,2,:].shape)


def stp(x,dt,col,wall,gravc,gravz=0):
    c = 0.5 * np.ones(x[:,0,:].shape)
    xl = x[:,0,:]
    x[:,0,:] = x[:,0,:] + x[:,1,:] * dt
    for i in range(len(s)):

    ### walls  ###

        if wall == True:
            if x[i,0,0] > 1 or x[i,0,0] < 0:
                x[i,1,0] = -x[i,1,0]
            if x[i,0,1] > 1 or x[i,0,1] < 0:
                x[i,1,1] = -x[i,1,1]

    ### collisions  ###

        if col == True:
            for j in range(i+1,len(s)):
                dl = np.sqrt(np.sum((xl[i,:]-xl[j,:])**2))
                d = np.sqrt(np.sum((x[i,0,:]-x[j,0,:])**2))
                if (d < 0.05 and dl > 0.025):
                    print('yyyyyyyyyyyy')
                    a,b = x[i,1,0] ,x[i,1,1]
                    x[i,1,0] ,x[i,1,1]= x[j,1,0],x[j,1,1]
                    x[j,1,0] ,x[j,1,1]= a,b
    x[:,1,:] = x[:,1,:] + x[:,2,:] * dt

    ####     cenntral gravity     ###

    if gravc == True:
        d3 = np.sum((c - x[:,0,:]) ** 2 + 0.0001,1)**(3/2)
        D3 = np.ones(s[:,0,:].shape)
        D3[:,0] , D3[:,1] = d3,d3
        x[:,2,:] = (c-s[:,0,:])/D3
    if gravz != 0:
        x[:,2,1] = -gravz
        
    '''for i in range(len(x)):
        for j in range(i+1,len(x)):
            #print(x[i,0,:] - x[j,0,:])
            d3 = (np.sum((x[i,0,:] - x[j,0,:])**2)+0.0001)**(3/2)
            #print(d3)
            x[i,2,:] += (x[i,0,:] - x[j,0,:]) / d3
            x[j,2,:] += -x[i,2,:]'''
    return x


###  put your physics here ###
    
wall = False
col = False
gravc = True
gravz = 0

###     produce time steps  ###

S =np.zeros((itr,N,3,2))
for i in range(itr):
    print(i)
    S[i,...] = stp(s,dt,col,wall,gravc,gravz)

###     make slides for movie   ###
    
if movie == True:
    import matplotlib
    matplotlib.use('Agg')
    for i in range(itr):
        print(i)
        if (i %100 == 0):
            fig = plt.figure()
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.plot(0.5,0.5,'.')
            for j in range(N):
                plt.plot(S[i,j,0,0],S[i,j,0,1],'.')
            #plt.plot(S[i,0,0,0],S[i,0,0,1],'.',S[i,1,0,0],S[i,1,0,1],'.',S[i,2,0,0],S[i,2,0,1],'.')#,s[3,0,0],s[3,0,1],'.',s[4,0,0],s[4,0,1],'.',s[5,0,0],s[5,0,1],'.')  
            plt.savefig('./plot/plt'+str(i).zfill(len(str(itr))))
            plt.clf


###     make slides for movie   ###
    
if anim == True:
    fig, ax = plt.subplots()
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ims =[]
    for i in range(itr):
        print(i)
        if (i %100 == 0):
            #fig = plt.figure()
            
            plt.plot(0.5,0.5,'.')
            for j in range(N):
                im, = ax.plot(S[i,j,0,0],S[i,j,0,1],'.r')
                ims.append([im])
            #plt.plot(S[i,0,0,0],S[i,0,0,1],'.',S[i,1,0,0],S[i,1,0,1],'.',S[i,2,0,0],S[i,2,0,1],'.')#,s[3,0,0],s[3,0,1],'.',s[4,0,0],s[4,0,1],'.',s[5,0,0],s[5,0,1],'.')  
            #plt.savefig('./plot/plt'+str(i).zfill(len(str(itr))))
            #plt.clf
    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)
    plt.show()
###    graph trojactories  ###

if traj == True:
    for i in range(N):
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.plot(0.5,0.5,'.')
        #plt.plot(np.mean(S[:,i,0,0]),np.mean(S[:,i,0,1]),'.')
        plt.plot(S[:,i,0,0],S[:,i,0,1],'-')
        
        
         
