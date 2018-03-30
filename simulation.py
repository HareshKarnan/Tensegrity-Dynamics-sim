from visual import *
from csv import *
from visual.controls import *
import matlab.engine
import numpy as np
import time
eng = matlab.engine.start_matlab()
#ret = eng.dynamics_sim()
#Replace the ___ in eng.___ with the function cum file name of the matlab file that runs the dynamics simulation
# and returns the simulation data over time.

#ret = eng.class_k_prism_reducedorder() # object invoking the matlab file
ret = eng.class_k_prism_control()
numofframes = len(ret["S"][1][1])
nb = len(ret["B"][1]) # no of bars
ns = len(ret["S"][1]) # no of strings
nf = len(ret["B"][1][1][:])  # number of frames in the dynamics simulation
rod_b = list() # empty list of bar objects in simulation
rod_s = list() # empty list of string objects in simulation
nd = list() # empty list for node black balls in simulation


def simulate():
    #ground = box(pos=(0,-0.1,0),length=100,height=0.1,width=100,opacity=1,color=color.yellow,material=materials.wood)
    #lamp = distant_light(direction=(0,20,0), color=color.yellow,intensity = 0.1)
    for dt in range(0,nf):
        if dt>1:
            for objs in rod_b:
                objs.visible = False
                del objs
            for objs in rod_s:
                objs.visible = False
                del objs
            for objs in nd:
                objs.visible = False
                del objs
        for i in range(0,nb):
            # drawing bars
            b_com = vector(ret["RB"][0][i][dt],ret["RB"][1][i][dt],ret["RB"][2][i][dt])
            b_vect = vector(ret["B"][0][i][dt],ret["B"][1][i][dt],ret["B"][2][i][dt])
            rod_b.append(cylinder(pos=b_com-b_vect/2,axis=b_vect, radius=0.08,color=(42/255, 54/255, 59/255),opacity=1,material = materials.wood,scene=scene))
        for i in range(0,ns):
            # drawing strings
            s_com = vector(ret["RS"][0][i][dt],ret["RS"][1][i][dt],ret["RS"][2][i][dt])
            s_vect = vector(ret["S"][0][i][dt], ret["S"][1][i][dt], ret["S"][2][i][dt])
            rod_s.append(cylinder(pos=s_com - s_vect / 2, axis=s_vect, radius=0.01,color=color.red,material=materials.emissive,scene=scene))
        for i in range(0,nb*2):
            # drawing nodes
            n = vector(ret["N"][0][i][dt],ret["N"][1][i][dt],ret["N"][2][i][dt])
            nd.append(sphere(pos=n,radius=0.025,color=color.black,scene=scene))
        if dt<nf-1:
            # regulate the rate of animation
            rate(1/(ret["dt"][0][dt+1]-ret["dt"][0][dt]))
            #print ret["dt"][0][0]
            #rate(1/0.5)
            pass
    for objs in rod_b:
        objs.visible = False
        del objs
    for objs in rod_s:
        objs.visible = False
        del objs
    for objs in nd:
        objs.visible = False
        del objs
    b.text = 'Restart Again'
w = window(menus=True, title='Tensegrity World', width=1020, height=728)
green = (0.184,0.5764,0.5843)
scene = display(window=w, width=1020, height=728,background=(0.118,0.565,1))
#text(window=w,text='My text is\ngreen',align='center', depth=-0.3, color=color.green)
b = button(scene = scene,pos=(0,0), width=100, height=100,text='Restart', action=lambda: simulate() )
simulate()
