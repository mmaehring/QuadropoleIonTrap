# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 14:10:37 2021

@author: marcu
"""

import skimage as ski
import skimage.draw as drw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image

path = r"C:\Users\marcu\OneDrive\Desktop\PraktikumIII\QuadropoleIonTrap\plots\FEM_Correction_alpha.png"

img = image.imread(path)

plt.figure(dpi=340)
plt.imshow(img)

# img[790:801, 1000:1151] = [255, 0, 0, 1]

# plt.imshow(img)

zero_point_graph = np.array([84, 699]) # (x, y)

# 1 pixel in d/r
conv_fac_x = 101/2.5 # 101 pixels from 0 to 2.5
conv_fac_y = 116/0.5 # 116 pixels from 0 to 0.5


# dr_val = np.float(input("Input the d/r value: "))
# dr_val = 1.33
dr_val = 7.5

orange_val = np.array([1., 0.49803922, 0.05490196, 1.])
blue_val = np.array([0.12156863, 0.46666667, 0.7058824 , 1.])
black_val = np.array([0., 0., 0., 1.])


start_point = zero_point_graph + np.array([dr_val * conv_fac_x, 0]) 
start_point = int(start_point)

start_point[0], start_point[1] = start_point[1], start_point[0]

temp_point = start_point

rr, cc = drw.disk(start_point, 5)
img[rr, cc] = np.array([0, 255, 0, 1])

plt.figure(dpi=250)
plt.imshow(img)

crossed_orange = False

# while ((img[temp_point, 0:4][1,1] != blue_val).any()) and (abs(temp_point) < img.shape[0:2]).any():
#     if (img[temp_point, 0:4][1,1] == orange_val).all():
#         crossed_orange = True
#     temp_point -= np.array([1, 0])
#     rr, cc = drw.disk(temp_point, 2)
#     img[rr, cc] = np.array([255, 0, 0, 1])
#     print(temp_point)

# print(temp_point)    
    
