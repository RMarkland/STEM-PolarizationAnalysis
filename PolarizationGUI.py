#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 16:35:11 2023

@author: reidmarkland
"""

import tkinter as tk
import hyperspy as hs
import cv2
from PIL import Image, ImageTk
# from tkinter import *
from tkinter import ttk, Tk, Canvas, Frame, Button, LEFT, RIGHT, BOTTOM, TOP, filedialog
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from PolarizationFunctionality import *
from PolarizationMap import *

class ImageCropper():
    
    def __init__(self, root):
        self.root = root
        self.root.title("Lattice Displacement Map")
        self.img = None
        self.crop_area = None
        self.line = None
        self.image_select = None
    
        self.tabControlMain = ttk.Notebook(self.root)
        
        self.tab1 = ttk.Frame(self.tabControlMain, width=600, height=600)
        self.tabControlMain.add(self.tab1, text='Cropper')
        self.tabControlMain.pack(expand=1, fill='both')
        
        self.canvas = Canvas(self.tab1, width=600, height=600)
        self.canvas.pack(side=LEFT)
        self.canvas.bind("<Button-1>", self.on_mouse_down)
        self.canvas.bind("<Double-1>", self.on_double_click)
        
        control_frame1 = Frame(self.tab1)
        control_frame1.pack(side=RIGHT, padx=10)
        
        open_button = Button(control_frame1, text='Open Image', command=self.open_image)
        open_button.pack(pady=10)
        
        crop_button = Button(control_frame1, text='Crop Image', command=self.crop_image)
        crop_button.pack(pady=10)
        
        lattice_button = Button(control_frame1, text = 'Generate Lattice', command = self.generate_lattice)
        lattice_button.pack(pady=10)
        
        
    def open_image(self):
        file_path = filedialog.askopenfilename(title = 'Select file',\
                                               filetypes=[('Image Files', '*.png')])
        if not file_path:
            return
        self.canvas.delete('all')
        self.img = cv2.imread(file_path, 0)
        # self.img = cv2.cvtColor(self.img, cv2.COLOR_BGR2RGB)
        height, width = self.img.shape
        self.scale_factor = min(600/width, 600/height)
        resized_img = cv2.resize(self.img, (int(width*self.scale_factor), int(height*self.scale_factor)))
        self.img_tk = ImageTk.PhotoImage(image=Image.fromarray(cv2.flip(resized_img,1)))
        self.canvas.config(width=self.img_tk.width()+10, height=self.img_tk.height()+10)
        self.canvas.create_image(5, 5, anchor='nw', image=self.img_tk)
        
        self.lineY = [5, height]
        self.crop_count = 0
        self.map_count = 0
        self.cropDict = {'1': self.img}

        
    def on_mouse_down(self, event):
        if self.img is not None:
            if self.line is not None:
                self.canvas.delete(self.line)
            self.line = self.canvas.create_line(0, event.y, self.img_tk.width(), event.y, fill='red')
            self.crop_area = (0, int(event.y/self.scale_factor), int(self.img_tk.width()/self.scale_factor),\
                              int(self.img_tk.height()/self.scale_factor))
            
                
    def crop_image(self):
        if self.img is not None and self.crop_area is not None:
            self.cropDict = {}; self.crop_tk = {}
            self.crop_count += 1
            if self.crop_count < 2:
                self.lineY.append(self.crop_area[1]); self.lineY.sort()
            else:
                self.lineY.append(self.crop_area[1]-int(10/self.scale_factor)); self.lineY.sort()
                            
            self.canvas.delete('all')
            self.canvas.config(width=self.img_tk.width()+10,height=self.img_tk.height()+self.crop_count*10+10)
            
            for i in range(self.crop_count+1):
                self.cropDict[str(i)] = self.img[self.lineY[i]:self.lineY[i+1], self.crop_area[0]:]
                height, width = self.cropDict[str(i)].shape
                resised_crop = cv2.resize(self.cropDict[str(i)], (int(width*self.scale_factor),\
                                                        int(height*self.scale_factor)))
                self.crop_tk[str(i)] = ImageTk.PhotoImage(image=Image.fromarray(cv2.flip(resised_crop,1)))
                self.canvas.create_image(5, int(self.lineY[i]*self.scale_factor)+10*i+5,\
                                         anchor='nw', image=self.crop_tk[str(i)])
            
            
    def on_double_click(self, event):
        if self.crop_count is not None:
            if self.image_select is not None:
                self.canvas.delete(self.image_select)
            self.canvas.delete(self.line)
            self.select = event.y
            selected_image, self.selected_crop = determine_crop(self)
            
            if self.crop_count==0:
                self.cropDict['0'] = self.img
            self.image_select = self.canvas.create_rectangle(5,(selected_image[0]+5)*self.scale_factor,\
                                    self.img_tk.width()+5,(selected_image[1]+5)*self.scale_factor,\
                                        outline='red',width=5)
                
                
    def generate_lattice(self):
        if self.selected_crop is not None:    
            self.map_count += 1            
            self.tab2 = ttk.Frame(self.tabControlMain, width = 600, height = 600)
            self.tabControlMain.add(self.tab2, text='Site Map'+str(self.map_count))
            self.tabControlMain.pack(expand=True, fill='both')
                                                
            pol_img = hs.signals.Signal2D(self.cropDict[str(self.selected_crop)])
            sublattice_A, image_noA = GetLatticeA(pol_img, separation=20) # Add separation tool
            sublattice_B = GetLatticeB(image_noA, separation=24)
            self.atom_lattice = am.Atom_Lattice(image = pol_img, name = 'test', sublattice_list=\
                                            [sublattice_A, sublattice_B])
                
            self.control_frame2 = Frame(self.tab2)
            self.control_frame2.pack(side=BOTTOM, anchor='sw', pady=5, padx=5, fill='x')
            
            map_button = Button(self.control_frame2, text='Generate Map', command=self.generate_map)
            map_button.pack(pady=10)
            
            self.atom_fig = plot_atomap(self.atom_lattice)
            self.site_canvas = FigureCanvasTkAgg(self.atom_fig, master=self.tab2)
            self.site_canvas.draw()
            self.site_canvas.get_tk_widget().pack(side=LEFT, expand=True, fill='y')
                                                
            
            self.root.update()
            self.tabControl.select(self.tab2)
            
            
    def generate_map(self):
        if self.atom_lattice is not None:
            pol_img = hs.signals.Signal2D(self.cropDict[str(self.selected_crop)])
            self.fig, ax = PolarizationMap(self.atom_lattice, pol_img)
            
            self.tab2.config(width=850,height=400)
            self.control_frame2.pack(side=BOTTOM, pady=5, padx=5, fill='x')
            self.site_canvas.get_tk_widget().pack(side=LEFT,expand=True, fill='y')
            
            if self.pol_canvas is not None:
                self.pol_canvas.delete('all')
            self.pol_canvas = FigureCanvasTkAgg(self.fig, master=self.tab2)
            self.pol_canvas.draw()
            self.pol_canvas.get_tk_widget().pack(side=LEFT, expand=True,fill='y')
            
            self.root.update()
            # self.site_canvas.get_tk_widget().config(width=300, height=300)
            
            
            
            
            
            
    def close_root(self):
        if self.img is not None:
            self.root.destroy()
            self.root.quit()
    
        
if __name__ == '__main__':
    root = Tk()
    app = ImageCropper(root)
    root.mainloop()



# TASK LIST
#
# GET POL MAP WORKING FOR ONE TAB XX
#
# CREATE TAB LIST. POSSIBLE SEPARATE TABS FOR EACH CROPS, EACH WITH ITS OWN POL MAP
#
# CREATE SEPARATION MODIFICATION
#
# CREATE SELECTION/DELETION OF ATOMIC POSITIONS IN TAB
#
# POL MAP MODS
#
# 'UNCROP' BUTTON (DELETES CONTRADICTING TABS)
#
# FINALIZE BUTTON: STICHES POL MAPS TOGETHER (MAYBE STITCHING ATOMAPS TOGETHER, 
#                   THEN RE-CALCULATE POL MAP)
#
# LOADING NOTIFIER
#
# BUG PATCHING BUG PATCHING BUG PATCHING
#
#
        
        
        