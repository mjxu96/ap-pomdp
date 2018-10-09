#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 15:17:18 2018

@author: minjun
"""
# In[]
import numpy as np
import matplotlib.pyplot as plt


# In[]
class avPOMDP:
    
    def __init__(self, scan_parameter, gaussian_sigma=5, scan_resolution=1):
        # scan_resolution should be fixed to 1 in the current version
        # the function of changing resolution "might" be added 
        # in the next version :)
        (scan_range, max_range) = scan_parameter
        self.resolution = scan_resolution
        self.scan_range = scan_range
        self.max_range = max_range
        self.x = np.arange(max_range).reshape((max_range, 1))
        self.y_pred = np.ones((max_range, 1)) / max_range
        self.y_gt = np.zeros((max_range, 1))
        self.y_scan = np.ones((max_range, 1)) / max_range
        self.points = []
        self.now_time = -1
        self.scan_angle = []
        self.gaussian_sigma = gaussian_sigma
        
    def update(self, input_data, scan_angle):
        self.now_time += 1
        self.scan_angle.append(scan_angle)
        if (input_data != None):
            (alpha, distance) = input_data
            self.addPoint(alpha, distance, self.now_time)
        self.updatePoints()
        self.updatePredictionDistribution()
        self.updateScanAngelDistribution()
        
    def draw(self):
        self.drawPredictionDistribution()
        self.drawGroundTrueth()
        self.drawScanAngelDistribution()
        
    def updatePoints(self):
        for point in self.points:
            point[2] += 1
        right_scan_bound = self.scan_angle[-1] + self.scan_range/2
        left_scan_bound = self.scan_angle[-1] - self.scan_range/2
        for point in self.points:
            if ((right_scan_bound >= self.max_range) or 
                (left_scan_bound < 0)):
                left_scan_bound = left_scan_bound % self.max_range
                right_scan_bound = right_scan_bound % self.max_range
                if (point[0] >= left_scan_bound or 
                    point[0] <= right_scan_bound):
                    point[2] = self.now_time
                    point[3] = 1
            else:
                if (point[0] >= left_scan_bound and 
                    point[0] <= right_scan_bound):
                    point[2] = self.now_time
                    point[3] = 1
        return
        
    def updatePredictionDistribution(self, sigma_expand_cof=1.1):
        self.y_pred = np.ones((self.max_range, 1)) / self.max_range
        for point in self.points:
            if (point[3] == 1):
                sigma = self.gaussian_sigma * np.power(sigma_expand_cof, 
                                                       (self.now_time-
                                                        point[2]))
                mu = point[0]
                for x in range(mu-3*int(np.round(sigma)), 
                               mu+3*int(np.round(sigma))):
                    x = x % self.max_range
                    coff = 1/point[1]
                    self.y_pred[x] += self.calculateGaussian(x, mu, 
                               sigma, coff)
        self.y_pred = self.y_pred / np.sum(self.y_pred)
        return
    
    def updateScanAngelDistribution(self):
        self.y_scan = np.ones((self.max_range, 1)) / self.max_range
        for x in range(self.max_range):
            left_bnd = int(x - self.scan_range/2)
            right_bnd = int(x + self.scan_range/2)
            if (left_bnd < 0 or right_bnd >= self.max_range):
                left_bnd = left_bnd % self.max_range
                right_bnd = right_bnd % self.max_range
                self.y_scan[x] += (np.sum(self.y_pred[0: right_bnd]) + \
                           np.sum(self.y_pred[left_bnd: self.max_range]))
            else:
                self.y_scan[x] += np.sum(self.y_pred[left_bnd: right_bnd])
        self.y_scan = self.y_scan / np.sum(self.y_scan)
        return
    
    def addPoint(self, alpha, distance, points_added_time):
        self.points.append([alpha, distance, points_added_time, 0])
        points_num = np.sum(self.y_gt > 0)
        if (points_num == 0):
            self.y_gt[alpha] = 1.0
        else:
            self.y_gt[alpha] = 1.0/points_num
            self.y_gt = self.y_gt * points_num / (points_num+1)
        return
    
    def drawPredictionDistribution(self):
        plt.figure()
        plt.plot(self.x, self.y_pred)
        plt.title("s-distribution")
        plt.show()
        return
    def drawGroundTrueth(self):
        plt.figure()
        plt.plot(self.x, self.y_gt, "ro")
        plt.title("ground truth")
        plt.show()
    
    def drawScanAngelDistribution(self):
        plt.figure()
        plt.plot(self.x, self.y_scan)
        plt.title("alpha-distribution")
        plt.show()
        return
    
    def calculateGaussian(self, x, mu, sigma, coff):
        y = coff * np.exp(-(x-mu)**2/(2*sigma**2)) / (np.sqrt(2*np.pi)*sigma)
        return y
    
# In[]
if __name__ == "__main__":
    scan_para = [60, 360]
    pdp = avPOMDP(scan_para)
    '''
    input_data = [180, 1]
    pdp.update(input_data, 100)
    input_data = [100, 1]
    pdp.update(input_data, 0)
    input_data = [120, 1]
    pdp.update(input_data, 100)
    pdp.update(None, 150)
    y_gt = pdp.y_gt
    y_pred = pdp.y_pred
    y_scan = pdp.y_scan
    pdp.draw()
    '''
    while (True):
        wait = input("PRESS ENTER TO CONTINUE.")