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
        self.points_x_pos = [0]
        self.points_y_pos = [0]
        self.x_scan_range_1 = []
        self.y_scan_range_1 = []
        self.x_scan_range_2 = []
        self.y_scan_range_2 = []
        self.now_time = -1
        self.scan_angle = []
        self.scan_angle_goal = []
        self.gaussian_sigma = gaussian_sigma
        
        self.drawInitialization()
        
    def update(self, input_data, scan_angle_goal):
        self.now_time += 1
        #self.scan_angle.append(scan_angle_goal)
        if (len(self.scan_angle) == 0):
            self.scan_angle.append(90)
        else:
            last_angle = self.scan_angle[-1]
            if (scan_angle_goal >= last_angle):
                if (scan_angle_goal-last_angle > self.max_range-\
                    (scan_angle_goal-last_angle)):
                    self.scan_angle.append((last_angle-5)%self.max_range)
                else:
                    self.scan_angle.append((last_angle+5)%self.max_range)
            else:
                if (last_angle-scan_angle_goal > self.max_range-\
                    (last_angle-scan_angle_goal)):
                    self.scan_angle.append((last_angle+5)%self.max_range)
                else:
                    self.scan_angle.append((last_angle-5)%self.max_range)
        self.scan_angle_goal.append(scan_angle_goal)
        if (input_data != None):
            (alpha, distance) = input_data
            self.addPoint(alpha, distance, self.now_time)
        self.updatePoints()
        self.updatePredictionDistribution()
        self.updateScanAngelDistribution()
        self.updateScanRangeDisplay()
        
    def draw(self):
        self.drawGroundTrueth()
        self.drawPredictionDistribution()
        self.drawScanAngelDistribution()
        self.fig.canvas.draw()
        
    def sampleScanAngle(self):
        alpha = np.random.choice(self.x[:,0], p=self.y_scan[:,0])
        return alpha
        
    def updatePoints(self):
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
        
    def updatePredictionDistribution(self, sigma_expand_cof=1.0):
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
                    coff = 0.1/point[1]
                    self.y_pred[x] += self.calculateGaussian(x, mu, 
                               sigma, coff)
                    
        # Newly added version
        
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
        m = 1.1
        n = 1.0
        epsilon = 1e-3
        last_angle_goal = self.scan_angle_goal[-1]
        last_angle = self.scan_angle[-1]
        
        for i in range(self.max_range):
            distance = min(abs(i-last_angle_goal), abs(self.max_range-(i- \
                           last_angle_goal)))
            coef = min(abs(i-last_angle), abs(self.max_range-(i- \
                           last_angle)))
#==============================================================================
#             coef = 1            
#==============================================================================
            self.y_scan[i] = coef**20 * ((self.y_scan[i]*100)**20) / \
                           ((distance + epsilon)**3)
        
        self.y_scan = self.y_scan / np.sum(self.y_scan)
        return
    
    def updateScanRangeDisplay(self):
        left_bnd = (self.scan_angle[-1] - self.scan_range/2) % self.max_range
        right_bnd = (self.scan_angle[-1] + self.scan_range/2) % \
                     self.max_range
        left_bnd = left_bnd/180 * np.pi
        right_bnd = right_bnd/180 * np.pi
        self.x_scan_range_1 = np.linspace(0, 10*np.cos(left_bnd), 100)
        self.y_scan_range_1 = np.linspace(0, 10*np.sin(left_bnd), 100)
        self.x_scan_range_2 = np.linspace(0, 10*np.cos(right_bnd), 100)
        self.y_scan_range_2 = np.linspace(0, 10*np.sin(right_bnd), 100)
        return
    
    def addPoint(self, alpha, distance, points_added_time):
        self.points.append([alpha, distance, points_added_time, 0])
        self.points_x_pos.append(distance*np.cos(alpha/180*np.pi)) 
        self.points_y_pos.append(distance*np.sin(alpha/180*np.pi))
        points_num = np.sum(self.y_gt > 0)
        if (points_num == 0):
            self.y_gt[alpha] = 1.0
        else:
            self.y_gt[alpha] = 1.0/points_num
            self.y_gt = self.y_gt * points_num / (points_num+1)
        return
    
    def drawInitialization(self):
        plt.ion()
        self.fig = plt.figure()
        self.fig_gt = self.fig.add_subplot(221)
        self.fig_pred = self.fig.add_subplot(222)
        self.fig_scan = self.fig.add_subplot(223)
        self.fig_gt.set_title("Truth")
        self.fig_pred.set_title("Prob of prediction")
        self.fig_scan.set_title("Prob of scan angle")
        self.line_pred, = self.fig_pred.plot(self.x, self.y_pred)
        self.line_gt, = self.fig_gt.plot(self.points_x_pos, 
                                         self.points_y_pos,
                                         ms=3,color='k',marker='o',ls='')
        self.line_range_1, = self.fig_gt.plot(self.x_scan_range_1,
                                              self.y_scan_range_1,
                                              color='g')
        self.line_range_2, = self.fig_gt.plot(self.x_scan_range_2,
                                              self.y_scan_range_2,
                                              color='g')
        self.line_scan, = self.fig_scan.plot(self.x, self.y_scan)
        return
    
    def drawPredictionDistribution(self):
        self.line_pred.set_ydata(self.y_pred)
        self.fig_pred.set_xlim(0, self.max_range)
        self.fig_pred.set_ylim(0, max(self.y_pred) * 1.1)
        #self.fig_pred.title("s-distribution")
        return
        
    def drawGroundTrueth(self):
        self.line_gt.set_data(self.points_x_pos, self.points_y_pos)
        self.line_range_1.set_data(self.x_scan_range_1, self.y_scan_range_1)
        self.line_range_2.set_data(self.x_scan_range_2, self.y_scan_range_2)
        display_range = max(max(max(self.points_x_pos), 
                                -min(self.points_x_pos)), 
            max(max(self.points_y_pos), -min(self.points_y_pos)))
        self.fig_gt.set_xlim(-display_range-1, display_range+1)
        self.fig_gt.set_ylim(-display_range-1, display_range+1)
        '''
        self.fig_gt.set_xlim(-max(max(self.points_x_pos), 
                                  -min(self.points_x_pos)) - 1, 
                               max(max(self.points_x_pos), 
                                  -min(self.points_x_pos)) + 1)
        self.fig_gt.set_ylim(-max(max(self.points_y_pos), 
                                  -min(self.points_y_pos)) - 1, 
                               max(max(self.points_y_pos), 
                                  -min(self.points_y_pos)) + 1)
        '''
        #self.fig_gt.title("ground truth")
        return
    
    def drawScanAngelDistribution(self):
        self.line_scan.set_ydata(self.y_scan)
        self.fig_scan.set_xlim(0, self.max_range)
        self.fig_scan.set_ylim(0, max(self.y_scan) * 1.1)
        #self.fig_scan.title("alpha-distribution")   
        return
    
    def calculateGaussian(self, x, mu, sigma, coff):
        y = coff * np.exp(-(x-mu)**2/(2*sigma**2)) / (np.sqrt(2*np.pi)*sigma)
        return y

            