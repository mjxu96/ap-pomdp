#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 14:04:40 2018

@author: minjun
"""

import pomdp

if __name__ == "__main__":
    scan_para = [60, 360]
    pdp = pomdp.avPOMDP(scan_para)
    while (True):
        input_string = input(\
        "Input the angle and distance of a point or just type ENTER.\n")
        if (input_string == ""):
            input_data = None
        else:
            input_tmp = input_string.split()
            if (len(input_tmp) != 2):
                print("Invalid input.")
                continue;
            else:
                input_data = [int(input_tmp[0]), float(input_tmp[1])]
        alpha = pdp.sampleScanAngle()
        pdp.update(input_data, alpha)
        pdp.draw()
        print("present scan angle is:")
        print(alpha)
            