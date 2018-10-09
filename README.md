# Toy demo of POMDP model for Active Perception
This is a toy demo of POMDP model for Active Perception

## Get Started
The following packages will be needed.
```
numpy
matplotlib
```

## Run
Type in 
```
$ python ./main.py to run this toy demo
```

Then you will see the following:

```
Input the angle and distance of a point or just type ENTER.
```
If you want a new point, which denotes a new object in the view, to appear, 
then type in its polar coordinate parameters. Otherwise just press ENTER.
For example, if you want a point, the polar coordinates of which are 
(200 degree, 3 meters), to appear, just type
```
200 3
```
then you will see in the upper-left graph, there is a black point appearing.

Everytime you input a point or just press ENTER, the demo will sample a scan angle
from the distrition of scan angle. If some points are in the scan range, then the
distribution of prediction of these points will update, shown in the upper-right
graph. After the update, the distribution of scan angle will also be updated in 
the lower-left graph to maximize the requirement of detecting as much close 
points as possible.

## Author
Minjun Xu  minjunxu@andrew.cmu.edu
