# Coordinates

Origin at the lower mounting hole of the lower fiber coupler.
X axis is right, Y axis is down.

# Raw measurements

When a distance is specified as two numbers, they are the inner and outer measurements
so the center distance should be the average of the two. All units are inch.

* Picture taken with 2.5x zoom which for iPad mini 6 rear camera is 73mm of 35mm equivalent focal length. Image size is 4032x3024 pixels.
* The fiber launcher mounting hole was recessed by 0.230.
* Distance between the two mounting screws on a mirror mount: 0.365, 0.478
* X distance between two mirror mounts (point D/E or point F/G in picture): 0.681, 0.794
* Distance between fiber coupler mounting holes (point A/B): 0.590, 0.766
* Distance between the two mounting screws on the fiber coupler mount: 1.735, 1.964
* Width of the long slot: 0.130
* Distance between the slot center and the closer one of the mounting holes in fron of where the 399 laser should be: 2.635, 2.855

# Points

|Point name|Pixel position (x,y)|
|----------|----------------|
|O (Origin)|(2085, 1608)|
|A|(2094, 1084)|
|B|(2283, 1088)|
|C|(2294, 562)|
|D|(1431, 1037)|
|E|(1642, 1043)|
|F|(1429, 1158)|
|G|(1640, 1163)|
|H|(2029, 750)|
|I|(2027, 870)|
|J|(1411, 858)|
|K|(1291, 855)|
|L|(959, 888)|
|M|(574, 877)|
|N|(1575, 1626)|
|P|(1573, 1748)|
|Q|(978, 1699)|
|R|(302, 1689)|
|S|(973, 1661)|
|T|(305, 1650)|


# Constraints

* `O = (0, 0, -0.230)`
* `A = (0, -1.850, -0.230)`
* `B = (0.678, -1.850, -0.230)`
* `C = (0.678, -3.700, -0.230)`
* `D - E = (-0.738, 0)`
* `D - F = (0, -0.422)`
* `F - G = (-0.738, 0)`
* `H - I = (0, -0.422)`
* `K - J = (-0.422, 0)`
* `Y(L) - (Y(T) + Y(S) + Y(R) + Y(Q)) / 4 = -2.745`
* `Y(M) - Y(L) = 0`
* `Y(T) - Y(S) = 0`
* `Y(R) - Y(Q) = 0`
* `Y(T) + Y(S) - Y(R) - Y(Q) = -0.260`
