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

The output beam needs to go through point (1479, 2738)

# Constraints

1. `O = (0, 0, -0.230)`
2. `A = (0, -1.850, -0.230)`
3. `B = (0.678, -1.850, -0.230)`
4. `C = (0.678, -3.700, -0.230)`
5. `D - E = (-0.738, 0)`
6. `D - F = (0, -0.422)`
7. `F - G = (-0.738, 0)`
8. `H - I = (0, -0.422)`
9. `K - J = (-0.422, 0)`
10. `Y(L) - (Y(T) + Y(S) + Y(R) + Y(Q)) / 4 = -2.745`
11. `Y(M) - Y(L) = 0`
12. `Y(T) - Y(S) = 0`
13. `Y(R) - Y(Q) = 0`
14. `Y(T) + Y(S) - Y(R) - Y(Q) = -0.260`

# Processing

The frame size of the 35mm film is 36x24mm. So with a equivalent focal length of 73mm,
and a image size of 4032x3024, (assuming the equivalence is on the diagonal of the image),
the effective focal length should be 8503 pixels.

The transformation from the coordinate on the board, to the coordinate on the picture is
an arbitrary 3D rotation and translation followed by a point projection.
For the actual fit, the center of the image needs to be relabelled as pixel 0, 0.

Since we know all the image coordinates, it's easier/more consistent to reverse transform
the points (given the known Z position in real space) back to the expected x, y position
in real space to do our fit using the constraints above. The fit shows that
the tilt of the camera was (in degree)

* Yaw: -5.19
* Pitch: 4.65
* Roll: 1.01

The camera center is shifted from the origin I defined by

* X: 0.258 inch
* Y: 0.339 inch

And the height of the cemera was 29.8 inch.

# Machining related positions (inch)

* Bottom left mounting hole (the top right hole along the slot): (-3.556, -0.148)
* Bottom right mounting hole (point N): (-1.801, 0.090)
* Top left mounting hole (point L): (-3.978, -2.431)
* Top right mounting hole (point G): (-1.603, -1.528)

* Horizontal beam goes through (center of the lower fiber coupler): (0, -0.925)
* Vertical beam goes through (slightly shifted from a screw hole
  next to the mirror we want to hit): (-2.066, 3.976)
