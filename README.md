# 3D SIFT keypoints utilities
Manipulate and analyse 3D SIFT-Rank keypoints with python.
 
3D SIFT rank keypoints were introduced in http://www.matthewtoews.com/papers/matt_MIA2013.pdf and provide a scale invariant compact representation of 3D images. The tool to produce 3D keypoints is available at http://www.matthewtoews.com/fba/featExtract1.3.zip.

The utility contains:
1. read/write utilities for 3D keypoints and medical volumes
2. filtering for keypoints
3. analysing and manipulating functions for keypoints

This code was written for my thesis and is a compilation of the code most usefull for anyone working with 3D keypoints.

The Python_visualize_keypoints directory contains code to visualize keypoints with 3D Slicer. Use the following commands on 3D Slicer (with python2 interpreter) to visualize keypoints:
<pre><code>
import sys
sys.argv=['',path_of_keypoint_file]
execfile(path_of_visualizeFeatures.py)
</code></pre>

If you have a python 3 interpreter in slicer, use visualizeFeatures3.py instead. It will generate a vtk file that you can then visualize. 

Don't forget to put a r in front of the path to indicate that it is a raw string. ( r"S:\3D-SIFT-keypoints-utilities\Python_visualize_keypoints\visualizeFeatures.py")


File calculateRm.py contains code necessary to calculate the Rm for any keypoint in a brain image. Rm is the ratio of data encoded in a keypoint originating from the mask (brain mask in this example). This file is a good example of how to use the utilities.
