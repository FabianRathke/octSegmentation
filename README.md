This package provides a toolbox for segmenting OCT retina scans, 2-D or 3-D. The toolbox is written in Matlab and C. Installation instructions and a detailed documentation may be found in documentation.pdf. The corresponding paper can be found here: http://www.sciencedirect.com/science/article/pii/S1361841514000449.

----------------
INSTALLATION
-----------------

1) Unzip the code into any directory

2) In Matlab set the bash variable OCT_CODE_DIR to the directory from step 1) via

setenv('OCT_CODE_DIR',xxx);

where xxx is a string pointing to that directory

3) call function compileMex to compile all C-Functions


-----------
EXAMPLES
-----------

The package provides two models (datafiles/modelFiles) for circular scans and 3-D volumes, that where trained using our labeled ground truth.
Unfortunately we are not allow to publish any of these scans. For 2-D we have one demoscan provided, for 3-D we have demo scripts that use an external dataset.

2-D
-----

1) run useExampleScan2D.m to obtain a segmentation for the circular scan that is included in the package
2) run octGUI: Load the matfile circularScan.mat (datafiles/exampleScans) and the model file (datafiles/modelFiles) and run the model by clicking on 'Segment'


3-D
------

1) Download the dataset of Srinivasan et al. from http://people.duke.edu/~sf59/Srinivasan_BOE_2014_dataset.htm
2) Run the script predSrinivasan.m for a demonstration how to segment volumes (uses the model provided in the package)
3) Run the script trainSrinivasan.m for a demonstration how to train a model for volumes using as ground truth the predictions from 2)


If you use this software in your publication, please cite the following publication:

"Probabilistic Intra-Retinal Layer Segmentation in 3-D OCT Images Using Global Shape Regularization"
F Rathke, S Schmidt, C Schnörr, Medical Image Analysis

The MIT License (MIT)

Copyright (c) 2014, Fabian Rathke

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

