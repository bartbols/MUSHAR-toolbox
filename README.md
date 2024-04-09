## *MUSHAR*:
A MATLAB toolbox for three-dimensional reconstruction and analysis of *MU*scle *SH*ape and *AR*chitecture.

<img src="img/overview.png" alt="MUSHAR-toolbox"
	title="MUSHAR-toolbox" width="100%"/>
    
    
## Features
1. Establish point-to-point correspondence on surfaces and inside volumes through non-rigid registration of distance maps.
1. Reconstruction of group-averaged muscle shape and muscle fibre orientations from magnetic resonance imaging and diffusion tensor imaging data.
    1. Includes code for averaging and interpolating diffusion tensors in the [log-Euclidean domain](https://doi.org/10.1002/mrm.20965)
1. Statistical analysis of local changes in shape and fibre orientations.
1. Visualization of changes in shape and fibre orientations.

## Installation
* Install Matlab (developed and tested in version R2019b)
* Add the MUSHAR-toolbox to the Matlab path.

The following software tools should be installed and made available on the command line:
* [Shapeworks](http://sciinstitute.github.io/ShapeWorks/) (tested in version 6.2.1 - older versions may not be compatible)
* [Elastix](https://elastix.lumc.nl/) (developed and tested in version 4.7)
* [Convert3D](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D)

## Getting started
* Run the demo scripts to  guide you through the main steps.
* Modify the scripts to set up your own analysis.

## Citation
Please cite the following paper when using this toolbox:

Bolsterlee, B., 2022. A new framework for analysis of three-dimensional shape and architecture of human skeletal muscles from in vivo imaging data. Journal of Applied Physiology 132, p712-725.
 [link](http://doi.org/10.1152/japplphysiol.00638.2021)

(Or read the preprint [here](https://doi.org/10.1101/2021.09.08.459536).)






