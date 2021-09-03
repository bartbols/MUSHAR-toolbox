## MUSHAR: A MATLAB toolbox for three-dimensional reconstruction and analysis of muscle shape and architecture.
<img src="img/overview.png" alt="MUSHAR-toolbox"
	title="MUSHAR-toolbox" width="100%"/>
    
    
## Features
1. Establish point-to-point correspondence on surface and inside volumes through non-rigid registration of distance maps.
1. Reconstruction of group-averaged muscle shape and muscle fibre orientations from magnetic resonance imaging and diffusion tensor imaging data.
    1. Includes code for averaging and interpolating diffusion tensors in the [log-Euclidean domain](https://doi.org/10.1002/mrm.20965)
1. Statistical analysis of local changes in shape and architecture (pennation angles).
1. Visualization of changes in shape and fibre orientations.

## Installation
* Install Matlab (developed and tested in version R2019b)
* Add the MUSHAR-toolbox to the Matlab path.

The following software tools should be installed and made available on the command line:
* [Shapeworks](http://sciinstitute.github.io/ShapeWorks/) (developed and tested in version 6.0.0-RC9)
* [Elastix](https://elastix.lumc.nl/) (developed and tested in version 4.7)

## Getting started
The DEMO scripts will guide you through the main steps involved. Modify these scripts to customize your own analysis.




