Lattice-Boltzmann Method Using GPGPU CUDA Platform
==================================================
---------------------------------------------------------------------------------------------
Description
-----------
The project is an open-source [GPGPU] implementation of [Lattice-Boltzmann Method] (LBM), a [computational fluid dynamics] (CFD) method for fluid simulation, that solves a [lid-driven cavity problem] on a D3Q19 lattice grid.


Technical details
-----------------
Due to a local nature of computations performed during lattice grid update the method is highly parallelizable and can be scaled to virtually the same amount of compute units as the number of cells being used for the domain. Since modern [GPUs] have thousands of execution cores and show a tendecny of increasing that number, they are the perfect candidate for LBM parallelized code to run on. This project utilizes CUDA platform due to the fact that comparing to it's competitor, OpenCL, it has a wider functionality and higher data transfer rate with virtually the same computaitonal throughput. However CUDA comes with a cost of locking developers to use NVidia GPUs, which is irrelevant for the purposes of this project.


Goals
-----
During the project implementation next goals were accounted for:
- high performance and efficiency of LBM solver
- high scalability of code to various NVidia GPU architectures
- maintaiability and clarity of the code


Technical Prerequisites
-----------------------
In case that the reader is not familiar with the GPGPU programming models or the innerworkings of GPU hardware it is highly recommended to skim through [nvidia programming guide] and [NVidia GPU architectures]. It is also recommended to have a general understanding of LBM solver principles.


Implementation
--------------
The project was implemented in [C] and utilized [CUDA 5.5 Toolkit]. GPU code is decoupled from CPU code and is included in files with '_gpu.cu' or '_gpu.cuh' endings. 


//TODO:include file structure and etc


Within those files GPU code is united so that all the processing is done within a GPU function. 

Compatability
-------------
The code is compatible with GPUs of [compute capability] 2.0 and higher and NVidia CUDA Toolkits of version 4.0 and higher. 

The project was tested against NVidia GeForce GTX 660M (CC 3.0) and GeForce GTX 460 (CC 2.1). Development was performed solely on Linux system, however, there should be no major problems with running it on windows.


Building and running
--------------------
These instructions are aimed at linux users who have [CUDA enabled GPUs] with compute capability 2.0+ and who have already [installed] and enabled gpu device drivers. It is also expected that the reader went through [NVidia getting started guide] and installed CUDA 4.0 Toolkit or newer. 

1. Clone the project from gihub repository:

    git clone https://github.com/nyxcalamity/lbm-gpu.git <project-dir>


2. Navigate to `<project-dir>` directory and run:

    make

3. Adjust grid size or physical properties of the problem in the configuration file located in `<project-dir>/data/lbm.dat`.
4. Run the project using next command:

    <project-dir>/build/lbm-sim -help

5. Read the help message and run the actual simulation as follows:

    <project-dir>/build/lbm-sim <project-dir>/data/lbm.dat -gpu


Known-Issues
------------
There are several known issues with the project which do not affect it's performance or the resulting simulation:
- due to optimization of boundary treatment code we reduced 57 checking branches to just 22 at a cost of exchanging probability distribution functions between boundary cells at the edges


[computational fluid dynamics]:http://en.wikipedia.org/wiki/Computational_fluid_dynamics
[lid-driven cavity problem]:http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem
[nvidia programming guide]:http://docs.nvidia.com/cuda/cuda-c-programming-guide/
[NVidia GPU architectures]:https://developer.nvidia.com/key-technologies
[C]:https://en.wikipedia.org/wiki/C_(programming_language)
[CUDA 5.5 Toolkit]:https://developer.nvidia.com/cuda-toolkit-55-archive