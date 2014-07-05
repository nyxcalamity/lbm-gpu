Lattice-Boltzmann Method Using GPGPU CUDA Platform
==================================================
---------------------------------------------------------------------------------------------
Description
-----------
The project is an open-source [GPGPU] implementation of [Lattice-Boltzmann Method] (LBM), a [computational fluid dynamics] (CFD) method for fluid simulation, that solves a [lid-driven cavity problem] on a D3Q19 lattice grid.

[GPGPU]:http://en.wikipedia.org/wiki/General-purpose_computing_on_graphics_processing_units
[Lattice-Boltzmann Method]:http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods
[computational fluid dynamics]:http://en.wikipedia.org/wiki/Computational_fluid_dynamics
[lid-driven cavity problem]:http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem


Technical details
-----------------
Due to a local nature of computations performed during lattice grid update the method is highly parallelizable and can be scaled to virtually the same amount of compute units as the number of cells being used for the domain. Since modern GPUs have thousands of execution cores and show a tendecny of increasing that number, they are the perfect candidate for LBM parallelized code to run on. This project utilizes [CUDA] platform due to the fact that comparing to it's competitor, [OpenCL], it has a wider functionality and higher data transfer rate with virtually the same computaitonal throughput. However CUDA comes with a cost of locking developers to use NVidia GPUs, which is irrelevant for the purposes of this project.

[CUDA]:http://en.wikipedia.org/wiki/CUDA
[OpenCL]:http://en.wikipedia.org/wiki/OpenCL


Goals
-----
During the project implementation next goals were accounted for:
- high performance and efficiency of LBM solver
- high scalability of code to various NVidia GPU architectures
- maintaiability and clarity of the code


Technical Prerequisites
-----------------------
In case that the reader is not familiar with the GPGPU programming models or the innerworkings of GPU hardware it is highly recommended to skim through [NVidia programming guide] and [NVidia GPU architectures]. It is also recommended to have a general understanding of LBM solver principles.

[NVidia programming guide]:http://docs.nvidia.com/cuda/cuda-c-programming-guide/
[NVidia GPU architectures]:https://developer.nvidia.com/key-technologies


Implementation
--------------
The project was implemented in [C] utilizing [CUDA 5.5 Toolkit] and consists of two aligned implementations of the LBM solver: CPU and GPU. GPU functionality is decoupled from CPU code and is enclosed in files with `_gpu.cu` or `_gpu.cuh` endings. General project structure is as follows:

- `main.c` - main funciton that triggers simulation routines
- `lbm_model.h` - problem/LBM specific constants and validation methods
- GPU:
    - `initialization_gpu.h` - GPU memory initialization and freeing
    - `lbm_solver_gpu.h` - LBM solver that encompasses streaming, collision and boundary treatment
    - `cell_computation_gpu.cuh` - decoupled local cell computations
- CPU:
    - `initialization.h` - CLI and config file parsing
    - `streaming.h` - streaming computations
    - `collision.h` - collision computations
    - `boundary.h` - boundary treatment
    - `cell_computation.h` - local cell computations
- `visualization.h` - visualization of fields

[C]:https://en.wikipedia.org/wiki/C_(programming_language)
[CUDA 5.5 Toolkit]:https://developer.nvidia.com/cuda-toolkit-55-archive


Compatability
-------------
The code is compatible with GPUs of [compute capability] `2.0` and higher and NVidia CUDA Toolkits of version `4.0` and higher. 

The project was tested against NVidia [GeForce GTX 660M] (CC `3.0`) and [GeForce GTX 460] (CC `2.1`). Development was performed solely on Linux system, however, there should be no problems with running it on windows.

[compute capability]:http://docs.nvidia.com/cuda/cuda-c-programming-guide/#compute-capability
[GeForce GTX 660M]:http://www.geforce.com/hardware/notebook-gpus/geforce-gtx-660m/specifications
[GeForce GTX 460]:http://www.geforce.com/hardware/desktop-gpus/geforce-gtx-460/specifications


Building and running
--------------------
These instructions are aimed at linux users who have [CUDA enabled GPUs] with compute capability 2.0+ and who have already [installed] and enabled gpu device drivers. It is also expected that the reader went through [NVidia getting started guide] and installed CUDA Toolkit `4.0` or newer. 

Other dependencies:

- [gcc] version `4.8.2+`
- [GNU Make] version `3.81+`
- [git] version `1.9.1+`

1. Clone the project from gihub repository:

        git clone https://github.com/nyxcalamity/lbm-gpu.git <project-dir>


2. Navigate to `<project-dir>` directory and run:

        make

3. Adjust grid size or physical properties of the problem in the configuration file located in `<project-dir>/data/lbm.dat`.
4. Run the project using next command:

        <project-dir>/build/lbm-sim -help

5. Read the help message and run the actual simulation as follows:

        <project-dir>/build/lbm-sim <project-dir>/data/lbm.dat -gpu

[CUDA enabled GPUs]:https://developer.nvidia.com/cuda-gpus
[installed]:https://help.ubuntu.com/community/BinaryDriverHowto/Nvidia
[NVidia getting started guide]:http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/
[gcc]:https://gcc.gnu.org/
[GNU Make]:http://www.gnu.org/software/make/
[git]:http://git-scm.com/


Known-Issues
------------
There are several known issues with the project which do not affect it's performance or the resulting simulation:

- due to optimization of boundary treatment code we reduced 57 checking branches to just 22 at a cost of exchanging probability distribution functions between boundary cells at the edges
- an unknown rounding error happens during visualization which might change a minority of values by not more than 0.000001