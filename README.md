# OrganoidChaste

This chaste user project allows the modelling of three-dimensional spherical organoids. 
It implements a three-dimensional vertex model in Chaste, where each cell is modelled as a polyhedron with one apical, one basal and multiple lateral faces.

To install it, you first need to install [Chaste](https://chaste.cs.ox.ac.uk/trac) and then may run this as a user project.
See the [User Projects](https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects) guide page on the Chaste wiki for more information.

The project is compatible with Chaste release 2019.1 and 2021.1.

## Citing us
This code has been developed for the simulation of organoids and has been used both in

*Drozdowski and Schwarz, Phys. Rev. Res. 6, L022045 (2024); DOI: 10.1103/PhysRevResearch.6.L022045*

and in the form of the bubbly vertex model in

*Drozdowski, Kocamese, Boonekamp, Boutros, and Schwarz, arXiv:2411.07141 (2024); DOI: 10.48550/arXiv.2411.07141*

If you use this code please cite the second manuscript (preprint), especially if you use the export to SurfaceEvolver for the bubbly vertex model. 

## Quick start

We use [Chaste Docker](https://github.com/Chaste/chaste-docker) to install Chaste and run all the code inside docker.

1. Install [Docker](https://www.docker.com/).
2. Build the current release version of the Chaste docker image by using on the terminal
  ```
  docker build -t chaste --build-arg TAG=2021.1 https://github.com/oliverdrozdowski/chaste-docker.git
  ```
3. Now you may create a testoutput folder for the simulations and a project folder in which we want to manage the source code of **OrganoidChaste**
  ```
  mkdir testoutput
  mkdir projects
  cd projects
  git clone https://github.com/oliverdrozdowski/OrganoidChaste.git
  ```
4. Now we may start Chaste and bind mount the two folders inside docker by running
  ```
  cd ..
  docker run -it --init -v chaste_data:/home/chaste -v $(pwd)/testoutput:/home/chaste/testoutput -v $(pwd)/projects:/home/chaste/src/projects --name chaste chaste
  ```
5. Inside the container we need to install the computer geometry algorithms library [CGAL](https://www.cgal.org/), as we use it in **OrganoidChaste**. For this we can easily install it via
  ```
  sudo apt-get update
  sudo apt-get install libcgal-dev
  ```
  where we use the password `chaste` to answer the prompt

6. We now move to the `lib` folder and configure CMake to run Chaste with **OrganoidChaste**. In the folder run `ccmake .` press `c`, then `e` again and finally generate the CMake structure by pressing `g`.

7. **OrganoidChaste** should now be set up correctly. You may check this by running a simple Test file from the **OrganoidChaste** test suite, e.g.
  ```make TestMonolayerVertexElement```
  and then
  ```ctest -R TestMonolayerVertexElement```

8. We may leave the Docker container via `exit` and may restart it via `docker start -i chaste`

For more infos on how to create/run projects in Chaste consult the official site with [tutorials and documentations](https://chaste.cs.ox.ac.uk/trac/wiki).

## Creating the Doxygen API
To get the API of **OrganoidChaste**, check the [tutorial](https://chaste.github.io/docs/dev-guides/cmake-build-guide/#other-useful-targets) in the Chaste wiki.

Before, however, you need to teach the original Chaste doxygen configuration to also include **OrganoidChaste**. Assuming **OrganoidChaste** is included in the projects folder (like in the above tutorial), we need to modify the `Doxyfile`in `chaste/src` to also look for **OrganoidChaste**.

In the `Doxyfile` replace

```
#---------------------------------------------------------------------------
# configuration options related to the input files
#---------------------------------------------------------------------------
INPUT                  =    heart/src \
                            global/src \
                            io/src \
                            linalg/src \
                            mesh/src \
                            ode/src \ 
                            pde/src \
                            continuum_mechanics/src \
                            cell_based/src \
                            crypt/src \
                            notforrelease/src \
                            notforrelease_cell_based/src \
                            notforrelease_lung/src 
```

with 
```
#---------------------------------------------------------------------------
# configuration options related to the input files
#---------------------------------------------------------------------------
INPUT                  =    heart/src \
                            global/src \
                            io/src \
                            linalg/src \
                            mesh/src \
                            ode/src \ 
                            pde/src \
                            continuum_mechanics/src \
                            cell_based/src \
                            crypt/src \
                            notforrelease/src \
                            notforrelease_cell_based/src \
                            notforrelease_lung/src \
                            projects/OrganoidChaste/src
```
and use the usual `ccmake` and then `make doxygen` to create the Chaste doxygen API, including **OrganoidChaste**.

## Example application: spherical shells
To create spherical shells, which can be imported into surface evolver for the bubbly vertex model analysis, the application *RandomizeHomogeneousSphere* can be used.

For this change to the `lib` folder and run `make MinimizeHomogeneousSphere`. After compilation you may run `projects/OrganoidChaste/apps/MinimizeHomogeneousSphere -h` to recieve help on the possible command line arguments that configure the simulation.