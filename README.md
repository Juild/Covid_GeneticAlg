# Covid GeneticAlg
Fitting a complex model of Covid-19 with Genetic Algorithms and Runge Kutta

## Install and build

To build the executable you will need to install the gsl and gsl-blas library: `sudo apt install libgsl-dev`.

Maybe you'll also need cmake to compile the files: `sudo apt install cmake`.

Then, simply navigate to the root directory of the project and run:

    cmake .
    make covidga

You should see a lot of new files being created, the only important one is the executable file `covidga`.
