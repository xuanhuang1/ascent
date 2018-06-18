# Ascent Ultroids Example

The ultroids application generates a synthetic data set based on Craig Reynolds boids simulations.

The utlroids application comes with a serial (ultroids_ser) and distributed memory parallel version (ultroids_par) that will automatically distribute a uniform data set across ranks.

# Options
 - `--time_steps=10` The number of time steps to generate
 - `--time_delta=0.5` The amount of time to advance per time step.
 - `--dims=x,y,z` The total spatial dimensions of the data set



Sample Runs:

./utlroids_ser

mpiexec -n 2 ./utroids_par

For more info, please visit:

http://ascent.readthedocs.io/en/latest/ExampleIntegrations.html
