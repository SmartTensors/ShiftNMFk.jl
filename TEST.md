Module ShiftNMFk

The ShiftNMFk functions are defined in the files in `src`.
The ShiftNMFk functions perform the entire ShiftNMFk analysis and find the locations of the sources.

Two `examples` are provided.
Both feature a grid of 16 detectors (4x4) and 3 sources.
In one example the sources are inside the grid while in the other they are outside.
These examples are the same found in our ShiftNMFk paper.

To install the required packages. (In the future this will be done by the REQUIRE file).

``` julia
include("deps/build.jl")
reload("ShiftNMFk")
ShiftNMFk.test()
```

To test ShiftNMFk:

``` julia
reload("ShiftNMFk")
ShiftNMFk.test()
```

Then run one of the Test.jl files in OutsideGrid and InsideGrid.

``` julia
include("examples/OutsideGrid/Test.jl")
include("examples/InsideGrid/Test.jl")
```

The Test.jl files are commented in a way that explains the ShifNMFk functions and their inputs.
Initially only a small number of runs and iterations are performed to see if the code runs alright.
These variables can be changed within the Test.jl files for running the substantial runs needed to reproduce the examples in our paper.
Running the example creates several folders with results and plots the Silhouette and Norm plot as well as the locations of the sources.

Run plotAnsPanel.jl for a graph that compares the ShiftNMFk results with the original signals used to create the example.
