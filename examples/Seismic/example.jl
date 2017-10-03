@everywhere reload("ShiftNMFk")
@everywhere ShiftNMFk.setdir()

#Inputing the observation matrix of the desired example Xnm-> 20 x 40000; where 20 is the number of sensors (genomes)  and 40,000  is the number of the time points (type of mutations)
X = readcsv("examples/Seismic/Observations_508_Normed_cut.csv")
X = X';
# The coordinates of the detectors in the grid
micPos = readcsv("examples/Seismic/MicPosition_508_cut.csv")

maxSource  = 3;								# max number of sources, or [] specific number of sources
globalIter = 30000;								        # 1000 NMF runs for each guess of a sources
nmfIter    = 20;								# 80,000 max number of iterations for each source.
locIter    = 10;								# 1000 minimizations are performed to find the location
numT       = size(X,2); 							# The number of sample points that make up our signals (time samples)
nd         = size(micPos,1);							# The number of detectors in our grid

ShiftNMFk.execute(X, maxSource, globalIter, nmfIter)

# ShiftNMFk.ShiftNMFk_CosClust(X, maxSource)

Sil, Norm     = ShiftNMFk.Plot(X, maxSource)

aic_min, nopt = ShiftNMFk.AIC( Norm, Sil, numT, nd)

W,H,T,Tstd    = ShiftNMFk.ResultsForNumSources(nopt)

SourcePositions = ShiftNMFk.FindLocations(W, T, Tstd, micPos, locIter)

include("../src/PlotReconstruction.jl") # Including function that plots sensors
PlotReconstruction(X,nopt)
