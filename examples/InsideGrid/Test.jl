@everywhere reload("ShiftNMFk")
@everywhere ShiftNMFk.setdir()

cd(joinpath(ShiftNMFk.shiftnmfkdir, "examples", "InsideGrid"))

X = readcsv("./InputInGrid/Observation.csv")		#Inputing the observation matrix of the desired example
micPos = readcsv("./InputInGrid/MicPosition.csv") # The coordinates of the detectors in the grid
maxSource = 5										#5 max number of sources
globalIter = 	10								# 1000 NMF runs for each guess of a sources
nmfIter = 100									# 80,000 max number of iterations for each source.
locIter = 10										# 1000 minimizations are performed to find the location
numT = 180 										#The numbers of sample points that make up our signals (time samples)
nd = 16											# The number of detectors in our grid

srand(2015)

ShiftNMFk.execute(X, maxSource, globalIter, nmfIter)

Sil, Norm = ShiftNMFk.Plot(X, maxSource)                        #Plots the Norm and Silhouette Value graph and returns both

aic_min, nopt = ShiftNMFk.AIC( Norm, Sil, numT, nd)       # nopt returns the correct number of sources in the experiment.

W,H,T,Tstd = ShiftNMFk.ResultsForNumSources(nopt)           #returns the mixing(W), signal(H), and delay(T) matrices for the number
                                                    #of sources nopt. Also returns sigma of the delays(Tstd)

SourcePositions = ShiftNMFk.FindLocations(W, T, Tstd, micPos, locIter)  #Returns the coordinates of the sources

include("plotAnsPanel.jl")

cd(ShiftNMFk.shiftnmfkdir)