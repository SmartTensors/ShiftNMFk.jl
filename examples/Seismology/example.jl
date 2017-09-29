using Gadfly
using Cairo
using Compose

addprocs()

@everywhere import JSON
@everywhere import Clustering 
@everywhere import Stats
@everywhere import Gadfly 
@everywhere import Compose 
@everywhere import NLopt


@everywhere function setdir(dir)	
	if isdir(dir)
		cd(dir)
	end
end

@everywhere function setdir()
	dir = remotecall_fetch( ()->pwd(), 1)
	setdir(dir)
end


@everywhere setdir()
@everywhere @show pwd()

@everywhere include("../src/ShiftNMFkMD.jl");   # including the module on all cores
@everywhere using ShiftNMFk    			# using the module ShiftNMFk on all cores

X      = readcsv("Observations_508_Normed_cut.csv");	#Inputing the observation matrix of the desired example Xnm-> 20 x 40000; where 20 is the number of sensors (genomes)  and 40,000  is the number of the time points (type of mutations)
X      = X';
micPos = readcsv("MicPosition_508_cut.csv");            # The coordinates of the detectors in the grid

# When you run the Example.jl in the NewEx folder, the script creates a Reconstruction Plots folder where each one of the sensors are 
# plotted( original observation and reconstruction) with the weights of each source that make up the seen waveform in a bar plot.
# Additionally you can set the maxSource  to a vector that will test only the guesses of sources you need.
# Example:  maxSource =3 will run it for sources 1, 2, 3.
# maxSource = [ 3 5 6] will run it for sources 3, 5, 6.




maxSource  = 3;								# max number of sources, or [] specific number of sources
globalIter = 30000;								        # 1000 NMF runs for each guess of a sources								
nmfIter    = 20;								# 80,000 max number of iterations for each source.
locIter    = 10;								# 1000 minimizations are performed to find the location
numT       = size(X,2); 							# The numbes of sample points that make up our signals (time samples) 
nd         = size(micPos,1);							# The number of detectors in our grid

# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

shiftNMFk(X, maxSource, globalIter, nmfIter);                            # 1		

#ShiftNMFk.ShiftNMFk_CosClust(X, maxSource);	             		 # 2	                # if stops before we can continue from here		 

Sil, Norm     = ShiftNMFk.Plot(X, maxSource);                		 # 3			# Plots the Norm and Silhouette Value graph and returns both 

aic_min, nopt = ShiftNMFk.AIC_final( Norm, Sil, numT, nd)    		 # 4			# nopt returs the correct number of sources in the experiment. 

W,H,T,Tstd    = ShiftNMFk.ResultsForNumSources(nopt);       		 # 5			#returns the mixing(W), signal(H), and delay(T) matricies for the number
								       				#of sources nopt. Also returs sigma of the delays(Tstd)

SourcePositions = ShiftNMFk.FindLocations(W, T, Tstd, micPos, locIter);	# 6     		#Returns the coordinates of the sources. 

using Gadfly
include("../src/PlotReconstruction.jl"); # Including function that plots sensors
PlotReconstruction(X,nopt);


exit()

