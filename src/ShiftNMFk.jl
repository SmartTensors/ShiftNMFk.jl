module ShiftNMFk
	using JSON #TODO switch to import
	using Clustering
	using Distances
	using Stats
	#using Gadfly
	#using Compose
	using NLopt

	include("ShiftNMF2.jl");
	include("Parallel_ShiftNMF.jl");
	include("ParseTrials.jl");
	include("Plot.jl");
	include("ShiftNMFk_CosClust.jl");
	include("ShiftNMFk_using_functionsCos.jl");

	include("TriangulatePos.jl");
	include("Parallel_Triangle.jl");
	include("LocCluster.jl");

	include("ParseLoc.jl");
	include("AIC_final.jl");
end