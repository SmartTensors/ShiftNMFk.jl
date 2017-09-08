module ShiftNMFk
	using JSON #TODO switch to import
	using Clustering
	using Distances
	using Stats
	#using Gadfly
	#using Compose
	using NLopt

	const shiftnmfkdir = splitdir(splitdir(Base.source_path())[1])[1]

	cd(shiftnmfkdir)

	include("execute.jl")
	include("execute_parallel.jl")
	include("execute_using_functionsCos.jl")
	include("ParseTrials.jl")
	include("Plot.jl")
	include("CosClust.jl")

	include("TriangulatePos.jl")
	include("Parallel_Triangle.jl")
	include("LocCluster.jl")

	include("ParseLoc.jl")
	include("AIC.jl")
	include("setdir.jl")

	function test()
		include(joinpath(shiftnmfkdir, "examples", "InsideGrid", "Test.jl"))
		include(joinpath(shiftnmfkdir, "examples", "OutsideGrid", "Test.jl"))
	end
end