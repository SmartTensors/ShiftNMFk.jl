import Gadfly

function PlotReconstruction(Observation, nSources)
	if !isdir("./Results")
		warn("Results directory does not exist!")
		mkdir("./Results")
	end
	W = readcsv("./Results/W$nSources.csv")
	for i=1:size(W, 1)
			sumW=sum(W[[i], :])
			W[[i], :] = W[[i], :]./sumW
	end
	W = 100*W
	X = Observation
	Rec = readcsv("./Results/Reconstruction$nSources.csv")
		#normalize X & Rec
		for i=1:size(X, 1)
			sumX=sum(X[[i], :])
			sumRec=sum(Rec[i, :])
			X[[i], :] = X[[i], :]./sumX
			Rec[[i], :] = Rec[[i], :]./sumRec
		end
	if isdir("./ReconstructionPlots") == false
		mkdir("./ReconstructionPlots")
	end
	x = 1:size(X, 2)
	for nMics = 1:size(X, 1)
		RecCorr = cor(X[[nMics], :]', Rec[[nMics], :]')
		ResultPlot = Gadfly.plot(
		    Gadfly.layer(x=x, y=X[[nMics], :], Gadfly.Geom.line, Gadfly.Theme(default_color=parse(Colors.Colorant, "green"))),
			Gadfly.layer(x=x, y=Rec[[nMics], :], Gadfly.Geom.point, Gadfly.Theme(default_color=parse(Colors.Colorant, "red"))),
			Gadfly.Guide.xlabel("Sample Points"), Gadfly.Guide.ylabel("Normalized Signal Amplitude"), Gadfly.Guide.title("Sensor: $nMics With $RecCorr Correlation"), Gadfly.Guide.manual_color_key("Legend",
			["Reconstruction", "Original Observatons"], ["red", "green", ])
		)
		BarPlot = Gadfly.plot(x=1:nSources, y = W[[nMics], :], Gadfly.Geom.bar, Gadfly.Guide.xlabel("Sources"), Gadfly.Guide.ylabel("Percentage of source in Sensor Observation"), Gadfly.Guide.title("Percentages of each source seen by the sensor"))
		Both = Gadfly.hstack(ResultPlot, BarPlot)
		Gadfly.draw(Gadfly.SVG("./ReconstructionPlots/Sensor$nMics.With$RecCorr.Correlation.svg", 40Gadfly.cm, 20Gadfly.cm), Both)
	end
	return Rec
end
