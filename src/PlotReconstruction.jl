function PlotReconstruction(Observation, nSources)
	if !isdir("./Results")
		println("ERROR: Results directory does not exist!")
	else
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
			ResultPlot = plot(
			 layer(x=x, y=X[[nMics], :], Geom.line, Theme(default_color=parse(Colors.Colorant, "green"))),
			 layer(x=x, y=Rec[[nMics], :], Geom.point, Theme(default_color=parse(Colors.Colorant, "red"))),
			 Guide.xlabel("Sample Points"), Guide.ylabel("Normalized Signal Amplitude"), Guide.title("Sensor: $nMics With $RecCorr Correlation"), Guide.manual_color_key("Legend",
				["Reconstruction", "Original Observatons"], ["red", "green", ])
			)
			BarPlot = plot(x=1:nSources, y = W[[nMics], :], Geom.bar, Guide.xlabel("Sources"), Guide.ylabel("Percentage of source in Sensor Observation"), Guide.title("Percentages of each source seen by the sensor"))
			Both = hstack(ResultPlot, BarPlot)
			draw(SVG("./ReconstructionPlots/Sensor$nMics.With$RecCorr.Correlation.svg", 40cm, 20cm), Both)
		end
	end
	return Rec
end
