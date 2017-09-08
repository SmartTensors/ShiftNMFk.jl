function AIC_final(RECON, SILL_AVG, numT, nd)
	siluet= SILL_AVG
	Norm = RECON
	ndata = numT*nd	# nd is the number of detectors, numT is the number of time points (80 in diffNMF example)
	# Finds the optimal number of sources by minimizing AIC
	# Takes as inputs two vectors with the same size (norm and siluet) and a
	# scalar ndata (total number of data points)
	aic_values = zeros(size(Norm))
	for i = 1:length(Norm)
		if siluet[i] > 0.7
			aic_values[i] = 2*i + ndata*log(Norm[i]/ndata)
		else
			aic_values[i] = NaN
		end
	end
	aic_min, nopt = findmin(aic_values)
	#nopt = find(aic_values .== aic_min)
	#nopt = nopt[1]
	return aic_min, nopt
end
