function ResultsForNumSources(N);
	N = N[1];
	if isdir("./Results") == false
		println("ERROR: Results directory does not exist!")
		throw(ArgumentError("ERROR: Results directory does not exist!"))
	else
		W = readcsv("./Results/W$N.csv");
		H = readcsv("./Results/H$N.csv");
		T = readcsv("./Results/T$N.csv");
		Tstd = readcsv("./Results/Tstd$N.csv");
	end
	return W, H, T, Tstd;
end
# Function definitions
function cluster_NMF_solutions(allProcesses, clusterRepeatMax)
	numberOfPoints = size(allProcesses, 1);
	numberOfProcesses = size(allProcesses, 2);
	globalIter =  size(allProcesses, 3);
	centroids = allProcesses[:, :, 1];
	idx = zeros(Int, numberOfProcesses, globalIter);
	for clusterIt = 1 : clusterRepeatMax
		#println(clusterIt)
		for globalIterID = 1 : globalIter
			processesTaken = zeros(numberOfProcesses, 1);
			centroidsTaken = zeros(numberOfProcesses, 1);
			for currentProcessID = 1 : numberOfProcesses
				distMatrix = ones(numberOfProcesses, numberOfProcesses) * 100;
				for processID = 1 : numberOfProcesses
					for centroidID = 1 : numberOfProcesses
						if ((centroidsTaken[centroidID] == 0) && (processesTaken[processID] == 0))
							distMatrix[processID, centroidID] = cosine_dist(allProcesses[:, processID, globalIterID], centroids[:,centroidID]);
						end
					end
				end
				minProcess,minCentroid = ind2sub(size(distMatrix), indmin(distMatrix));
				processesTaken[minProcess] = 1;
				centroidsTaken[minCentroid] = 1;
				idx[minProcess, globalIterID] = minCentroid;
			end
		end
		centroids = zeros(numberOfPoints, numberOfProcesses);
		for centroidID = 1 : numberOfProcesses
			for globalIterID = 1 : globalIter
				centroids[:, centroidID] = centroids[:, centroidID] + allProcesses[:, findin(idx[:, globalIterID], centroidID), globalIterID];
			end
		end
		centroids = centroids ./ globalIter;
	end
	return idx;
end
function final_processes_and_mixtures(allProcesses, allMixtures, allDelays, idx)
	numberOfSamples = size(allMixtures,2);
	numberOfPoints = size(allProcesses, 1);
	numberOfProcesses = size(allProcesses, 2);
	globalIter =  size(allProcesses, 3);
	idx_r = vec(reshape(idx, numberOfProcesses * globalIter, 1));
	allProcesses_r = reshape(allProcesses, numberOfPoints, numberOfProcesses * globalIter);
	#allMixtures_r = reshape(allMixtures, numberOfProcesses * globalIter, numberOfSamples);
	allMixtures_r = Array{Float64}(numberOfProcesses * globalIter, numberOfSamples);
	allDelays_r = Array{Float64}(numberOfProcesses * globalIter, numberOfSamples);
	for i=1:size(allMixtures,3)
		if i==1
			allMixtures_r = allMixtures[:,:,1];
		else
			allMixtures_r = [allMixtures_r ; allMixtures[:,:,i]];
		end
	end
	for i=1:size(allDelays,3)
		if i==1
			allDelays_r = allDelays[:,:,1];
		else
			allDelays_r = [allDelays_r ; allDelays[:,:,i]];
		end
	end
	allProcessesDist = pairwise(CosineDist(), allProcesses_r);
	stabilityProcesses = silhouettes(idx_r, vec(repmat([globalIter], numberOfProcesses, 1)), allProcessesDist);
	avgStabilityProcesses = zeros(numberOfProcesses, 1);
	processes = zeros(numberOfPoints, numberOfProcesses);
	mixtures = zeros(numberOfProcesses, numberOfSamples);
	delays = zeros(numberOfProcesses, numberOfSamples);			#NEW for ShiftNMF
	for i = 1 : numberOfProcesses
		avgStabilityProcesses[i] = mean(stabilityProcesses[findin(idx_r,i)]);
		processes[:, i] = mean(allProcesses_r[:, findin(idx_r,i)],2);
		mixtures[[i], :] = mean(allMixtures_r[findin(idx_r,i),:],1);
		delays[[i], :] = mean(allDelays_r[findin(idx_r,i),:],1);
	end
	return processes, mixtures, delays, avgStabilityProcesses
end
function ShiftNMFk_CosCluster(X, maxSource)
	if typeof(maxSource) == Int;
		Trials = collect(1:maxSource);
	else
		Trials = maxSource;
	end
	inputMatrix = X;
	inputMatrix=inputMatrix';
	numberOfPoints = size(inputMatrix, 1);
	numberOfSamples = size(inputMatrix, 2);
	VarChange = Array{Any}(length(Trials), 6);
	ALLCOST = Array{Any}(length(Trials));
	for z = 1:length(Trials)
		# DON'T FORGET the input matrix needs to be transposed as well as the output W and H matricies
		W, H, T, cost, Vars = ParseTrials(Trials[z]);
		globalIter = size(H,3);
		nmfIter = 100000;
		numberOfProcesses = size(H,1);
		numberOfPoints   = size(inputMatrix, 1);
		numberOfSamples = size(inputMatrix, 2);
		allProcesses = zeros(numberOfPoints, numberOfProcesses, globalIter);
		allMixtures  = zeros(numberOfProcesses, numberOfSamples, globalIter);
		allDelays = zeros(size(allMixtures));
		allProcessesOG = zeros(numberOfPoints, numberOfProcesses, globalIter);
		allMixturesOG  = zeros(numberOfProcesses, numberOfSamples, globalIter);
		allDelaysOG = zeros(size(allMixtures));
		allCostOG = Array{Float64}(globalIter);
		for i =1:globalIter
			allProcessesOG[:, :, i] = H[:,:,i]';					# Again transposed
			allMixturesOG[:, :, i] = W[:,:,i]';
			allDelaysOG[:,:,i] = T[:,:,i]';
			allCostOG[i] = cost[i];
		end
		# Getting Rid of the worst solutions
		GoodTrials = ones(globalIter);
		maxVar = 0.08
		for f = 1:10
			for i=1:globalIter
				if Vars[i] > maxVar
					GoodTrials[i] = 0;
				end
				if minimum(W[:,:,i]) < .001
					GoodTrials[i] = 0;
				end
			end
			if length(findin(GoodTrials,1)) > 10
				break;
			else
				maxVar = maxVar + 0.1;
				GoodTrials = ones(globalIter);
			end
			println("The Norm cut off for $z is now $maxVar");
		end
		sizeGood = size(find(GoodTrials.==1));
		println("Good Trials after norm cut: $sizeGood");
		if length(findin(GoodTrials,1)) < 1    # So that the script doesnt stop if GoodInd=0
			GoodTrials = ones(globalIter);
		end
		sizeGood = size(find(GoodTrials.==1));
		println("Good Trials after regen cus goodtrials  = 0: $sizeGood");
		GoodInd = find(GoodTrials.== 1);
		allProcesses = allProcessesOG[:,:,GoodInd];
		allMixtures = allMixturesOG[:,:,GoodInd]
		allDelays = allDelaysOG[:,:,GoodInd];
		allCost = allCostOG[GoodInd];
		# clustering extracted processes
		idx = cluster_NMF_solutions(allProcesses, 1000);
		#idx = cluster_NMF_solutions(T[:,:,GoodInd], 1000);
		#making ordered 3D matricies
		GoodTrials = ones(size(allProcesses,3));
		orderedPros = Array{Float64}(size(allProcesses));
		orderedMix = Array{Float64}(size(allMixtures));
		orderedDel = Array{Float64}(size(allDelays));
		for k=1:size(allProcesses, 3)
			for i=1:size(allProcesses,2)
				orderedPros[:,i, k] = allProcesses[:, findin(idx[:,k],i), k];
				orderedMix[[i],:, k] = allMixtures[findin(idx[:,k],i), :, k];
				orderedDel[[i],:, k] = allDelays[findin(idx[:,k],i), :, k];
			end
		end
		# ///////////////// Clustering T ////////////////////
		orderedT = Array{Float64}(size(orderedDel,2), size(orderedDel,1), size(orderedDel,3));
		for d =1:size(orderedDel,3)
			for c =1:size(orderedDel,1)
				orderedT[:,c,d] = orderedDel[[c],:,d];
			end
		end
		maxSpM = 0.8;
		GoodTrialsOLD = ones(size(GoodTrials));
		GoodTrialsOLD[:] = GoodTrials[:];
		for jk = 1:10
			for f = 1:size(orderedT,2)
				maxClust = 6;
				Eval = Array{Float64}(maxClust-1);
				TClusters = Array{Int64}(size(orderedT,3), maxClust-1);
				TCenters = Array{Any}(maxClust-1);
				for i = 2:maxClust
					Dt = reshape(orderedT[:,f,:], size(orderedT,1), size(orderedT,3));
					R = kmeans(Dt, i; maxiter=200);
					TClusters[:,i-1] = R.assignments;
					TCenters[i-1] = R.centers;
					Dist = pairwise(CosineDist(), Dt);
					Eval[i-1] = mean(silhouettes(R, Dist));
				end
				ClustTrue = TClusters[:,findin(Eval,maximum(Eval))];
				CenterTrue = TCenters[findin(Eval,maximum(Eval))];
				StdTrue	 = Array{Float64}(size(CenterTrue[1],2));
				Tmean = Array{Float64}(size(CenterTrue[1],2));
				for i=1:size(CenterTrue[1],2)
					Tpos = abs.(CenterTrue[1][:,i]);
					Tmean[i] = mean(Tpos);
					StdTrue[i] = std(Tpos);
				end
				GoodClust = find((StdTrue./Tmean).<maxSpM)
				for k=1:size(allProcesses,3)
					if !in(ClustTrue[[k],1], GoodClust)
						GoodTrials[k] = 0;
					end
				end
			end
			if length(findin(GoodTrials,1)) > 10
				break;
			else
				maxSpM = maxSpM + 0.1;
				GoodTrials[:] = GoodTrialsOLD[:];
			end
			println("The SpM cut off for $z is now $maxSpM");
		end
		# //////////////////////////////////////////////////////////
		sizeGood = size(find(GoodTrials.==1));
		println("Good Trials after maxSpM cut: $sizeGood");
		println(length(GoodInd))
		#don't forget to get rid of terms in idx also
		GoodInd = find(GoodTrials.== 1);
		println(length(GoodInd))
		idx = idx[:,GoodInd];
		averagedH = mean(orderedPros[:,:,GoodInd],3);
		averagedW = mean(orderedMix[:,:,GoodInd],3);
		averagedT = mean(orderedDel[:,:,GoodInd],3);
		Tstd = std(orderedDel[:,:,GoodInd],3);
		averagedT = averagedT[:,:,1]';
		averagedW = averagedW[:,:,1]';
		averagedH = averagedH[:,:,1]';
		Tstd = Tstd[:,:,1]';
		allProcesses = allProcesses[:,:,GoodInd];
		allMixtures = allMixtures[:,:,GoodInd]
		allDelays = allDelays[:,:,GoodInd];
		allCost = allCost[GoodInd];
		ALLCOST[z] = allCost[:];
		# calculating stability and final processes and mixtures
		processes, mixtures, delays, avgStabilityProcesses = final_processes_and_mixtures(allProcesses, allMixtures, allDelays, idx);
		if length(avgStabilityProcesses) == 1
			avgStabilityProcesses = 1.0;
		end
		# Generate Observation Matrix by taking into account Shift
		usecomp = 1:numberOfProcesses;
		H =  processes';
		W =  mixtures';
		T=  delays';
		X=Array{Float64}(size(W,1),size(H,2));
		Hf=fft(H,2);
		Hf=Hf[:,1:Int(floor(size(Hf,2)/2))+1];
		N=size(H,2);
		f=im*2*pi*collect(0:N-1)'/N;
		#f=f[1:size(Hf,2)]'*(-1);
		if VERSION < v"0.6"
			f=f[1:size(Hf,2)]'*(-1);
		else
			f=f[1:size(Hf,2)]*(-1);
		end
		for i=1:size(W,1)
			Hft=Hf[usecomp,:].*exp.(T[i:i,usecomp]'*f);
			if mod(size(X,2),2)==0
				Hft=[Hft conj(Hft[:,end-1:-1:2])];
			else
				Hft=[Hft conj(Hft[:,end:-1:2])];
			end
			Ht=real(ifft(Hft,2));
			X[i:i,:]=W[i:i,usecomp]*Ht;
		end
		X[X.<0]=0;
		dataRecon = X'; #processes * mixtures and delays;
		dataReconCorr = zeros(numberOfSamples, 1);
		for i = 1 : numberOfSamples
			dataReconCorr[i] = cor(inputMatrix[:,i], dataRecon[:, i]);
		end
		if isdir("./Results") == false
			mkdir("./Results");
		end
		writecsv("./Results/dataReconCorr$numberOfProcesses.csv",dataReconCorr);
		writecsv("./Results/T$numberOfProcesses.csv",T);
		writecsv("./Results/W$numberOfProcesses.csv",W);
		writecsv("./Results/H$numberOfProcesses.csv",H);
		writecsv("./Results/Tstd$numberOfProcesses.csv",Tstd);
		writecsv("./Results/Reconstruction$numberOfProcesses.csv",X);
		writecsv("./Results/avgStability$numberOfProcesses.csv",avgStabilityProcesses);
		writecsv("./Results/Cost$numberOfProcesses.csv",allCost);
		VarChange[z,1] = "NMF runs for $z Sources are taken below this % of reconstruction Error: ";
		VarChange[z,2] = maxVar*100;
		VarChange[z,3] = "For $z Sources MaxSpM: ";
		VarChange[z,4] = maxSpM;
		VarChange[z,5] = "For $z Sources we have this many trials after SpM cut: ";
		VarChange[z,6] = sizeGood;
	end
	writedlm("ErrorLog.txt", VarChange);
	Plot(inputMatrix', maxSource);
end
