function Pexecute(p, numiter, Obs, noc, opts)
	X=Obs
	noc=noc
	opts=opts
	W, H, T, varexpl, cost = execute(X, noc, opts)
end

function sendto(p; args...)
	for i in p
		for (nm, val) in args
			@spawnat(i, eval(Main, Expr(:(=), nm, val)))
		end
	end
end

function Parallel_execute(pnum, iter, Obs, num_sources, opts)
	cores=nprocs()
	sendto(1:cores, iterations=iter) #TODO i do not think sendto is needed
	sendto(1:cores, Obs=Obs)
	sendto(1:cores, noc=num_sources)
	sendto(1:cores, opts=opts)
	Try = pmap(a->Pexecute(a, iter, Obs, num_sources, opts), 1:pnum)
	return Try
end