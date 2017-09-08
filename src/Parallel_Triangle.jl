function PTri(pnum, T, Tstd, W, micPos, xtol)
	positions, speed, init, minf = Triangulate(T, Tstd, W, micPos, xtol)
end

function Parallel_Tri(pnum, T, Tstd, W, micPos, xtol)
	cores=nprocs()
	sendto(1:cores, T=T)
	sendto(1:cores, Tstd=Tstd)
	sendto(1:cores, W=W)
	sendto(1:cores, micPos=micPos)
	sendto(1:cores, xtol=xtol)
	Try = pmap(pnum->PTri(pnum, T, Tstd, W, micPos, xtol), 1:pnum)
	return Try
end