function setdir(dir)
	if isdir(dir)
		cd(dir)
	end
end

function setdir()
	dir = remotecall_fetch(()->pwd(), 1)
	setdir(dir)
end