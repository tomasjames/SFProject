function rebinkernel,inpsf,indx,refdx

;	maxlocin = wheremax(inpsf)
	maxlocin = array_indices(inpsf,where(inpsf eq max(inpsf)))

	insize = size(inpsf,/dimen)
	refsize = (insize*indx)/refdx
	if floor(refsize[0]) mod 2 eq 0 then refsize=floor(refsize)+1 else refsize=floor(refsize)
        
        ;stop

	maxlocref = refsize/2

	xvecin = dindgen(insize[0])*indx
	xvecin = xvecin - xvecin[maxlocin[0]]
	yvecin = dindgen(insize[1])*indx
	yvecin = yvecin - yvecin[maxlocin[1]]
	xarrin = rebin(xvecin,insize[0],insize[1])
	yarrin = rebin(reform(yvecin,1,insize[1]),insize[0],insize[1])

	xvecref = dindgen(refsize[0])*refdx
	xvecref = xvecref - xvecref[maxlocref[0]]
	yvecref = dindgen(refsize[1])*refdx
	yvecref = yvecref - yvecref[maxlocref[1]]
	xarrref = rebin(xvecref,refsize[0],refsize[1])
	yarrref = rebin(reform(yvecref,1,refsize[1]),refsize[0],refsize[1])

	newpsf = myregrid(inpsf,xarrin,yarrin,xarrref,yarrref)

	; renormalize
	newpsf = newpsf/total(newpsf)

        ;stop
	return,newpsf
        
end
