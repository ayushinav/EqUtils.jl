
"""
Narrow Band Ambiguity Function only for τ=0
The shift in frequency Δf is given by the lag in 
the cross-correlation in the frequency domain below.
Input spectra of two signals: a, b 
It returns the correlation as a function of frequency shift
Just do a correlation in the frequency domain
"""
function narrow_band_AF(a,b;glags=[size(a,1)-1, 0])
	@assert size(a,1)==size(b,1)
	@assert size(b,2)==1
	n=size(a,1)
	p_conv=Conv.P_conv(eltype(a),dsize=[n,size(a,2)], ssize=[n], 
		    gsize=[sum(glags)+1,size(a,2)],glags=glags,
		    dlags=[n-1,0], slags=[n-1,0])
	Conv.mod!(p_conv, d=a, s=b, Conv.G());
	return p_conv.g
end



"""
Wide band ambiguity using an Interpolation object
"""
function wide_band_AF(x1,x2,X; ny1=length(x1), ny2=length(x2))
	knots=(x1,x2)
	Y=Interpolations.interpolate(knots, X, Gridded(Linear()));
	y2=range(minimum(x2)+1e-5, stop=maximum(x2), length=ny2)





	αv=range(0.1, stop=2.0, length=300)
	corr=zeros(ny2,length(αv))
	y11=range(minimum(x1)+1e-5, stop=maximum(x1)*0.5, length=ny1)
	a=Y[y11,y2[50]]

	for (iα,α) in enumerate(αv)
		y12=range(minimum(x1)+1e-5, stop=maximum(x1)*0.5*α, length=ny1)
		@showprogress for i1 in 1:ny2
			b=Y[y12,y2[i1]]
			corr[i1,iα]=dot(a,b)
		end
	end
	return αv, y2, corr

end

