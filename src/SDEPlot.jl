module SDEPlot
using SDE
using Winston
export visualize_posterior


function graph(X::MvPath)
    p = FramedPlot()
    if ndims(X.yy) == 1 || (ndims(X.yy) == 2 && size(X.yy,1) == 1)
        setattr(p, "xrange", SDE.hrange(X.tt))
        setattr(p, "yrange", SDE.hrange(X.yy[:]))
        add(p, Curve(X.tt,X.yy[:], "color","black", "linewidth", 0.5))
    
    elseif ndims(X.yy) == 2 && size(X.yy,1) == 2
        setattr(p, "xrange", SDE.hrange(X.yy[1,:]))
        setattr(p, "yrange", SDE.hrange(X.yy[2,:]))
        add(p, Curve(X.yy[1,:],X.yy[2,:], "color","black", "linewidth", 0.5))
    
    end
    p
end


#%  .. function:: visualize_posterior(post[, truedrift])
#%  	
#%  	Plot 2r*se wide marginal credibility bands, where ``post`` is the result of 
#%  	bayes_drift and truedrift the true drift (if known :-) ).
#% 


function visualize_posterior(post, truedrift, r)
	x = post[:,1] 
	bhat = post[:,2] 
	c = post[:,3].*r #standard errors method 1
	xfull = linspace(x[1], x[end], 1000)
	b = map(truedrift,xfull)



	p = Winston.FramedPlot()
	Winston.setattr(p, "aspect_ratio", 1.)
	Winston.setattr(p.frame, "tickdir", +1)
	Winston.setattr(p.frame, "draw_spine", false)
	Winston.add( p, Winston.FillBetween(x, bhat+c, x, bhat-c, "color", "grey") )
	Winston.add( p, Winston.Curve(x, bhat, "type", "solid", "color", "blue") )
	Winston.add( p, Winston.Curve(xfull, b, "type", "dashed", "color", "red") )

	return(p)
end
function visualize_posterior(post, r)
	x = post[:,1] 
	bhat = post[:,2] 
	c = post[:,3] .* r
	 

	p = Winston.FramedPlot()
	Winston.setattr(p, "aspect_ratio", 1.)
	Winston.setattr(p.frame, "tickdir", +1)
	Winston.setattr(p.frame, "draw_spine", false)
	Winston.add( p, Winston.FillBetween(x, bhat+c, x, bhat-c, "color", "grey") )
	Winston.add( p, Winston.Curve(x, bhat, "type", "solid", "color", "blue") )
 
	return(p)
end

#plot (2d) sample path

function pl(x::Array{ FloatingPoint,2})
	p = FramedPlot()
	#compute range
	R1 = range(x[1,:]) 
	R2 = range(x[2,:])

	#widen a bit	
	R1 = (R1[1] -0.1*(R1[2]-R1[1]), R1[2] + 0.1*(R1[2]-R1[1]))
	R2 = (R2[1] -0.1*(R2[2]-R2[1]), R2[2] + 0.1*(R2[2]-R2[1]))

	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

	add(p, Curve(x[1,:],x[2,:]))
	Winston.display(p)

	p	
end


end