module SDEPlot
using Winston
export visualize_posterior

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


end