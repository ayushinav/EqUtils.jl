
"""
th is in degrees
"""
function vecavg(th)
	mean(broadcast(x->complex(sin(deg2rad(x)), cos(deg2rad(x))), th))
end

"""
"""
function circ_mean(th)
	m=vecavg(th)
	s=real(m); c=imag(m)
	if(s>0 && c>0)
		return rad2deg(atan(s*inv(c)))
	elseif(c < 0)
		return rad2deg(atan(s*inv(c))) + 180
	elseif(s<0 && c>0)
		return rad2deg(atan(s*inv(c))) + 360
	else
		return 0
	end
end



