function alpha = tau2alpha(tau, dt)
%TAU2ALPHA Converts the average time of an infection to the recovery rate.

alpha = dt ./ tau;

end

