function tau = alpha2tau(alpha, dt)
%ALPHA2TAU Converts the infection rate to the average time of an infection.

tau = dt ./ alpha;

end

