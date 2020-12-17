function nr_of_deaths = total_deaths(M,l)
%TOTAL_DEATHS 
% This function takes two variables as its input, M and l
% M is an instance of the ModelResult class
% l is a vector with lethality percentage after contracting covid, per age
% group

nr_of_deaths = dot(l , sum(M.StoI * M.DeltaT, 2) );

end

