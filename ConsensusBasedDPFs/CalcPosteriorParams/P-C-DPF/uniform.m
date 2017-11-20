function rNum = uniform(a,b)
% Function returns a random number rNum drawn from the uniform
% disgtribution U(a,b). We assume b > a.

rNum = a + rand(1)*(b-a);

end