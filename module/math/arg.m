function arg = arg(z)
% return value belongs [ 0, 2pi ]
arg = angle(z); 
if arg < 0 % since angle returns value between -pi and pi
    arg = arg + 2 * pi;
end