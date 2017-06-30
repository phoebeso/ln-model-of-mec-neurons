function y = gaussian(x,sig,c)
% Calculates gaussian of inputed x values given sigma and c values 

y = exp((-(x-c).^2)/(2*(sig^2)));

end