function [T_BI] = mrp2trafo(sigma)

T_BI = eye(3) + (8*skew(sigma)^2 - 4*(1-sigma'*sigma)*skew(sigma)) / ((1+sigma'*sigma)^2);

end