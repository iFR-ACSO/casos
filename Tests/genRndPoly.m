% ========================================================================
%
%
% Generate Polynomials for testing purposes.
%
% Description: Generate n different polynomials that differ in min/max
% degree, no. of indeterminates
%
% Procedure:   In several for-loop different polynomials are generated
% using the multipoly-toolbox. The polynomials and their information (e.g.
% max degree) are stored in a.mat file. This .mat file is then read in in
% test scripts.
%                
% Date: 02/06/2024 
%
%
% ========================================================================



clear
clc

% we on
m = 1;

testPolyStruct = struct();

for n = 1:10    % generate different kinds of indetermininates 

          if mod(n,2) && n < 8
              
                x = mpvar('x',n,m);

                indet = x;

          elseif n >= 8

                x = mpvar('x',n,m);
                t = mpvar('t',1,1);
                w = mpvar('w',2,m);

                indet = [x;t;w];
          else
                x = mpvar('x',n,m);
                t = mpvar('t',1,1);
        
                indet = [x;t];
          end
        
         for minDeg = 0:3
             for maxDeg = minDeg:6

                 mono = monomials(indet,minDeg:maxDeg);

                 % here comes the randomnes
                 poly = rand(length(mono),1)'*mono;

             end
         end

end