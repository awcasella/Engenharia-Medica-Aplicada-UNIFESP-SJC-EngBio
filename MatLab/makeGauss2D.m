function [P, T] = makeGauss2D(x, y, medias, covariancia, varargin)
%MAKEGAUSS2D Summary of this function goes here
%   Detailed explanation goes here

condicao = false;
if nargin == 6
    f = varargin{1};
    tipo = varargin{2};    
    if isequal(lower(f),'plot')
        condicao = true;
    end
end
% x = x(:);
% y = y(:);

n = length(medias);
norm = 1/(sqrt((2*pi)^n)*(sqrt(det(covariancia))));

% P = zeros(length(x), length(y));

alpha = x - medias(1);
betha = y - medias(2);
T = [alpha; betha];

P = norm*exp(-0.5*(T')*(inv(covariancia))*(T));
% z = [x(c); y(d)] - medias;
% P(c,d) = norm*exp(-0.5*(z')*(inv(covariancia))*(z));

if condicao
    try
        figure;
        eval([tipo,'(x,y,P);']);
    catch
        error('Verifique o tipo de plotagem que deseja');
    end
end


return