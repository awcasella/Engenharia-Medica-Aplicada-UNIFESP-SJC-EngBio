function [VETORNORM] = normalizacao(VETOR, metodo, r)
% Normaliza os valores de uma determinada caracteristica.
%
% INPUT:
% - VETOR: valores de uma caracteristica calculada em N padroes e C classes,
%   Vetor com dimensao 1 x N*C
% - metodo ='std' : normalizacao linear (padrao)
%       = 'mmx': limitada entre -1 e 1
%       = 'sfm': rescala nao linear no intervalo 0 a 1
% - r = parametro do metodo sfm (padrao =1)
%
% OUTPUT:
% - VETORNORM = 1 x N*C: vetor com os valores normalizados da caracteristica
% 

if nargin<2
    metodo='std';
end

if nargin<3
    r=1;
end
[M, N] = size(VETOR);
VETOR = VETOR(:);
switch metodo
    case 'std'
        VETORNORM=VETOR-mean(VETOR);
        VETORNORM=VETORNORM./(std(VETOR));
    case 'mmx'
        VETORNORM=2*VETOR./(max(VETOR)-min(VETOR));
        VETORNORM=VETORNORM-(min(min(VETORNORM))+1);       
    case 'sfm'
        Y=VETOR-mean(VETOR);
        Y=Y./(r*std(VETOR));
        VETORNORM=1./(1+exp(-Y));
    otherwise
        errordlg('Método desconhecido');
        VETORNORM=[];      
end

VETORNORM = reshape(VETORNORM, [M, N]);