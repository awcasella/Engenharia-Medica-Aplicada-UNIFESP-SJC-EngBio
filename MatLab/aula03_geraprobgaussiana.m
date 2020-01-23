function [P]=aula03_geraprobgaussiana(x,y,medias,covariancia)
% Essa fun��o calcula a probabilidade gaussiana bivariada em um 
% determinado n�mero de pontos do plano x-y
% INPUT:
% x = vetor com pontos no eixo x
% y = vetor com pontos no eixo y
% medias = 2 x 1 (m�dias de cada vari�vel da distribui��o)
% covariancia = 2 x 2 (matriz de covariancia das vari�veis)
% OUTPUT:
% P= probabilidade para cada valor de (x,y): P=P(x,y|medias,covariancia)

% if size(medias,1)~=2 || size(covariancia,1)~=2 || size(covariancia,2)~=2
%     error('Confira a dimens�o dos seus dados');
% end
L=size(medias,1);
Norm = 1/(sqrt(((2*pi)^L)*det(covariancia))); %normaliza��o
P=zeros(length(x),length(y));
for i=1:1:length(x)
    for j=1:1:length(y)
        z=[x(i) ; y(j)]-medias;
        P(i,j)=Norm*exp(-0.5*z'*inv(covariancia)*z);
    end
end


