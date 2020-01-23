function [P]=aula03_geraprobgaussiana(x,y,medias,covariancia)
% Essa função calcula a probabilidade gaussiana bivariada em um 
% determinado número de pontos do plano x-y
% INPUT:
% x = vetor com pontos no eixo x
% y = vetor com pontos no eixo y
% medias = 2 x 1 (médias de cada variável da distribuição)
% covariancia = 2 x 2 (matriz de covariancia das variáveis)
% OUTPUT:
% P= probabilidade para cada valor de (x,y): P=P(x,y|medias,covariancia)

% if size(medias,1)~=2 || size(covariancia,1)~=2 || size(covariancia,2)~=2
%     error('Confira a dimensão dos seus dados');
% end
L=size(medias,1);
Norm = 1/(sqrt(((2*pi)^L)*det(covariancia))); %normalização
P=zeros(length(x),length(y));
for i=1:1:length(x)
    for j=1:1:length(y)
        z=[x(i) ; y(j)]-medias;
        P(i,j)=Norm*exp(-0.5*z'*inv(covariancia)*z);
    end
end


