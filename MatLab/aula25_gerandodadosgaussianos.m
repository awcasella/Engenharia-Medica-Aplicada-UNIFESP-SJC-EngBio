function [dadossim,classessim]=aula25_gerandodadosgaussianos(medias,covariancias, N,priors,plotar,seed)
% Essa funcao gera um conjunto de dados simulados representando um
% determinado numero de caracteristicas em um determinado numero de classes 
% As classes possuem medias distintas e covariancias distintas. Os dados
% seguem uma distribuicao gaussiana.
% INPUT:
% medias =  caracteristicas x classes (vetores com as medias das caracteristica para cada classe)
% covariancia =  caracteristicas x caracteristicas x classes (matrizes de covariancia para cada classe)
% N = numero de padroes a serem gerados
% priors = classes x 1 (prior de cada classe: probabilidade de um padrao pertencer a cada classe)
% plotar = true (exibe graficos) ou false (nao exibe graficos)
% seed = controle do seed na geracao de dados aleatorios
% OUTPUT:
% dadossim=caracteristicas x padroes: dados simulados
% classessim= vetor 1 x N contendo o numero da classe (de 1 ate C) de cada padrao

if nargin <5
    plotar=1;
end
if nargin<6
    seed=0;
end
if nargin<4
    priors=ones(size(medias,2),1)./size(medias,2);
end

M=size(medias,2);%numero de classes;
L=size(medias,1);%numero de caracteristicas;
%% TESTES DOS INPUTS:
if size(covariancias,3)~=M || size(covariancias,2)~=L || size(covariancias,1)~=L
    error('Confira a dimensao dos seus dados de input');
end
if length(priors)~=M
    error('Confira a dimensao dos priors');
end
if sum(priors)~=1
    error('Confira os valores dos priors');
end    

Ni=zeros(M,1);
for i=1:1:M %classes
    [T,P] = cholcov(squeeze(covariancias(:,:,i)),0);
    if P~=0
        error('Nao eh uma matriz de covariancia valida');
    end
    Ni(i)=round(priors(i)*N);
end

%% SIMULACAO
randn('seed', seed); 
dadossim=zeros(L,sum(Ni));
classessim=zeros(1,sum(Ni));
for i=1:1:M
        dadossim(:,sum(Ni(1:i-1))+1:sum(Ni(1:i)))=mvnrnd(medias(:,i),squeeze(covariancias(:,:,i)),Ni(i))';
        classessim(1,sum(Ni(1:i-1))+1:sum(Ni(1:i)))=i;
end

%% GRAFICOS
if plotar
if size(dadossim,1)<=3 && size(dadossim,1)>=2 && M<=6
    figure('Color','white');
    hold on
    cores={'b','r','g','m','c','y'};
    for i=1:1:M
        ps=find(classessim==i);
        try, 
            plot3(dadossim(1,ps),dadossim(2,ps),dadossim(3,ps),['.',cores{i}],'DisplayName',['Classe ',num2str(i)]);
            zlabel('Caract. 3');
            view([25 34])
        catch,
            plot(dadossim(1,ps),dadossim(2,ps),['.',cores{i}],'DisplayName',['Classe ',num2str(i)]);
        end
        xlabel('Caract. 1');
        ylabel('Caract. 2');
    end
    box on
end
end
