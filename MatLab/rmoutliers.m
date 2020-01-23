function [dados, outliers,indexes] = rmoutliers(dados,p,g)
% Remove outliers de um conjunto de valores utilizando como limiar um numero
% p de desvios padroes em relacao a mediana.

%INPUT:
% x = vetor com dados
% p = numero de desvios padroes em relacao a mediana acima do qual as
%    amostra sao consideradas outliers (padrao = 3)
% g = plotar (1) ou nao plotar (0).

if nargin<3
    g = 0;
end
if nargin<2
    p = 3;
end

sup=find(dados>(median(dados)+ p*std(dados)));
inf=find(dados<(median(dados)- p*std(dados)));
indexes=union(sup,inf);
outliers=dados(indexes);
if g
    figure('Name','Deteccao de Outliers','Color','w');
    plot(1:length(dados),dados,'*');
    hold on
    plot(indexes,outliers,'or');
    xa=line([1 length(dados)],[median(dados), median(dados)]);
    set(xa,'Color','k')
    xa=line([1 length(dados)],[median(dados)+ p*std(dados), median(dados)+ p*std(dados)]);
    set(xa,'Color','k','LineStyle','--')
    xa=line([1 length(dados)],[median(dados)- p*std(dados), median(dados)- p*std(dados)]);
    set(xa,'Color','k','LineStyle','--')
    title(['Deteccao de Outliers:  (mediana \pm ',num2str(p),'*std)'])
    xlabel('Amostras')
    xlim([1 numel(dados)])
end

dados(indexes) = [];