function [classificacao,erro]=aula34_SVMclass(coefs,svms,kernel,kpar1,kpar2,w0,X,classesX)
%essa funcao parte de parametros obtidos apos o treinamento de uma SVM e
%classifica os padroes da matriz X.
% INPUTS:
%- svms sao os vetores que dao suporte à SVM:
%   smvs(:,i) = padrao i dos dados de treinamento que possui multiplicador
%   de lagrange diferente de zero.
%- coefs sao os coeficientes correspondentes aos vetores suporte:
%   coefs(1,i)= (multiplicador do padrao i)*(classe do padrao i)
%- kernel = metodo de kernel usado na funcao aula34_SVM
%- kpar1, kpar2 = parametros dos kernels usados na funcao aula 34_SVM
%- w0 = threshold do classificador
%- X = dados a classificar (L x N)
%- classesX = classes dos dados X (caso se conheca)

if size(X,1)~=size(svms,1)
    error('Confira a dimensao dos dados')
end
if length(coefs)~=size(svms,2)
    error('Confira a dimensao dos dados')
end
if nargin==8
    if ~isequal(size(classesX),[1,size(X,2)])
        error('Confira a dimensao dos dados')
    end
else
    classesX=[];
end

%% classificacao
classificacao=zeros(1,size(X,2));
for i=1:size(X,2)
    t=sum(coefs.*aula34_CalcKernel(svms',X(:,i)',kernel,kpar1,kpar2))+w0;
    classificacao(i)=sign(t);
end

%% erro
if ~isempty(classesX)
    erro=100*sum(classesX~=classificacao)/size(X,2);
else
    erro=[];
end


%% grafico
if size(X,1)==2 %plotar em caso bidimensional
    %criando grid:
    xmin=min(X(1,:));
    xmax=max(X(1,:));
    ymin=min(X(2,:));
    ymax=max(X(2,:));
    [x,y] = meshgrid(xmin:(xmax-xmin)/50:xmax,ymin:(ymax-ymin)/50:ymax);
    z = w0*ones(size(x));
    wh = waitbar(0,'Construindo meshgrid...');
    for x1 = 1 : size(x,1)
        for y1 = 1 : size(x,2)
            ponto =[x(x1,y1) y(x1,y1)];
            for i = 1 : length(coefs)
                    z(x1,y1) = z(x1,y1) + coefs(i)*aula34_CalcKernel(ponto,svms(:,i)',kernel,kpar1,kpar2);
            end
        end        
        waitbar((x1)/size(x,1)) ;
        drawnow
    end
    delete(wh)
    drawnow
    
    figure('Color','w');
    if ~isempty(classesX)
        plot(X(1,(classesX==1)),X(2,(classesX==1)),'+r');
        hold on
        plot(X(1,(classesX==-1)),X(2,(classesX==-1)),'ob');
    else
        plot(X(1,:),X(2,:),'+r');
        hold on        
    end
    xlabel('Caract. 1')
    ylabel('Caract. 2')   
    if isempty(classesX)
       title({['SVM ',kernel,'(',num2str(kpar1),',',num2str(kpar2),')'];...
            [num2str(size(svms,2)), ' Support vectors']})
    else
       title({['SVM ',kernel,'(',num2str(kpar1),',',num2str(kpar2),')'];...
            [num2str(size(svms,2)), ' Support vectors'];...
            ['erro = ',...
            num2str(erro,'%.2f'),'%']})
    end
    contour(x,y,z,[0 0],'k')
    contour(x,y,z,[-1 -1],'b:')
    contour(x,y,z,[1 1],'r:')
    set(gca,'DataAspectRatio',[1 1 1],'XLim',...
        [min(X(1,:)) max(X(1,:))],'Ylim',[min(X(2,:)) max(X(2,:))])
end

