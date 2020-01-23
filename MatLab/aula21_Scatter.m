function [J1, J2, J3, combinacoes] = aula21_Scatter(classes, n)
% AULA21 

% The first dimension of any element of cell "classes" is the number of
% characteristcs "L". It is the same for all classes. However, the second
% element is the number of patterns in each class, it is not constant which
% means that each class has its own number of patterns.
[L, ~] = size(classes{1});

% Computing Pc = Nc/N assuming that all the classes have a different
% ammount of patterns.
N = 0;
for s = 1:length(classes)
    Nc(s) = size(classes{s},2); % Number of patterns of each class.
    N = N + size(classes{s},2); % Counting all the patterns of all classes.
    mu(:,s) = sum(classes{s},2)/Nc(s); 
end

Pc = Nc/N; % Probability of patterns in each class. (VETOR: 1xS)
Pc = Pc(:);

Q = 1;
for s = combnk(1:L, n)'
    posicoes = s';
    Sw = 0;
    dados = [];
    for t = 1:length(classes)
        mat = classes{t}(posicoes,:);
        P = Pc(t);                              % Pego sua respectiva probabilidade
        M = mu(:,t);                            % Pego sua respectiva m�dia
        Sw = Sw + P*cov(mat',1);                % Calculo a Sw como sum(P*MC)
        dados = [dados,classes{t}(posicoes,:)]; % montando o conjunto de todas as classes
        
    end
    Sm = cov(dados',1);

    J1(Q) = sum(diag(Sm))/sum(diag(Sw));
    J2(Q) = det(Sw\Sm);
    J3(Q) = sum(diag(Sw\Sm))/n;
    combinacoes(Q,:) = s(:)';
    Q = Q + 1;
end
