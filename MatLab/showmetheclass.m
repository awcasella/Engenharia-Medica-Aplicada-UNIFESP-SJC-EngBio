function [] = showmetheclass(classes)
% SHOWMETHECLASS is a function which receives a cell of matrixes, in which
% matrix represents a different class.


% L = size(classes{1},1);
for n = 1:length(classes)
    L(n) = size(classes{n}, 1);
end
if sum(L - mean(L)) ~= 0
    error(['Watch out for the dimensions of each class on the function ', mfilename]);
end

L = L(1); % How many characteristics I got
if L > 3 || L < 1
    error('Watch out! I cannot plot a multidimensional space, only a 1D, 2D or 3D')
end
cores={'b','r','k','g','m','c','y'};

figure;
for n = 1:length(classes)
    matriz = classes{n};
    switch L
        case 1
            plot(matriz, [cores{n},'.'], 'DisplayName', ['Classe ', num2str(n)])
            hold on
            xlabel('Carac 1 e 2')
        case 2
            plot(matriz(1, :), matriz(2, :), [cores{n},'.'], 'DisplayName', ['Classe ', num2str(n)])
            hold on
            xlabel('Carac 1')
            ylabel('Carac 2')
        case 3
            plot3(matriz(1, :), matriz(2, :), matriz(3, :), [cores{n},'.'], 'DisplayName', ['Classe ', num2str(n)])
            hold on
            xlabel('Carac 1')
            ylabel('Carac 2')
            zlabel('Carac 3')
    end
end

title(['Espaco com ',num2str(L),' Caracteristicas'])

box on

