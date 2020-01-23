function [varargout] = rocmebabe(classe1, classe2, varargin)
% ROCMEBABE is a function which receives two arrays corresponding to two
% different classes of a characteristics array, the it calculates the ROC
% (Receiver Operating Characteristic) curve, plot it and return the area 
% under the curve (AUC) by the Riemman sum.

% Input control
if nargin > 3 || nargin == 1
    error('Verifique o n�mero de inputs');
elseif nargin == 2
    plotar = false;
elseif nargin == 3
    plotar = varargin{1};
end

% Input adustment
if classe1(end) > classe2(end)
    aux = classe1;
    classe1 = classe2;
    classe2 = aux;
end

% Let's do the work
ini = 0;
fim = 1;
classe1 = classe1(:);
classe2 = classe2(:);
distribuicoes = [union(classe1, classe2); 0];

VP = zeros(size(distribuicoes));
FP = zeros(size(distribuicoes));

for n = length(distribuicoes):-1:1
    s = length(distribuicoes) - n + 2;
    if any(distribuicoes(n) == classe1) && ~any(distribuicoes(n) == classe2)
        VP(s) = VP(s-1);
        FP(s) = FP(s-1) + 1;
    elseif ~any(distribuicoes(n) == classe1) && any(distribuicoes(n) == classe2)
        VP(s) = VP(s-1) + 1;
        FP(s) = FP(s-1);
    elseif any(distribuicoes(n) == classe1) && any(distribuicoes(n) == classe2)
        VP(s) = VP(s-1) + 1;
        FP(s) = FP(s-1) + 1;
    end
end

% Normalized data
VP = VP/max(VP);
FP = FP/max(FP);

% Riemann Sum
AUC = sum(((VP(2:end))+(VP(1:end-1)))/2.*diff(FP));

if nargout == 0 || plotar
    figure;
    plot(FP, VP, 'Color', 'r', 'LineWidth', 3);
    hold on;
    line([ini, fim], [ini, fim], 'color', 'k', 'LineWidth', 3);
    gordurinha = ((fim-ini)/20);
    xlim([ini - gordurinha, fim + gordurinha]);
    ylim([-0.05, 1.05]);
    grid on
%     title('ROC curve, bitch!');
    ylabel('%VP(1-\beta)');
    xlabel('%FP(\alpha)');
end
if nargout == 1
    varargout{1} = AUC;
end