function [autovalores, autovetores, y, erro] = klt(x, c, m)
% KLT 

% Input control
% if nargin ~= 2
%     error(['Verifique o n�mero de inputs da fun��o ', mfilename]);
% elseif nargin == 2
%     if ~isscalar(m)
%         error(['Verifique o input da fun��o ', mfilename]);
%     end
% end

% Removing the mean of the data in order to center the data in it.
x = x - repmat(mean(x,2), [1, size(x,2)]);

Rx = (x*x')/size(x, 2); % Tanto faz usar o cov ou x*x' desde que a meadia seja zero.
% Rx = cov(x'); % Tanto faz usar o cov ou x*x' desde que a meadia seja zero.
[eigvec, eigval] = eig(Rx);
Ry = eigvec'*Rx*eigvec;

classe1 = x(:, c==c(1));
classe2 = x(:, c==c(end));
for n = 1:size(x, 1)
    resultadoAUC(n) = rocmebabe(classe1(n, :), classe2(n, :));
end
variancias = var(x');
if resultadoAUC(variancias == max(variancias)) ~= 1
    disp('AUC diferente de 1 para a caracteristica com maior variancia, PCA pode não ser a melhor opcao.');
end

% Error and sorting "Verificar se est� certa a ordena��o"
[eigval, eigvec] = ordenar(eigval, eigvec);

% Defining two of the outputs
autovalores = diag(eigval);
autovetores = eigvec;

% Calculating the sum of eigenvalues
eigval = diag(eigval);
total = sum(eigval);

% Selecting only the principal components
A = eigvec(:, 1:m);

% Rotating the data
y = A'*x;

% Calculating the error associated with dimension reduction.
erro = sum(eigval(m+1:end))/total;

return

% Function to order the eigenvalues and eigenvectors
function [eigvalord, eigvecord] = ordenar(eigval, eigvec)

eigvalord = fliplr(sort(diag(eigval))');
eigvecord = zeros(size(eigvec));
for n = 1:length(eigvalord)
    locs = find(eigvalord(n) == diag(eigval), 1);
    eigvecord(:, n) = eigvec(:, locs);
end
eigvalord = diag(eigvalord);

return

