function [ordem, M] = SelecaoEscalar(Mcorr, criterios, varargin)
% SELECAOESCALAR

if length(criterios) ~= size(Mcorr,1)
    error('Atencao! Verifique se os criterios correspondem a matris de correlacao.')
end

if nargin == 1 || nargin > 4
    error('Atencao! Passe o numero certo de inputs.');
elseif nargin == 2
    a1 = 0.5; % Default
    a2 = 0.5; % Default
elseif nargin == 3
    a1 = varargin{1};
    a2 = 1 - a1;
elseif nargin == 4
    a1 = varargin{1};
    a2 = varargin{2};
end

% Making it a column vector
criterios = criterios(:);

% The first one is special, let's calculate it outside the loop
M(1) = max(criterios);
ordem(1) = find(criterios == M(1));
Mcorr = abs(Mcorr);
% let's do the work
for k = 2:length(criterios)
    % getting the correlation of the previous characteristc with all the
    % others
    fator = Mcorr(:, ordem); 
    s = length(ordem);
    % Calculating all the possible correlation factors.
    fator = (sum(fator, 2)/(s));
    % Let's test the new possible criterias for all the possible
    % correlation factors
    teste = a1*criterios - a2*fator;
    % Let's put them all together...
    matriz = [1:length(criterios); teste'];
    % ... and delete the ones we've got already.
    matriz(:,ordem) = [];
    % Here we find the value of the maximum new "Mk"
    M(k) = max(matriz(2,:));
    % and here we find where he's hinding.
    locs = matriz(2,:) == M(k);
    nova = matriz(1,:);
    % and we place the next best characteristc in the vector "ordem".
    ordem(k) = nova(locs);
    
    % and here we go again in the loop
end



