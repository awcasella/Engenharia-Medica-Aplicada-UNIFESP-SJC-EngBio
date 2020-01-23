function [S, DEPar, f] = featureExtraction(sinal, Fs, bandas)
%CEXTRACTION Summary of this function goes here
%   Detailed explanation goes here


[M,N] = size(sinal);

% 1.1 Average 
S.mu0 = sum(sinal, 2)/N;

% 1.2 Rectified Average
S.muR = sum(abs(sinal), 2)/N;

% 1.3 Variance
media = repmat(S.mu0, [1, N]);
S.var0 = sum((sinal - media).^2, 2)/(N-1);
% ca.var0 = var(sinal, 0, 2);

% 1.4 Mean of the first order diferences
S.mu1 = sum(diff(sinal,1,2), 2)/(N-1);

% 1.5 Variance of the first order derivative
media = repmat(S.mu1, [1, N-1]);
S.var1 = sum((diff(sinal,1,2)-media).^2, 2)/(N-2);
% ca.var1 = var(diff(sinal, 1, 2), 0, 2);

% 1.6 Mean of second order differences
S.mu2 = sum(diff(sinal, 2, 2), 2)/(N-2);

% 1.7 Variance of second order differences
media = repmat(S.mu2, [1, N-2]);
S.var2 = sum((diff(sinal, 2, 2) - media).^2, 2)/(N-3);
% ca.var2 = var(diff(sinal, 2, 2), 0, 2);

% 1.8 Statistical Mobility
S.M = (S.var1)./(S.var0);

% 1.9 Statistical Complexity
S.C = sqrt((S.var2)./(S.var1) - S.M);

% Aqui estimamos a densidade de potencia espectral (PSD).
% Como estamos trabalhando com variabilidade da frequencia cardiaca (HRV),
% precisamos de uma resolução muito alta no dominio da frequencia pois as
% informações pertinentes estão em bandas de frequencia muito estreitas. 
% 
% Assim, optou-se por utilizar um metodo autoregressivo para estimar o PSD.
% A ordem utilizada para o modelo foi 11, com uma frequencia de amostragem
% (interpolada no script) de 4 Hz.
ordem = 11;
% [DEPar, f] = pburg(sinal', ordem, 2^ordem, Fs);

[DEPar, f] = pwelch(sinal', 256, 50, [],Fs);
% [pxx, f] = pwelch(y2, 256, [], 256, Fi);
DEPar = DEPar';
f = f';

% 2.1 Normalized Power in the frequency bands
f1 = bandas(:, 1);
f2 = bandas(:, 2);

for z = 1:size(bandas, 1)
    ini = find((f-f1(z)) == min(abs(f - f1(z))));
    fim = find((f-f2(z)) == min(abs(f - f2(z))));
    eval(['S.PER', num2str(z),' = sum(DEPar(:, ini:fim),2)./sum(DEPar,2);']);
end


% 2.2 Central Frequency
S.fc = sum(DEPar.*repmat(f,[size(DEPar, 1), 1]), 2)./sum(DEPar, 2);

% 2.3 PSD on the Central Frequency
S.Pfc = DEPar(abs(f-S.fc) == min(abs(f-S.fc), [], 2));

% 2.4 Band Width Index
numerador = sum(((f - S.fc).^2).*DEPar, 2);
denominador = sum(DEPar,2);
S.ILB = sqrt(numerador./denominador);

% 2.5 Edge Frequency

soma = cumsum(DEPar, 2)./repmat(sum(DEPar,2),[1, size(DEPar,2)]);
for n = 1:size(DEPar,1)
    aux = abs(soma(n,:)-0.9);
    locs = find(aux == min(aux));
    vet(n) = f(locs);
end
S.EdgeFreq = vet';

return

