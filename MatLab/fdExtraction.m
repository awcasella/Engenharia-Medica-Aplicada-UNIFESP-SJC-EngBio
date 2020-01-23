function [S] = fdExtraction(sinal, Fs, f1, f2)
%FDEXTRACTION Summary of this function goes here
%   Detailed explanation goes here

[pxx, f] = pwelch(sinal', 256, 50, 2*Fs, Fs);
pxx = pxx';
f = f';
% figure;
% plot(f, pxx(:,1:5));

% 2.1 Normalized Power in the frequency bands
ini = find((f-f1) == min(abs(f - f1)));
fim = find((f-f2) == min(abs(f - f2)));

S.PER = sum(pxx(:, ini:fim),2)./sum(pxx,2);

% 2.2 Central Frequency
S.fc = sum(pxx.*repmat(f,[size(pxx, 1), 1]), 2)./sum(pxx, 2);

% 2.3 PSD on the Central Frequency
% locs = ((f - S.fc) == min(abs(f - S.fc)));
% % S.Pfc = pxx(find(min(abs(f-ca.fc))), :)';
% S.Pfc = pxx(:, locs);
% % S.Pfc = pxx(:, (f - S.fc) == min(abs(f - S.fc)));
S.Pfc = pxx(abs(f-S.fc) == min(abs(f-S.fc), [], 2));

% 2.4 Band Width Index
numerador = sum(((f - S.fc).^2).*pxx, 2);
denominador = sum(pxx,2);
S.ILB = sqrt(numerador./denominador);

% 2.5 Margin Frequency



return

