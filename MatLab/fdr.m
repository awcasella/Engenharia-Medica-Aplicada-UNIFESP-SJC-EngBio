function [result] = fdr(classe1, classe2)
%FDR Summary of this function goes here
%   Detailed explanation goes here

result = ((mean(classe1)-mean(classe2))^2)/(var(classe1)+var(classe2));

return

