% INTRODU��O AOS PROCESSOS ESTOC�STICOS E MODELOS LINEARES
% Aulas 09 e 11: Engenharia M�dica Aplicada
% Prof. Dr. Adenauer G. Casali

clear all 

%% PROCESSO ESTOC�STICOS
% Um determinado sinal biol�gico {x_i}, i=1...n, �
% uma s�rie temporal discreta que decorre de algum processo que o gera,
% isto �, de alguma REGRA que determina os valores da s�rie ao longo do 
% tempo. Estas series, na maioria dos casos de interesse, n�o s�o sempre as
% mesmas, elas variam, mesmo se o PROCESSO que a gera for sempre o mesmo 
% (por exemplo, trechos do ECG gerados pelo sistema cardiovascular de um 
% determinado individuo saud�vel). Quando isto acontece dizemos que o 
% processo em quest�o � ESTOC�STICO.
%
% Podemos pensar em um processo estocastico como uma populacao estatistica
% formada por todas as possiveis series temporais que s�o geradas pela 
% regra. Denotaremos esta popula��o pela vari�vel X mai�scula. 
% Considere por exemplo o processo que 
% chamaremos de AR1 e que definiremos por: 
% X_(i+1) = A*X_(i) + alpha*R, 
% onde A e alpha s�o constantes e R um n�mero aleat�rio retirado de uma 
% distribui��o normal com vari�ncia 1 e m�dia zero.
% Este processo, entendido como a fam�lia de todas as poss�veis s�ries 
% com esta forma, pode ser AMOSTRADO: uma AMOSTRA do processo � uma 
% determinada s�rie {x_i}, i = 1...n resultante deste processo. Por exemplo, 
% encontremos uma amostra do processo AR1 com A = 0.8 e alpha = 1:

A=0.8;
alpha=1;
x1(1)=0;
for i=2:1:1000
    x1(i)=A*x1(i-1)+alpha*randn;
end
% O vetor x1 no seu workspace � uma amostra, com 1000 pontos, do
% processo AR1:
h=figure;
plot(x1);
xlabel('Pontos')
ylabel('x1');
title('Amostra do processo AR1')

% Podemos gerar outra s�rie desse tipo (outra amostra):
delete(h)
x2(1)=alpha*randn;
for i=2:1:1000
    x2(i)=A*x2(i-1)+alpha*randn;
end
clear i
h=figure;

subplot(2,1,1)
plot(x1);
xlabel('Pontos')
ylabel('x1');
ylim([-6 6])
title('Amostra do processo AR1')

subplot(2,1,2)
plot(x2);
xlabel('Pontos')
ylabel('x2');
ylim([-6 6])
title('Outra amostra do processo AR1')

% Cada vez que amostramos o processo AR1,encontraremos um vetor x 
% diferente (uma amostra diferente).Podemos ent�o falar da 
% probabilidades de se encontrar determinados valores em determinados
% pontos da nossa amostra. Vamos amostrar o processo AR1 500 vezes e
% armazenar cada uma dessas amostras na primeira dimens�o de uma vari�vel
% x. A segunda dimens�o dessa vari�vel � reservada para os pontos da amostra 
% (no tempo, se vc preferir pensar assim):

clear x1 x2 x
try, delete(h), catch end
x(:,1)=alpha*randn(500,1); % inicializa��o (primeiro ponto de cada amostra
for i=2:1:1000 % amostragem: 1000 amostras diferentes
    x(:,i)=A*x(:,i-1)+alpha*randn(500,1);
end
clear i h

% Vamos inspecionar as primeiras quatro amostras que obtivemos:
figure('Color','w');
for i=1:1:4
    subplot(4,1,i)
    plot(x(i,:));
    xlabel('Pontos')
    ylabel(['x_',num2str(i),'(t)']);
    title(['Amostra ',num2str(i),' do processo AR1'])
    ylim([-6 6])
end
clear i

% Ja que temos 500 diferentes amostras, cada uma delas com 1000 pontos, nos 
% podemos plotar histogramas com as probabilidades de obtermos
% determinados valores em cada ponto (no tempo) das amostras. Por exemplo,
% esta � a forma da distribuicao de probabilidade das amostras no ponto 15:

% figure('Color','w')
% cont = 1;
% for ponto= [15,17, 915];
for ponto = 15;
figure('Color','w')
[N,BIN]=histc(x(:,ponto),[min(x(:,ponto))-0.25:0.5:max(x(:,ponto))+0.25]);
% subplot(3,1,cont)
bar([min(x(:,ponto))-0.25:0.5:max(x(:,ponto))+0.25],100*N./sum(N),1);
xlim([-6,6])
xlabel('Valor de x')
ylabel('Probabilidade de encontrar x (%)')
title({['Distribui��o de probabilidade do ponto ',num2str(ponto),' do processo AR1 '];'(entre todas 500 amostras)'})
% mean(BIN)
% cont = cont + 1;
end

% e esta � a distribui��o de probabilidade no ponto 315:
ponto = 315;
figure('Color','w')
[N,BIN]=histc(x(:,ponto),[min(x(:,ponto))-0.25:0.5:max(x(:,ponto))+0.25]);
bar([min(x(:,ponto))-0.25:0.5:max(x(:,ponto))+0.25],100*N./sum(N),1);
xlabel('Valor de x')
ylabel('Probabilidade de encontrar x (%)')
title({['Distribui��o de probabilidade do ponto ',num2str(ponto),' do processo AR1 '];'(entre todas 500 amostras)'})
clear ponto BIN N

% Voc� deve ter observado que apesar das distribui��es terem muitas 
% semelhan�as, elas n�o s�o id�nticas. 

% Observe ainda que n�s podemos calcular distribui��es conjuntas. 
% Por exemplo, podemos perguntar por qual seria a probabilidade
% de encontrarmos x(:,15) em determinado intervalo e x(:,16) em outro 
% intervalo, isto �, a probabilidade conjunta das amostras nos pontos 15 e
% 16.
% Como uma ilustra��o desse processo, considere bins de largura = 1 e o 
% seguinte c�digo, escrito de modo pouco eficiente mas did�tico:
close all
bins=[-5:1:5];
ponto1=15;
ponto2=16;
COUNT=zeros(length(bins)+1,length(bins)+1); % contar quantos eventos est�o 
% em cada intervalo definido pela vari�vel bins
for i=1:1:length(bins)+1 %dividindo as amostras no primeiro ponto
    if i==1
        a1=find(x(:,ponto1)<bins(i)); %amostras com x menor que o limite inferior
    elseif i==length(bins)+1
        a1=find(x(:,ponto1)>=bins(i-1)); %amostras com x maior que o limite superior
    else
        a1=find(x(:,ponto1)<bins(i) & x(:,ponto1)>=bins(i-1)); %amostras entre dois limites
    end
    for j=1:1:length(bins) %dividindo as amostras no segundo ponto
        if j==1
            a2=find(x(:,ponto2)<bins(j));
        elseif j==length(bins)+1
            a2=find(x(:,ponto2)>=bins(j-1));
        else
            a2=find(x(:,ponto2)<bins(j) & x(:,ponto2)>=bins(j-1));
        end
        COUNT(i,j)=length(intersect(a1,a2)); %contando quantos valores estão nos limites especificados para cada ponto
    end
end
PROB=100*COUNT/sum(sum(COUNT)); %probabilidade conjunta entre x(:,15) e x(:,16)
clear COUNT a1 a2 i j

% A vari�vel PROB(1,1) � a probabilidade de encontrar um valor menor que -5
% no ponto 15 e um valor menor que -5 no ponto 16.
% PROB(1,2) � a probabilidade de encontrar um valor menor que -5 no ponto 15 
% e um valor entre -5 e -4 no ponto 16 (e sucessivamente para o intervalo 
% entre -5 e 5):
figure;
[E1,E2]=meshgrid([-5.5:1:5.5]);
surf(E1,E2,PROB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
clear E1 E2
axis tight
xlabel(['Valor de x no ponto ',num2str(ponto1)])
ylabel(['Valor de x no ponto ',num2str(ponto2)])
zlabel('Probabilidade (%)')
title(['Probabilidade Conjunta de x(:,',num2str(ponto1),') e x(:,',num2str(ponto2),')'])
view([-68 44])
rotate3d on

% Inspecione esta probabilidade, e note como h� uma importante rela��o 
% linear entre os dois pontos. Compare essa probabilidade conjunta com o 
% produto entre as probabilidades destes dois pontos:
close all
figure;
subplot(1,2,1)
[E1,E2]=meshgrid([-5.5:1:5.5]);
surf(E1,E2,PROB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
clear E1 E2
axis tight
xlabel(['Valor de x no ponto ',num2str(ponto1)])
ylabel(['Valor de x no ponto ',num2str(ponto2)])
zlabel('Probabilidade (%)')
a1=zlim;
title(['Probabilidade Conjunta de x(:,',num2str(ponto1),') e x(:,',num2str(ponto2),')'])
view([-68 44])
PRODPROB=100*repmat(sum(PROB,1)/100,size(PROB,1),1).*repmat(sum(PROB,2)/100,1,size(PROB,2));
subplot(1,2,2)
[E1,E2]=meshgrid([-5.5:1:5.5]);
surf(E1,E2,PRODPROB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
clear E1 E2
axis tight
xlabel(['Valor de x no ponto ',num2str(ponto1)])
ylabel(['Valor de x no ponto ',num2str(ponto2)])
zlabel('Probabilidade (%)')
zlim(a1)
title(['Produto das Probabilidades de x(:,',num2str(ponto1),') e x(:,',num2str(ponto2),')'])
view([-68 44])
rotate3d on

% A clara diferen�a entre essas duas distribui��es � um sinal de que as
% amostras 15 e 16 N�O s�o independentes. De fato, h� uma rela��o linear
% entre elas (conhecemos o processo!). 
clear PROB PRODPROB a1 bins ponto1 ponto2

% Vamos fazer o mesmo plot agora para as amostras 15 e 915:

bins=[-5:1:5];
ponto1=15;
ponto2=915;
COUNT=zeros(length(bins)+1,length(bins)+1);
for i=1:1:length(bins)+1
    if i==1
        a1=find(x(:,ponto1)<bins(i));
    elseif i==length(bins)+1
        a1=find(x(:,ponto1)>=bins(i-1));
    else
        a1=find(x(:,ponto1)<bins(i) & x(:,ponto1)>=bins(i-1));
    end
    for j=1:1:length(bins)
        if j==1
            a2=find(x(:,ponto2)<bins(j));
        elseif j==length(bins)+1
            a2=find(x(:,ponto2)>=bins(j-1));
        else
            a2=find(x(:,ponto2)<bins(j) & x(:,ponto2)>=bins(j-1));
        end
        COUNT(i,j)=length(intersect(a1,a2));
    end
end
PROB=100*COUNT/sum(sum(COUNT)); %probabilidade conjunta entree x(:,15) e x(:,915)
clear COUNT a1 a2 i j
figure;
subplot(1,2,1)
[E1,E2]=meshgrid([-5.5:1:5.5]);
surf(E1,E2,PROB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
clear E1 E2
axis tight
xlabel(['Valor de x no ponto ',num2str(ponto1)])
ylabel(['Valor de x no ponto ',num2str(ponto2)])
zlabel('Probabilidade (%)')
a1=zlim;
title(['Probabilidade Conjunta de x(:,',num2str(ponto1),') e x(:,',num2str(ponto2),')'])
view([-68 44])
subplot(1,2,2)
PRODPROB=100*repmat(sum(PROB,1)/100,size(PROB,1),1).*repmat(sum(PROB,2)/100,1,size(PROB,2));
[E1,E2]=meshgrid([-5.5:1:5.5]);
surf(E1,E2,PRODPROB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
clear E1 E2
axis tight
xlabel(['Valor de x no ponto ',num2str(ponto1)])
ylabel(['Valor de x no ponto ',num2str(ponto2)])
zlabel('Probabilidade (%)')
zlim(a1)
title(['Produto das Probabilidades de x(:,',num2str(ponto1),') e x(:,',num2str(ponto2),')'])
view([-68 44])
rotate3d on

% Voc� certamente deve ter notado que boa parte dessa correla��o foi perdida
% e que as distribui��es de probabilidade agora s�o semelhantes. Isso se d�
% porque na medida em que o tempo passa, o determinismo da rela��o entre
% amostras no tempo para nosso modelo AR1 vai sendo perdido pela influ�ncia 
% da aleatoriedade que afeta cada um dos pontos independentemente. 
% As amostras passam a se tornar INDEPENDENTES. 
% 
% Estes s�o os crit�rios de independ�ncia de dois 
% pontos (tempo) em um processo estoc�stico:
% -> X_1 e X_2 s�o independentes no sentido FORTE se
%     P(X_1,X_2) = P(X_1)P(X_2). 
% -> X_1 e X_2 s�o independentes no sentido FRACO se
%     E[X_1, X_2] = E[X_1]E[X_2].  (onde E aqui denota o valor esperado).

% Quando todos os pontos temporais de um processo forem independentes entre
% si (no sentido forte ou no sentido fraco) o processo � dito aleat�rio (no
% sentido forte ou no sentido fraco, respectivamente). O processo AR1 que
% definimos n�o � aleat�rio. Mas, claramente, tamb�m n�o � determin�stico.

clear PROB PRODPROB ponto1 ponto2 bins a1

%% ESTACIONARIEDADE DE PROCESSOS ESTOC�STICOS
% -> Um processo estoc�stico � dito estacion�rio no sentido forte se as
% distribui��es de probabilidade s�o INVARIANTES NO TEMPO:
%  
%    P(x_1, x_2, ....) = P(x_{1+k}, x_{2+k}, ...) PARA TODO k
% 
% -> Um processo � dito estacion�rio no sentido fraco se os valores
% esperados dos dois primeiros momentos estat�sticos (pelo menos) s�o
% invariantes no tempo:
% 
%  E[x(i)] = E[x(i+k)] para todo i e para todo k.
%  E[x(i) x(i+j)] = E[x(i+k) x(i+j+k)] para todo i, j e k.
% 
% Ou seja, a ideia de estacionariedade � a seguinte: o tempo (o ponto da
% amostra) n�o tem NADA de especial do ponto de vista estat�stico. 
% 
% Vamos explorar a estacionariedade no sentido fraco (suficiente para o
% nosso estudo) das amostras armazenadas na vari�vel x. Primeiro, vamos
% calcular o primeiro momento estat�stico entre todas as amostras para cada
% ponto:

Em1 = mean(x,1); % note, essa m�dia n�o � no tempo, mas entre as diferentes amostras
figure('Color','w');

subplot(2,1,1)
plot(x(1,:))
hold on
plot(x(2,:),'r')
plot(x(3,:),'k')
plot(x(4,:),'g')
xlabel('Pontos')
ylabel('Valores de x')
title('Algumas Amostras')
ylim([-6 6])

subplot(2,1,2)
plot(Em1)
xlabel('Pontos')
ylabel('Valor Esperado')
title('Primeiro Momento Estat�stico (valor esperado entre as 500 amostras)')
ylim([-6 6])
clear Em1

% Note que o valor esperado entre as amostras � praticamente zero para
% todos os pontos. De fato, podemos provar facilmente que, um processo AR1 
% do tipo estacion�rio tem valor esperado zero (voc� consegue?).

% Vamos agora calcular o nosso segundo momento estat�stico. Veja que este
% momento estat�stico � uma fun��o da delay (j) entre duas amostras. Fa�amos
% a seguinte defini��o do vetor rho(j,i) 
% rho(j,i) = E[X(i), X(i+j-1)]. 
% O primeiro elemento deste vetor (j=1) corresponde ao valor esperado 
% (sempre entre amostras) do quadrado do sinal (pot�ncia m�dia) em cada ponto.
% O segundo elemento (j=2) corresponde ao valor esperado do produto entre
% um determinado ponto e o ponto sucessivo (delay 1). E assim por diante.
% Vamos calcular esse produto com um delay que vai de 0 at� 99 pontos.


% rho � uma matriz com 100 linhas e 900 colunas. Cada linha representa uma
% amostra do processo AR1. As colunas representam os pontos dessas 
% sequencias. 
% 
% A variavel x tem 500 linhas e 1000 colunas, logo s�o 500 amostras do 
% processo AR1 onde cada amostra possui 1000 pontos do sinal. Nas 
% pr�ximas linhas de c�digo ser� feito o seguinte: 
% Calcularemos as autocorrela��es entre os pontos "atuais", i.e. pontos de 
% 1:900, com os pontos "futuros" de 1+j:900+j. Isso � feito multiplicando a
% matriz "atual" elemento a elemento com a matriz "pr�xima" e calculando o
% valor esperado em cada coluna. O vetor resultante � armazenado na linha
% da matriz rho.
% 
% Em cada coluna de rho teremos as autocorrela��es daquele ponto com suas 
% amostras sucessivas

rho = zeros(100,900);
for j=1:100 %delays
    atual = x(:,1:1000-100);
    proximo = x(:,1+j:1000+j-100);
    rho(j,:) = mean(atual.*proximo,1); %delay vai de 0 at� 99 amostras
end
% Observe algumas autocorrela��es:
figure;
for i=1:1:4
    subplot(4,1,i)
    a=round(900*rand);
    plot(0:99,rho(:,a))
    xlabel('Delay')
    ylim([-2 4])
    title(['Autocorrela��o (ponto ',num2str(a),' com suas amostras sucessivas)'])
end
clear a i j rho
% Claramente estas autocorrela��es s�o muito pouco dependentes do ponto em
% quest�o. O sinal gerado pelo processo AR1 � estacion�rio (no sentido fraco).

% Sinais n�o estacion�rios s�o aqueles em que coisas espec�ficas ocorrem em
% tempos espec�ficos. Por exemplo, considere o seguinte processo Y, amostrado 500
% vezes, e que poderia ser resultado de uma estimula��o de um nervo mediano :
y = 0.2*randn(500,1000);
y(:,300:599) = y(:,300:599)+repmat(sind(linspace(0,360,300)),500,1);
%Veja as primeiras quatro amostras deste processo:
figure;
for i=1:1:4
    subplot(4,1,i)
    plot(y(i,:));
    xlabel('Pontos')
    ylabel(['y',num2str(i)]);
    title(['Amostra ',num2str(i),' do processo Y'])
end
clear i

% Note agora que o c�lculo do primeiro momento indica uma clara 
% n�o-estacionariedade:
close all
Em1=mean(y,1);
figure('Color','w');

subplot(2,1,1)
plot(y(1,:))
hold on
plot(y(2,:),'r')
plot(y(3,:),'k')
plot(y(4,:),'g')
xlabel('Pontos')
ylabel('Valores de y')
title('Algumas Amostras')
a=ylim;

subplot(2,1,2)
plot(Em1)
xlabel('Pontos')
ylabel('Valor Esperado')
title('Valor Esperado do Primeiro Momento Estat�stico (entre as 500 amostras)')
ylim(a)
clear a Em1 ans

%% TEOREMA ERG�DICO:
close all
% Um sinal estacion�rio (no sentido fraco) � erg�dico (no sentido fraco) se
% podemos estimar os valores esperados dos primeiros dois momentos E[X_1], 
% E[X_1, X_2] por estat�sticas realizadas no tempo (ao longo dos pontos):
% 
% E[X1]= m�dia no tempo de X_1
% E[X1, X2]= convolu��o no tempo de X_1 e X_2

% O teorema erg�dico afirma que todo sinal estacion�rio no sentido fraco,
% que possui uma autocorrela��o que vai a zero para delays infinitos, pode
% ser considerado erg�dico. Este � o caso, por exemplo, do nosso processo 
% AR1. Ele � estacion�rio e a autocorrela��o � alta apenas com delay baixo,
% logo pode ser dito erg�dico (no sentido fraco). Isso significa que ao
% inv�s de calcularmos m�dias entre as amostras, podemos estimar essas
% m�dias calculando-as no tempo. Veja:

figure('Color','w');

subplot(2,1,1)
hist(mean(x(:,300:800),1))
title('M�dias entre as amostras para 500 pontos no tempo do processo X')
xlim([-2 2])

subplot(2,1,2)
hist(mean(x,2))
title('M�dias entre os tempos para 500 amostras do processo X')
xlim([-2 2])

% A m�dia entre as amostras � claramente melhor que a m�dia entre os tempos
% (o valor correto da m�dia desse processo � zero e portanto a m�dia entre
% amostras aproxima melhor esse valor do que a m�dia no tempo). Mas ambas 
% as distribui��es s�o bem parecidas (aproximadamente normais, centradas 
% em zero). Para reduzir a vari�ncia da distribui��o de m�dias no tempo, 
% basta aumentar o tamanho das amostras.
% 
% Note o que aconteceria se fiz�ssemos o mesmo para o nosso processo Y:

figure;

subplot(2,1,1)
hist(mean(y(:,300:800),1))
title('M�dias entre as amostras para 500 pontos no tempo do processo Y')
xlim([-2 2])

subplot(2,1,2)
hist(mean(y,2))
title('M�dias entre os tempos para 500 amostras do processo Y')
xlim([-2 2])

% Nesse caso, a m�dia no tempo � uma p�ssima aproxima��o para a m�dia entre 
% amostras, at� porque a m�dia entre as amostras VARIA NO TEMPO (n�o � 
% estacion�ria)

close all

% A Ergodicidade � extremamente importante para objetivos diagn�sticos. Ela
% possibilita que estimemos as propriedades de um processo estoc�stico a
% partir de uma �nica s�rie bem longa desse processo. Quanto mais longa a
% s�rie, melhor ela representar� o processo.

% Vamos inspecionar essas propriedades em um sinal real. Trata-se de um EEG
% espont�neo obtido em um canal ocipital de um sujeito saud�vel (e
% acordado). O sinal encontra-se na vari�vel aula09-EEG.mat:
load('aula09-EEG.mat')

% Voc� agora tem no seu workspace uma vari�vel EEG com 87 amostras, cada uma
% de 1 segundo (1000 pontos por segundo), do EEG. Vamos visualizar algumas
% dessas amostras e calcular a m�dia (entre amostras):
Em1=mean(EEG,1); %note, essa m�dia n�o � no tempo, mas entre as diferentes amostras
figure;

subplot(2,1,1)
plot(EEG(1,:))
hold on
plot(EEG(2,:),'r')
plot(EEG(3,:),'k')
plot(EEG(4,:),'g')
xlabel('Pontos')
ylabel('\mu V')
title('Algumas Amostras do EEG')
ylim([-50 50])

subplot(2,1,2)
plot(Em1)
xlabel('Pontos')
ylabel('Valor Esperado (\mu V)')
title('Valor Esperado do Primeiro Momento Estat�stico (entre as 87 amostras)')
ylim([-50 50])

clear Em1

% Note como a m�dia parece bem estacion�ria. Vamos agora calcular algumas
% correla��es (delays entre 0 e 299ms):
rho=zeros(300,700);
for j=1:300 % delays
    rho(j,:)=mean(EEG(:,1:1000-300).*EEG(:,1+j:1000+j-300),1); % delay vai de 0 at� 299 amostras
end
figure;
for i=1:1:4
    subplot(4,1,i)
    a=round(700*rand);
    plot(0:299,rho(:,a))
    xlabel('Delay (ms)')
    title(['Autocorrela��o (ponto ',num2str(a),' com suas amostras sucessivas)'])    
end
clear a i j rho

% Novamente, as correla��es s�o bastante semelhantes entre pontos
% distintos (em termos da escala de flutua��o do pr�prio EEG). 
% Este EEG parece suficientemente estacion�rio (ao menos no sentido fraco). 

% Mas notem agora como essas correla��es flutuam e n�o v�o rapidamente a
% zero como no modelo AR1. Na verdade essas correla��es v�o a zero 
% para delays grandes, mas muito lentamente, e � preciso um registro muito 
% longo para perceber isso. O sinal pode ser considerado um sinal erg�dico,
% (aproximadamente) mas sua autocorrela��o � uma fun��o que tem um
% claro conte�do espectral. Inspecionando o gr�fico anterior, voc� poder�
% ver que os picos da autocorrela��o est�o espa�ados por cerca de 130 a 140 
% ms.  Isto corresponde a uma frequ�ncia de pouco mais de 7Hz (entre as 
% bandas teta e alfa), tipicamente observada no EEG ocipital.

% De fato, em um sinal estacion�rio, a autocorrela��o, uma fun��o par que 
% depende apenas do delay entre as amostras envolvidas, est� intimamente
% relacionada � s�rie temporal original: mostramos em aula que a 
% transformada de Fourier dessa autocorrela��o corresponde � densidade 
% espectral da s�rie! Este � um dos interesses em se encarar um sinal
% biol�gico como sendo resultante de um processo estoc�stico: encontramos 
% outros m�todos para calcular o seu espectro!

%% MODELOS LINEARES
close all
clear x y A alpha
% Mas ainda temos um problema: a qualidade da estimativa das propriedades 
% de um processo estoc�stico atrav�s de estat�sticas no tempo depende de
% termos uma s�rie longa. E a impossibilidade de registrarmos uma s�rie 
% muito longa � justamente o que nos levou a considerar estes m�todos 
% estoc�sticos. Precisar�amos descobrir algum modo de construir uma s�rie 
% longa a partir de uma s�rie curta.
% A ideia � essa: fazer um modelo de predi��o. 

% Suponha por exemplo que o nosso processo seja tal que ele possa ser
% modelado pela seguinte rela��o:
% MODELO ARMA:
% x_t = soma_{k=1 at� p} A_k * x_{t-k} + soma_{k=0 at� q} B_k * r_{t-k} 

% onde A_k e B_k s�o par�metros (os par�metros do modelo) e r_i � um sinal
% aleat�rio (ru�do branco). Esse modelo � dito modelo ARMA de ordem 
% (p,q) (ordem = n�mero de par�metros).

% Se o modelo tiver A_{k}=0, para todo k, ele � dito modelo MA (moving
% average) de ordem q. Quando o modelo tiver B_{k}=0 para todo k exceto k=0,
% ele � dito modelo AR (autoregressivo) de ordem p.

% Note que o PROCESSO AR1 que definimos no in�cio desse arquivo � um
% processo AR simples de ordem 1. 

% Processos do tipo AR podem ser encarados como filtros IIR que possuem na
% na entrada um ruido branco r e na sa�da o sinal x. Ora, sabemos que em
% filtros IIR, o espectro da sa�da, X(omega), nada mais � do que o produto do
% espectro da entrada R(omega) com a fun��o transfer�ncia do filtro. Mas se
% a entrada for um ruido branco, o espectro da entrada � constante! 
% Portanto, se conseguirmos aproximar um sinal por um modelo AR, o espectro
% desse sinal pode ser aproximado pela fun��o transfer�ncia do filtro 
% associado ao modelo!

% Por exemplo, vamos gerar um filtro IIR de ordem 4 com os seguintes
% par�metros:
A = [1 -2.7607 3.8106 -2.6535 0.9238];
% onde A(1) = A_{k=0}, A(2) = A_{k=1},  etc... na express�o acima
% A fun��o de transfer�ncia desse filtro �
[H,w] = freqz(1,A);
figure;plot(w,20*log10(abs(H))); ylabel('Magnitude do espectro (dB)');xlabel('Frequencia normalizada');hold on;
% figure;plot(w,abs(H)); ylabel('Magnitude do espectro');xlabel('Frequencia normalizada');hold on;


% Como poder�amos estimar esse espectro a partir de uma amostra do processo
% AR correspondente? O m�todo Burg faz justamente isso: estima o espectro do
% modelo linear a partir de uma amostra desse processo.
% Para gerar uma amostra basta passar um ruido branco por esse filtro:
ruido = randn(1000,1); % entrada � o sinal aleatorio
x = filter(1,A,ruido); % saida � o sinal filtrado
% Agora observe como o resultado do m�todo BURG aproxima muito bem o espectro
% do processo:
[Pxx,F] = pburg(x,4);% estima do espectro 
plot(F,10*log10(Pxx))


%TEOREMA DE WOLD: TODA SEQU�NCIA ESTACION�RIA PODE SER REPRESENTADA POR UM
%MODELO ARMA, ou AR, ou MA. Veremos em aula a import�ncia desse resultado.

%Para encerrar, vamos calcular um modelo AR para o nosso EEG. Observe:
EEG=EEG-repmat(mean(EEG,2),1,size(EEG,2));
N=size(EEG,2); % numero de pontos do EEG
dt=1/size(EEG,2); % intervalo entre os pontos do EEG (em segundos)

trecho=10;

figure;

subplot(2,1,1)
Wfft=fft(EEG(trecho,:),N);
Pfft=Wfft.*conj(Wfft)/N;
fmax=ceil(N/2);
f=(0:fmax)/(N*dt);
plot(f,Pfft(1:length(f)))
xlim([0 40])
Pfft=Pfft(1:length(f));
ylabel('Densidade de Pot�ncia Espectral')
title('Espectro FFT para um trecho de 1 segundo do EEG')

subplot(2,1,2)
hold on

p=5; %ordem do modelo
[DEPar,f] = pburg(EEG(trecho,:),p,N*4,1/dt);
xa=plot(f,DEPar,'k','DisplayName','Ordem 5','LineWidth',2);

p=15; %ordem do modelo
[DEPar,f] = pburg(EEG(trecho,:),p,N*4,1/dt);
xa=plot(f,DEPar,'b','DisplayName','Ordem 15','LineWidth',2);

p=200; %ordem do modelo
[DEPar,f] = pburg(EEG(trecho,:),p,N*4,1/dt);
xa=plot(f,DEPar,'g','DisplayName','Ordem 200','LineWidth',2);

p=90; %ordem do modelo
[DEPar,f] = pburg(EEG(trecho,:),p,N*4,1/dt);
xa=plot(f,DEPar,'r','DisplayName','Ordem 90','LineWidth',2);

xlim([0 40])
ylabel('Densidade de Pot�ncia Espectral')
xlabel('Frequencia (Hz)')
legend show
title('Espectro do Modelo AR para um trecho de 1 segundo do EEG')

clear xa DEPar N Pfft Wfft dt f fmax p trecho ans