% Projete um filtro digital passa-baixas Chebyshev que permita a passagem 
% de sinais com até 10 rad/s com atenuação máxima de 2 dB e rejeite sinais 
% acima de  30 rad/s com atenuação mínima de 20 dB. Utilize invariancia ao
% impulso

clc
clear
close all

%%Parâmetros do Passa-Baixas
wp = 10;
ws = 30;
ap = 2;
as = 20;

% Passo 1 - Determina a ordem do Filtro de Chebyshev
K = ceil(acosh(sqrt((10^(as/10)-1)/(10^(ap/10)-1)))/acosh(ws/wp));

% Passo 2 - Determina os Polos da Funcao de Transferencia
epsilon = sqrt(10^(ap/10)-1);
ii = 1:K;
pk = -wp*sinh(asinh(1/epsilon)/K)*sin(pi*(2*ii-1)/(2*K))+...
     1j*wp*cosh(asinh(1/epsilon)/K)*cos(pi*(2*ii-1)/(2*K));
 
% Passo 3 - Determina a Função de Transferência
H0 = (mod(K,2)==1)+(mod(K,2)==0)/sqrt(1+epsilon^2);
B = H0*prod(-pk);
A = poly(pk);

% Passo 4 - Cria a funcao de transferência
sysc = tf(B,A);

% Passo 5 - Determina o Intervalo de Amostragem
Omega = linspace(20,100,200);
Ha_abs = abs(polyval(B,1j*Omega)./polyval(A,1j*Omega));
Omega_0 = Omega(find(Ha_abs<0.05,1));
Ts = pi/Omega_0;

% Passo 6 - Discretizacao
Hz = c2d(sysc,Ts,'impulse');

% Passo 7 - Constante de Normalizacao
Hz_1 = evalfr(Hz,1);
Ha_0 = abs(polyval(B,0)./polyval(A,0));
K = Ha_0/Hz_1;
Hz = K*Hz;

% Passo 8 - Diagrama de Polos e Zeros e Resposta em Frequência
b = cell2mat(Hz.numerator);
a = cell2mat(Hz.denominator);
[hz1, hp1, ht1] =  zplane(b,a);
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
xlabel('Re\{{\itz}\}')
ylabel('Im\{{\itz}\}')
set(gca,'FontName','Times')
set(gca,'FontSize',10)

% Passo 9 - Plota a Resposta em Frequência
[h,w] = freqz(b,a,'whole',2001);
omega = linspace(0,35,1000);
H = polyval(B,1j*omega)./polyval(A,1j*omega);

figure,
plot(w*1/Ts,abs(h),'k','Linewidth',1.0)
hold on
plot(omega,abs(H),'--k','Linewidth',1.0)
xlim([min(omega) max(omega)])
grid on
xlabel('\omega (rad/s)')
legend('Digital','Analógico')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
