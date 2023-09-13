% Projete um filtro passa-altas Chebyshev que permita a passagem de sinais 
% acima de 165 rad/s com atenuação máxima de 2 dB e rejeite sinais abaixo 
% de 100 rad/s com atenuação mínima de 20 dB.

clc
clear
close all

%%Parâmetros do Passa-Altas
wp1 = 165;
ws1 = 100;
ap = 2;
as = 20;

%%Passo 1 - Determina o filtro passa-baixos prototipo
wp = 1;
ws = wp1/ws1;

% Passo 2 - Determinar a ordem do filtro passa-baixa protótipo
K = ceil(acosh(sqrt((10^(as/10)-1)/(10^(ap/10)-1)))/acosh(ws/wp));

% Passo 3 - Determina os Polos do Filtro Protótipo
epsilon = sqrt(10^(ap/10)-1);
ii = 1:K;
pk = -wp*sinh(asinh(1/epsilon)/K)*sin(pi*(2*ii-1)/(2*K))+...
     1j*wp*cosh(asinh(1/epsilon)/K)*cos(pi*(2*ii-1)/(2*K));
 
% Passo 4 - Determina a Função de Transferência do Filtro Prototipo
H0 = (mod(K,2)==1)+(mod(K,2)==0)/sqrt(1+epsilon^2);
B = H0*prod(-pk);
A = poly(pk);

% Passo 5 - Determinar os polos do Filtro Passa-Altas
pi_pa = wp1./pk;
B_pa = H0*poly(zeros(K,1));
A_pa = poly(pi_pa);

%%Passo 6 - Plotar a resposta em magnitude
omega = linspace(0,1000,1000);
H = polyval(B_pa,1j*omega)./polyval(A_pa,1j*omega);
subplot(2,1,1)
plot(omega,abs(H),'k','Linewidth',1.0)
grid on
xlabel('\omega (rad/s)')
ylabel('|{\itH}(\omega)|')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(2,1,2)
plot(omega,20*log10(abs(H)),'k','Linewidth',1.0)
grid on
xlabel('\omega (rad/s)')
ylabel('|{\itH}(\omega)| (dB)')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
