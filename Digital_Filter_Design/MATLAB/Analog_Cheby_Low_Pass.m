% Projete um filtro passa-baixas Chebyshev que permita a passagem de sinais 
% com até 10 rad/s com atenuação máxima de 2 dB e rejeite sinais acima de 
% 30 rad/s com atenuação mínima de 20 dB.

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

%%Passo 4 - Plotar a resposta em magnitude
omega = linspace(0,35,1000);
H = polyval(B,1j*omega)./polyval(A,1j*omega);
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
