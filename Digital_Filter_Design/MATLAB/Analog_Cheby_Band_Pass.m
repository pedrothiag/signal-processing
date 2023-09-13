% Projete um filtro de Chebyshev passa-faixas que atenda as seguintes 
% especificações:
% - Permitir a passagem de sinais entre 1000 rad/s e 2000 rad/s com 
%   atenuação máxima de 1 dB;
% - Rejeitar sinais abaixo de 450 rad/s com atenuação mínima de 20 dB;
% - Rejeitar sinais acima de 4000 rad/s com atenuação mínima de 20 dB.

clc
clear
close all

%%Parâmetros do Passa-Faixas
wp1 = 1000;
wp2 = 2000;
ws1 = 450;
ws2 = 4000;
ap = 1;
as = 20;

%%Passo 1 - Determinar a frequencia da banda de rejeicao do passa-baixa
%%prototipo
ws_p_1 = (wp1*wp2 - ws1^2)/(ws1*(wp2 - wp1));
ws_p_2 = (ws2^2 - wp1*wp2)/(ws2*(wp2 - wp1));
ws_p = min(ws_p_1,ws_p_2);

% Passo 2 - Determinar a ordem do filtro passa-baixa protótipo
K = ceil(acosh(sqrt((10^(as/10)-1)/(10^(ap/10)-1)))/acosh(ws_p/1));

% Passo 3 - Determina os Polos do Filtro Protótipo
epsilon = sqrt(10^(ap/10)-1);
ii = 1:K;
pk = -1*sinh(asinh(1/epsilon)/K)*sin(pi*(2*ii-1)/(2*K))+...
     1j*1*cosh(asinh(1/epsilon)/K)*cos(pi*(2*ii-1)/(2*K));
 
% Passo 4 - Determina a Função de Transferência do Filtro Prototipo
H0 = (mod(K,2)==1)+(mod(K,2)==0)/sqrt(1+epsilon^2);
B = H0*prod(-pk);
A = poly(pk);

%%Passo 5 - Determinar os polos do filtro passa-faixa
a = 1;
b = -pk*(wp2 - wp1);
c = wp1*wp2;
pbar = [(-b + sqrt(b.^2 - 4*a*c))./(2*a),(-b - sqrt(b.^2 - 4*a*c))./(2*a)];

%%Passo 6 - Determinar funcao transferencia do passa-faixa
B = (H0)*((wp2-wp1)^K)*poly(zeros(K,1));
A = real(poly(pbar));

%%Passo 7 - Plotar a resposta em magnitude
omega = linspace(0,8000,8000);
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
