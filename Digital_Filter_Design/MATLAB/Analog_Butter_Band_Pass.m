% Projete um filtro de Butterworth passa-faixas que atenda as seguintes 
% especificações:
% - Permitir a passagem de sinais entre 10 rad/s e 20 rad/s com atenuação 
%   máxima de 2 dB;
% - Rejeitar sinais abaixo de 5 rad/s com atenuação mínima de 20 dB;
% - Rejeitar sinais acima de 35 rad/s com atenuação mínima de 20 dB.

clc
clear
close all

%%Parâmetros do Passa-Faixas
wp1 = 2*pi*30;
wp2 = 2*pi*120;
ws1 = 2*pi*10;
ws2 = 2*pi*200;
ap = 1;
as = 20;

%%Passo 1 - Determinar a frequencia da banda de rejeicao do passa-baixa
%%prototipo
ws_p_1 = (wp1*wp2 - ws1^2)/(ws1*(wp2 - wp1));
ws_p_2 = (ws2^2 - wp1*wp2)/(ws2*(wp2 - wp1));
ws_p = min(ws_p_1,ws_p_2);

%%Passo 2 - Determinar a ordem do filtro passa-baixa protótipo e a sua
%%frequencia de corte
wp_p = 1;
K = ceil(log10((10^(as/10)-1)/(10^(ap/10)-1))/(2*log10(ws_p/wp_p)));
wc_1 = wp_p/(10^(ap/10)-1)^(1/(2*K));
wc_2 = ws_p/(10^(as/10)-1)^(1/(2*K));
wc = (wc_1 + wc_2)/2;

%%Passo 3 - Determinar os polos do filtro passa-baixa protótipo
ii = 1:K;
pk = 1j*wc*exp(1j*pi/(2*K)*(2*ii-1));
Q = real(poly(pk));

%%Passo 4 - Determinar os polos do filtro passa-faixa
a = 1;
b = -pk*(wp2 - wp1);
c = wp1*wp2;
pbar = [(-b + sqrt(b.^2 - 4*a*c))./(2*a),(-b - sqrt(b.^2 - 4*a*c))./(2*a)];

%%Passo 5 - Determinar funcao transferencia do passa-faixa
B = ((wc^K)*(wp2-wp1)^K)*poly(zeros(K,1));
A = real(poly(pbar));

%%Passo 6 - Plotar a resposta em magnitude
omega = linspace(0,2*pi*300,1000);
H = polyval(B,1j*omega)./polyval(A,1j*omega);
subplot(2,1,1)
plot(omega./(2*pi),abs(H),'k','Linewidth',1.0)
grid on
xlabel('{\itf} (Hz)')
ylabel('|{\itH}({\itf})|')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(2,1,2)
plot(omega./(2*pi),20*log10(abs(H)),'k','Linewidth',1.0)
grid on
xlabel('{\itf} (Hz)')
ylabel('|{\itH}({\itf})| (dB)')
set(gca,'FontName','Times')
set(gca,'FontSize',10)