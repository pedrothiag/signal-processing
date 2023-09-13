% Projete um filtro de Butterworth rejeita-faixas que atenda as seguintes 
% especificações:
% - Rejeitar sinais entre 100 rad/s e 200 rad/s com atenuação mínima de 
%   20 dB;
% - Permitir a passagem de sinais abaixo de 60 rad/s com atenuação máxima 
%   de 2.2 dB;
% - Permitir a passagem de sinais acima de 260 rad/s com atenuação máxima 
%   de 2.2 dB.

clc
clear
close all

%%Parâmetros do Passa-Faixas
wp1 = 60;
wp2 = 260;
ws1 = 100;
ws2 = 200;
ap = 2.2;
as = 20;

%%Passo 1 - Determinar a frequencia da banda de rejeicao do passa-baixa
%%prototipo
ws_p_1 = (ws1*(wp2 - wp1))/(wp1*wp2 - ws1^2);
ws_p_2 = (ws2*(wp2 - wp1))/(ws2^2 - wp1*wp2);
ws_p = min(ws_p_1,ws_p_2);

%%Passo 2 -Determinar a ordem do filtro passa-baixa protótipo e a sua
%%frequencia de corte
wp_p = 1;
K = ceil(log10((10^(as/10)-1)/(10^(ap/10)-1))/(2*log10(ws_p/wp_p)));
wc_1 = wp_p/(10^(ap/10)-1).^(1/(2*K));
wc_2 = ws_p/(10^(as/10)-1).^(1/(2*K));
wc = (wc_1 + wc_2)/2;

%%Passo 3 - Determinar os polos do filtro passa-baixa protótipo
ii = 1:K;
pk = 1j*wc*exp(1j*pi/(2*K)*(2*ii-1));
Q = real(poly(pk));
 
%%Passo 4 - Determinar os polos do filtro rejeita-faixa
a = 1;
b = (wp2 - wp1)./(-pk);
c = wp1*wp2;
pbar = [(-b + sqrt(b.^2 - 4.*a.*c))./(2.*a), ...
    (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a)];

%%Passo 5 - Determinar os zeros do filtro rejeita-faixa
zersimple = [-sqrt(-wp1*wp2),sqrt(-wp1*wp2)];
zer = repmat(zersimple,1,K);

%%Passo 6 - Determinar funcao transferencia do passa-faixa
A = real(poly(pbar));
B = real(poly(zer));

%%Passo 7 - Plotar a resposta em magnitude
omega = linspace(0,300,3001);
H = polyval(B,1j*omega)./polyval(A,1j*omega);
subplot(2,1,1)
plot(omega,abs(H),'k','Linewidth',1.0)
grid on
xlabel('\Omega (rad/s)')
ylabel('|{\itH}({\itj}\Omega)|')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(2,1,2)
plot(omega,20*log10(abs(H)),'k','Linewidth',1.0)
grid on
xlabel('\Omega (rad/s)')
ylabel('|{\itH}({\itj}\Omega)| (dB)')
set(gca,'FontName','Times')
set(gca,'FontSize',10)