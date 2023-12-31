% sistema_atraso_grupo.m
% Esse código simula um sistema cujo o atraso de grupo não é constante,
% gerando, assim, uma distorção na fase. Para mais detalhes, consultar
% Cap. 5 do livro Processamento em Tempo Discreto de Sinais do Oppenheim
% & Schafer.

clc
clear
close all

% Vetor de amostras e janela (envelope)
n = 0:60;
w = 0.54 - 0.46*cos(2*pi*n/60);

% Sinais em banda passante
x1 = w.*cos(0.2*pi*n);
x2 = w.*cos(0.4*pi*n - pi/2);
x3 = w.*cos(0.8*pi*n + pi/5);

% Gera a sequencia de pulsos e plota o espectro
x = [x3 x1 x2 zeros(1,520)];
Xk = fftshift(fft(x));
N = length(x);
omega = linspace(-pi,pi,N);
figure,
set(gcf,'Position',[100 100 700 400])
subplot(211)
plot(x,'k','Linewidth',1.0)
xlim([0 300])
xlabel('\itn')
ylabel('{\itx}[{\itn}]')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(212)
plot(omega,abs(Xk),'k','Linewidth',1.0);
xticks([-pi -0.8*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.8*pi pi]);
xticklabels({'\pi', '-0,8\pi', '-0,4\pi', '-0,2\pi', '0', '0,2\pi', '0,4\pi', '0,8\pi', '\pi'})
xlim([0 pi])
xlabel('\omega')
ylabel('|{\itX}({\ite^{j\omega}})|')
set(gca,'FontName','Times')
set(gca,'FontSize',10)

% Gera o sistema com atraso de grupo não-constante. 
k = 1:4;
ck = 0.95*exp(1j*(0.15*pi + 0.02*pi*k));
z1 = 0.98*exp(1j*0.8*pi);
p1 = 0.8*exp(1j*0.4*pi);
zeros = [z1 conj(z1) 1./ck conj(1./ck) 1./ck conj(1./ck)];
poles = [p1 conj(p1) ck conj(ck) ck conj(ck)];
k = 0.1;
sys = zpk(zeros,poles,k,0.1);
[b,a] = zp2tf(zeros',poles',k);

% Mostra o diagrama de polos e zeros do sistema
figure,
set(gcf,'Position',[100 100 700 400])
[hz1, hp1, ht1] = zplane(b,a);
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
xlabel('Re\{{\itz}\}')
ylabel('Im\{{\itz}\}')
set(gca,'FontName','Times')
set(gca,'FontSize',10)

% Mostra a resposta em frequencia do sistema
[Hk,w] = freqz(b,a,'whole',2001);
tg = -1.0*diff(unwrap(angle(Hk)))*320;
figure,
set(gcf,'Position',[100 100 700 400])
subplot(311)
plot(w,abs(Hk),'k','Linewidth',1.0);
xticks([-pi -0.8*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.8*pi pi]);
xticklabels({'\pi', '-0,8\pi', '-0,4\pi', '-0,2\pi', '0', '0,2\pi', '0,4\pi', '0,8\pi', '\pi'})
xlim([0 pi])
xlabel('\omega')
ylabel('|{\itH}({\ite^{j\omega}})|')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(312)
plot(w,unwrap(angle(Hk)),'k','Linewidth',1.0);
xticks([-pi -0.8*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.8*pi pi]);
xticklabels({'\pi', '-0,8\pi', '-0,4\pi', '-0,2\pi', '0', '0,2\pi', '0,4\pi', '0,8\pi', '\pi'})
xlim([0 pi])
xlabel('\omega')
ylabel('\angle{\itH}({\ite^{j\omega}})')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(313)
plot(w(1:end-1),tg,'k','Linewidth',1.0)
xticks([-pi -0.8*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.8*pi pi]);
xticklabels({'\pi', '-0,8\pi', '-0,4\pi', '-0,2\pi', '0', '0,2\pi', '0,4\pi', '0,8\pi', '\pi'})
xlim([0 pi])
xlabel('\omega')
ylabel('{\itn_g}(\omega)')
set(gca,'FontName','Times')
set(gca,'FontSize',10)

% Gera a saída
y = filter(b,a,x);
figure,
set(gcf,'Position',[100 100 700 400])
subplot(211)
plot(y,'k','Linewidth',1.0)
xlim([0 300])
xlabel('{\itn}')
ylabel('{\ity}[{\itn}]')
set(gca,'FontName','Times')
set(gca,'FontSize',10)
subplot(212)
plot(x,'k','Linewidth',1.0)
xlim([0 300])
xlabel('{\itn}')
ylabel('{\itx}[{\itn}]')
set(gca,'FontName','Times')