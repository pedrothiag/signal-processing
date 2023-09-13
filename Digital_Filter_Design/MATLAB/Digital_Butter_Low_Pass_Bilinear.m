% Projete um filtro digital recursivo de forma a realizar um filtro de 
% Butterworth passa-baixas de primeira ordem com frequência de corte de 
% $\Omega_c = 10^5$ rad/s. Utilize a transformação bilinear.

clc
clear
close all

%Parametros do Filtro e Funcao de Transferencia
Omega_c = 1e5;
num = [Omega_c];
den = [1 Omega_c];
Ha = tf(num,den);

%Determina o intervalo de amostragem (para comparar, vamos utilizar o mesmo
%valor da invariância ao impulso)
Omega = linspace(1e5,5e6,1e6);
Ha_abs = abs(polyval(num,1j*Omega)./polyval(den,1j*Omega));
Omega_0 = Omega(find(Ha_abs<0.1,1));
T = pi/Omega_0;

%Discretiza o Sistema
Hz = c2d(Ha,T,'tustin');

%Diagrama de Polos e Zeros
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

%Resposta em Frequência
[h,w] = freqz(b,a,'whole',2001);
Omega = linspace(0,1e6,1e6);
H = polyval(num,1j*Omega)./polyval(den,1j*Omega);
figure,
plot(w*1/T,abs(h),'k','Linewidth',1.0)
hold on
plot(Omega,abs(H),'--k','Linewidth',1.0)
xlim([min(Omega) max(Omega)])
grid on
xlabel('\omega (rad/s)')
legend('Digital','Analógico')
set(gca,'FontName','Times')
set(gca,'FontSize',10)