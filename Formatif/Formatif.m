%% Problème 6
close all
clc

Num1 = [1 3 23];
Den1 = [1 5 22 7 9];

G = tf(Num1, Den1)
%% a)
%Nombre de mdoes
Pol = roots(Den1)
%% b)
wn = abs(Pol)
zeta = -real(Pol)./wn

[r, p, k] = residue(Num1,Den1)
poids = abs(r)./abs(real(p))

Phi = acosd(zeta)
Mp = 100*exp(-pi.*tand(Phi))
wa = wn.*sqrt(1-zeta.^2)
tp = pi./wa
ts = 4./(zeta.*wn)

%% c)
%Grand K
pol = rlocus(G, 5000)

figure
plot(real(pol), imag(pol), 'p')
hold on
rlocus(G)

ts = 4./real(pol)

%Regle 5
zero = roots(Num1)
temp = (sum(pol)-sum(zero))/2
ts2 = 4/abs(temp)

%% d)
%Graphiquement j'ai trouver 2.04 et 40

%% e)
figure
margin(G)

%Augementer de 15° à 3.51rad/s

%% f)
sys1 = tf([r(3)], [1 -p(3)])
sys2 = tf([r(4)], [1 -p(4)])


TF = sys1+sys2

figure
hold on
rlocus(G, 'r')
rlocus(TF, 'k')

%% Problème 7
close all

%% a)


%% Numero 8
clc
clear all
close all

Num1 = [1 3 5];
Den1 = [1 7 5 3 0];

TF = tf(Num1, Den1)
%% a)
mode = roots(Den1)

%% b)
figure
margin(TF)

K = db2mag(-8.5)

%% c)
TF2 = (K)*TF

figure
margin(TF2)

%% d)
% On diminue le K

%% e) 

[Gm,Pm,Wcg,Wcp] = margin(TF2)

ans = (Pm/Wcp)*(pi/180)

