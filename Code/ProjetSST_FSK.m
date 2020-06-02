clear;
clc;
close all;

%------Parametres------%
N=1000;    % Nbre de bits       
F=500;     % Frequence symbole
T=1/F;     % Periode symbole
Ttot=N/F;  % Periode total
Fp=6000;   % Frequence porteuse    
Fs=20000;  % Frequence d'echantillonnage
ts=1/Fs;   % Periode d'echantillonnage
kf=500;    % Selectivite frequentielle
Nb=Fs/F;   % Nombre d'echantillons par bit

%-----Generation d'une liste de -1 et 1 de longueur N-----%

bits = randi(2,1,N)-1;
bits_msg = bits*2-1;   

%-----Creation du message avec oversampling-----%
       %On a Nb echantillons par bits
       
for i=1:N                      
    for t=1:Nb                 
        m(Nb*(i-1)+t)=bits_msg(i);
    end
end

%-----Generation de l'enveloppe complexe-----%

es = exp(1j*2*pi*kf*[0 cumsum(m(1:end-1))]/Fs);
        % NB : on divise par Fs a cause de l'integrale

%-----Signal Radio Frequence-----%
t=0:ts:Ttot-ts;
s = real(es).*cos(2*pi*Fp*t)-imag(es).*sin(2*pi*Fp*t);

%-----Densite spectrale des signaux avec Welch-----%

[psdm,f] = welch(m,t,1024,500);
[psds,f] = welch(s,t,1024,500);
[psdes,f] = welch(es,t,1024,500);

figure("Name","PSD Emetteur")
plot(f,10*log10(psdm),"b")
hold on
plot(f,10*log10(psdes),"r")
hold on
plot(f,10*log10(psds),'Color','#EDB120')
legend("m(t)","es(t)","s(t)")
xlabel("f[Hz]")
hold off

%----Ajout du bruit----%

r = s;


%-----Calcul de la bande-passante (Carlson)-----%
    % Amplitude = 1 => Deltaf = kf*1
    % BT = 2*Deltaf+2*Fm

BT = 2*kf + 2*F;

%----Filtre PB-----%
    % La bande-passante BT valant 2000, la frequence de coupure filtre PB 
    % doit donc etre 1000 (on choisit un peu plus en pratique, 1100)
Fcutoff = 1100;
order = 64;
PB = fir1(order,Fcutoff/(Fs/2),'low');        % On a besoin de la frequence normalisee
                                           % mais on prend Nyquist en
                                           % compte

%----Convolution du filtre----%
r1 = r.*cos(2*pi*Fp*t);
r2 = r.*sin(2*pi*Fp*t);
er1=conv(PB,r1);
er2=conv(PB,r2);
    % NB : produit de convolution discret => le nombre de samples vaut N +
    % ordre du filtre
er1cut= er1(order/2+1:Nb*N+order/2);
er2cut= er2(order/2+1:Nb*N+order/2);
er=2*(er1cut+er2cut*j);
    % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
    % multiplier par 2
[psder,f] = welch(er,t,1024,500);
    % On peut voir que l'enveloppe er(t) ressemble bien a l'enveloppe es(t)
figure('Name',"PSD Recepteur")

plot(f,10*log10(psder))
        % On a un parasite au milieu du spectre qui est du à du repliement
        % spectral. On peut y remédier en mettant un filtre passe-bande
        % plutot qu'un filtre passe-bas ou en augmentant la fréquence
        % d'échantillonnage. Mais ici c'est juste pour l'exemple.
        % Ce repliement spectral est situé en 2*Fp=12000 (symétrie en 10000
        % ce qui donne 8000);
 xlabel("f[Hz]")

%----- Demodulation FM-----%

a=1/(4*pi*kf);
    
    % Pour faire la dérivée on fait la variation avec :
s1=a*(([(er(2:end)-er(1:end-1))./(t(2:end)-t(1:end-1)) 0])+j*pi*BT*er(1:end));
s2=-a*(([(er(2:end)-er(1:end-1))./(t(2:end)-t(1:end-1)) 0])-j*pi*BT*er(1:end));

s0 = abs(s2)-abs(s1);
figure("Name","Comparaison signal envoyé/démodulé")
plot(t,s0)
hold on 
plot(t,m)
hold off
legend("s0","m")

for i = 1:N
    sum=0;
    for j=1:Nb
        sum=sum+s0(Nb*(i-1)+j);
    end
    %mf(i)=round(sum/Nb);
    
    if (sum/Nb) >= 0
        mf(i) = 1;
    else
        mf(i) = 0;
    end
end

figure("Name","Interprétation du message reçu")
plot((1:N),bits,'r.')
hold on
plot((1:N),mf,"Color",'#4DBEEE')
hold off
ylim([-0.5 1.5])
legend("bits","mf")

%---------------------------------%

