clear;
clc;
close all;
%------Parametres------%

N=1000;    % Nbre de bits       
F=200;     % Frequence symbole
T=1/F;     % Periode symbole
Ttot=N/F;  % Periode totale
Fc=6000;   % Frequence porteuse    
Fs=20000;  % Frequence d'echantillonnage
OSF = Fs*T;% Oversampling factor
ts=1/Fs;   % Periode d'echantillonnage
kf=300;    % Selectivite frequentielle
noiseratio = 7;   %Rapport Eb/N0 en dB

%-----Generation d'une liste de -1 et 1 de longueur N-----%

bits = randi(2,1,N)-1;
bits_msg = bits*2-1;   

%-----Creation du message avec oversampling-----%
       %On a Nb echantillons par bits

m = zeros(1,OSF*N);
for i=1:N                      
    for t=1:OSF                
        m(OSF*(i-1)+t)=bits_msg(i);
    end
end

%-----Generation de l'enveloppe complexe-----%

es = exp(1j*2*pi*kf*[0 cumsum(m(1:end-1))]/Fs);
        % NB : on divise par Fs a cause de l'integrale

%-----Signal Radio Frequence-----%

t=0:ts:Ttot-ts;
s = real(es).*cos(2*pi*Fc*t)+imag(es).*sin(2*pi*Fc*t);
figure;
plot(t,s)
xlim([0 0.01])
hold off;
%-----Densite spectrale des signaux avec Welch-----%

[psdm,f] = welch(m,t,1024,500);
[psds,f] = welch(s,t,1024,500);
[psdes,f] = welch(es,t,1024,500);

figure("Name","PSD Emetteur")
%plot(f,10*log10(psdm),"b")
hold on
%plot(f,10*log10(psdes),"r")
hold on
plot(f,10*log10(psds),'Color','r')
%legend("m(t)","e_s(t)","s(t)")
legend("s(t)   (k_f=300Hz, F=200Hz)")
ylabel("PSD[dB]")
xlabel("f[Hz]")
hold off

%----Ajout du bruit blanc additif gaussien----%

   %On evalue la puissance Ps du signal sans bruit s(t)
 
Ps_in = sum(s.^2)/(N*OSF);
   %NB : en bande de base il faut diviser par 2

   %On evalue Eb (energie du bit) a partir du signal sans bruit s(t)

 Eb = Ps_in*T;

   %Le rapport Eb/N0 = 10dB
        
N0=Eb/(10^(noiseratio/10));

   %La puissance du bruit Gaussien est Pn=N0*Fs/2
       
Pn_in = N0*Fs/2;

   %Generation du bruit

n=sqrt(Pn_in)*randn(1,OSF*N);

r = s+n;

figure;
[psdr,f] = welch(r,t,1024,500);
plot(f,10*log10(psdr),"r");
hold on;
plot(f,10*log10(psds),"b:");
hold off;
ylabel("PSD[dB]")
xlabel("f[Hz]")
legend("PSD r(t)","PSD s(t)")
%-----Calcul de la bande-passante (Carlson)-----%
    % Amplitude = 1 => Deltaf = kf*1
    % BT = 2*Deltaf+2*Fm

BT = 2*kf + 2*F;

%----Filtre PB-----%
    % La bande-passante BT valant 2000, la frequence de coupure filtre PB 
    % doit donc etre 1000 (on choisit un peu plus en pratique, 1100)
    
Fcutoff = 1100;
order = 64;
PB = fir1(order,Fcutoff/(Fs/2),'low');     % On a besoin de la frequence normalisee
                                           % mais on prend Nyquist en
         % Spectre du filtre PB
    h=zeros(1,1024);
    [h(513:1024),w]=freqz(PB,1,512);
    h(1:512)=h(1024:-1:513);                                      % compte

%===== SIGNAL BRUITE ====%                                          
%----Convolution du filtre----%

r1 = r.*cos(2*pi*Fc*t);

r2 = r.*sin(2*pi*Fc*t);

er1=conv(PB,r1);
er2=conv(PB,r2);
    % NB : produit de convolution discret => le nombre de samples vaut N +
    % ordre du filtre
er1cut= er1(order/2+1:OSF*N+order/2);
er2cut= er2(order/2+1:OSF*N+order/2);
er=2*(er1cut+er2cut*j);
    % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
    % multiplier par 2
[psderf,f] = welch(er,t,1024,500);
[psder,f]=welch((r1+r2*j),t,1024,500);
    % On peut voir que l'enveloppe er(t) ressemble bien a l'enveloppe es(t)
figure('Name',"PSD Recepteur")
plot(f,10*log10(abs(h)),"color","[0, 0.75, 0.75]");hold on;
plot(f,10*log10(psderf),"color",'	r');
plot(f,10*log10(psder),'k--')
        % On a un parasite au milieu du spectre qui est du à du repliement
        % spectral. On peut y remédier en mettant un filtre passe-bande
        % plutot qu'un filtre passe-bas ou en augmentant la fréquence
        % d'échantillonnage. Mais ici c'est juste pour l'exemple.
        % Ce repliement spectral est situé en 2*Fc=12000 (symétrie en 10000
        % ce qui donne 8000);
 legend("H(f) LPF","PSD e_r(t) with LPF","PSD e_r(t) without LPF")
 xlabel("f[Hz]")

%---- Demodulation FM ----%
 
Ac = 1; % Amplitude du signal R*(t) en bande de base

a=1/(Ac*4*pi*kf);

    % Pour faire la derivee on fait la variation avec :
s1=a*([(er(2:end)-er(1:end-1)) 0]*Fs+j*pi*BT*er(1:end));
s2=-a*([(er(2:end)-er(1:end-1)) 0]*Fs-j*pi*BT*er(1:end));

mrn = abs(s1)-abs(s2);    %message receive noise


figure("Name","Comparaison signal envoyé/démodulé")
plot(t,mrn)
hold on 
plot(t,m,'r')
hold off
legend("m_r(t)","m_s(t)")
xlabel("t[s]")

%===== SIGNAL SANS BRUIT ====%                                          
%----Convolution du filtre----%

r1 = s.*cos(2*pi*Fc*t);
r2 = s.*sin(2*pi*Fc*t);
er1=conv(PB,r1);
er2=conv(PB,r2);
er1cut= er1(order/2+1:OSF*N+order/2);
er2cut= er2(order/2+1:OSF*N+order/2);
er=2*(er1cut+er2cut*j);

%---- Demodulation FM ----%

s1=a*([er(2:end)-er(1:end-1) 0]*Fs+j*pi*BT*er(1:end));
s2=-a*([er(2:end)-er(1:end-1) 0]*Fs-j*pi*BT*er(1:end));

mr = abs(s1)-abs(s2);  %message receive without noise

%===== Calcul du SNR =====%

SNRin = 10*log10(Ps_in/(N0*BT));

pure_noise = mrn-mr;    %bruit pur
Ps_out = trapz(mrn.^2)/(N*OSF);
Pn_out = trapz(pure_noise.^2)/(N*OSF);

SNRout = 10*log10(Ps_out /Pn_out);

    % Facteur de mérite 
FoM = SNRout-SNRin

    %Pour ameliorer le facteur de merite il faut filtrer le signal de 
    %facon a enlever les parasites qui apparaissent dans la PSD du 
    %signal recu
    
    %NB : il faut faire un plot Fom en fonction de Eb/N0 et pour chaque
    %point il faut faire une moyenne sur 100 FoM afin d'avoir la courbe la
    %plus précise
    
%==== Interpretation du message ====%
mf1=zeros(1,N);

for i = 1:N
    mf_sum=0;
    for j=1:OSF
        mf_sum=mf_sum+mrn(OSF*(i-1)+j);
    end
    if (mf_sum/OSF) >= 0
        mf1(i) = 1;
    else
        mf1(i) = 0;
    end
end

figure("Name","Interprétation du message reçu (FM)")
plot((1:N),bits,'r.')
hold on
plot((1:N),mf1,"Color",'#4DBEEE')
hold off
ylim([-0.5 1.5])
legend("bits TX","bits RX")

%==== Demodulation FSK ====%
s0 = exp(-1j*2*pi*kf*t);     
s1 = exp(1j*2*pi*kf*t);      


es0 = er.*conj(s0);
es1 = er.*conj(s1);

[psds0,f]=welch(es0,t,1024,500);
[psds1,f]=welch(es1,t,1024,500);

figure;
plot(f,10*log10(psds0));hold on;
plot(f,10*log10(psds1));hold off;
legend("PSD e_{r}(t)e'_{s0}(t)","PSD e_{r}(t)e'_{s}(t)");
xlabel("f[Hz]");
ylabel("PSD[dB]");

%Correlation
mf2 = zeros(1, N);
for k = 1 : N 
        sum_0 = abs(sum(es0(((k-1)*OSF)+1:OSF*k)));
        sum_1 = abs(sum(es1(((k-1)*OSF)+1:OSF*k)));
        if sum_0 < sum_1
            mf2(k) = 1;
        end
end

figure("Name","Interprétation du message reçu (FSK)")
plot((1:N),bits,'r.')
hold on
plot((1:N),mf2,"Color",'#4DBEEE')
hold off
ylim([-0.5 1.5])
legend("bits","mf2")

