clear;
clc;
close all;
slCharacterEncoding('UTF-8')

%-----Generation d'une liste de -1 et 1 de longueur N-----%
text_tx =  fileread('text.txt');
text_tx_int = uint8(double(text_tx)) ; % 8 bits/character
text_tx_bin = de2bi(text_tx_int,8) ; % decimal to binary
nbchar = size(text_tx_bin,1) ;
bits = double(reshape(text_tx_bin',1,nbchar*8)) ;

bits_msg = 2*bits-1;
%------Parametres------%
N=length(bits_msg);      
F=500;     % Frequence symbole
T=1/F;     % Periode symbole
Ttot=N/F;  % Periode totale
Fc=6000;   % Frequence porteuse    
Fs=20000;  % Frequence d'echantillonnage
OSF = Fs*T;% Oversampling factor
ts=1/Fs;   % Periode d'echantillonnage
kf=500;    % Selectivite frequentielle
noiseratio = 50;   %Rapport Eb/N0 en dB
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
                                           % compte

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
ern=2*(er1cut+er2cut*j);
    % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
    % multiplier par 2
[psder,f] = welch(ern,t,1024,500);
    % On peut voir que l'enveloppe er(t) ressemble bien a l'enveloppe es(t)
figure('Name',"PSD Recepteur")

plot(f,10*log10(psder))
        % On a un parasite au milieu du spectre qui est du à du repliement
        % spectral. On peut y remédier en mettant un filtre passe-bande
        % plutot qu'un filtre passe-bas ou en augmentant la fréquence
        % d'échantillonnage. Mais ici c'est juste pour l'exemple.
        % Ce repliement spectral est situé en 2*Fc=12000 (symétrie en 10000
        % ce qui donne 8000);
 legend("er(t)")
 xlabel("f[Hz]")

%---- Demodulation FM ----%
 
Ac = 1; % Amplitude du signal R*(t) en bande de base

a=1/(Ac*4*pi*kf);

    % Pour faire la derivee on fait la variation avec :
s1=a*([(ern(2:end)-ern(1:end-1)) 0]*Fs+j*pi*BT*ern(1:end));
s2=-a*([(ern(2:end)-ern(1:end-1)) 0]*Fs-j*pi*BT*ern(1:end));

mrn = abs(s1)-abs(s2);    %message receive noise

figure("Name","Comparaison signal envoyé/démodulé")
plot(t,mrn)
hold on 
plot(t,m,'r')
hold off
legend("mrn(t)","m(t)")
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

pure_noise = mrn-mr;    %bruit pure
Ps_out = sum(mrn.^2)/(N*OSF);
Pn_out = sum(pure_noise.^2)/(N*OSF);

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
legend("bits","mf1")

%==== Demodulation FSK ====%
s0 = exp(-1j*2*pi*kf*t);     
s1 = exp(1j*2*pi*kf*t);      

es0 = ern.*conj(s0);
es1 = ern.*conj(s1);

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

text_rx_bin = reshape(mf2,8,nbchar)' ; % nbchar : number of characters
text_rx_int = bi2de(text_rx_bin)' ; % binary to decimal
text_rx = char(text_rx_int)  % decimal to text