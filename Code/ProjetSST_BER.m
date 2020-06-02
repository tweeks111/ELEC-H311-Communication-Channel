clear;
clc;
close all;
%---------------------
%   Projet SST
%  
%   Théo LEPOUTTE
%
%------Parametres------

N=1000;    % Nbre de bits       
F=500;     % Frequence symbole
T=1/F;     % Periode symbole
Ttot=N/F;  % Periode totale
Fc=6000;   % Frequence porteuse    
Fs=20000;  % Frequence d'echantillonnage
OSF = Fs*T;% Oversampling factor
ts=1/Fs;   % Periode d'echantillonnage
EbtoN0final = 20;
scale = 0:0.5:EbtoN0final;
            % Rapport Eb/N0 en dB
                    % Si Eb/N0 = 1 -> [Eb/N0]dB= 0 dB
                    % Un recepteur parfait aurait un Eb/N0= +Inf
                    % Eb/N0 -> 0 dB => Bruit tres eleve
                    % Eb/N0 -> Inf dB => Bruit tres faible

%-----Generation d'une liste de -1 et 1 de longueur N-----%

bits_tx = randi(2,1,N)-1;
bits_msg = bits_tx*2-1;

%-----Creation du message avec oversampling-----%
       %On a OSF echantillons par bits

m = zeros(1,OSF*N);
for i=1:N                      
    for j=1:OSF                
        m(OSF*(i-1)+j)=bits_msg(i);
    end
end



BER_FM = zeros(3,length(scale));
BER_FSK = zeros(3,length(scale));
BER_FSK_COH = zeros(3,length(scale));
FOM = zeros(3,length(scale));

kfiter=1;
for kf = [400,500,600]
    
%-----Generation de l'enveloppe complexe-----%
    es = exp(1j*2*pi*kf*[0 cumsum(m(1:end-1))]/Fs);
        % NB : on divise par Fs a cause de l'integrale

%-----Signal en bande passante-----%

    t=0:ts:Ttot-ts;
    s = real(es).*cos(2*pi*Fc*t)+imag(es).*sin(2*pi*Fc*t);

%----Ajout du bruit blanc additif gaussien----%
   %On evalue la puissance Ps du signal sans bruit s(t)
 
    Ps_in = trapz(s.^2)/(N*OSF);

   %On evalue Eb (energie du bit) a partir du signal sans bruit s(t)
    Eb = Ps_in*T;

%----Filtre PB-----%
    % Calcul de la bande passante (Carson)
    BT = 2*kf + 2*F;
    % Frequence de coupure du filtre
    Fcutoff = 1.1*BT/2;
    order = 64;
    PB = fir1(order,Fcutoff/(Fs/2),'low');     % On a besoin de la frequence normalisee
                                               % On prend Nyquist en compte
                                               
    % Spectre du filtre PB
    h=zeros(1,1024);
    [h(513:1024),w]=freqz(PB,1,512);
    h(1:512)=h(1024:-1:513);

    iter=1;
    for EbtoN0 = scale
        BER_FM_SUM=0;
        BER_FSK_SUM=0;
        BER_FM_COH_SUM=0;
        BER_FSK_COH_SUM=0;
        FOM_SUM=0;
        
        % Moyenne du BER et de la FOM sur Nb repetitions
        repetition = 20;
        for i=1:1:repetition
            [N0,Pn_in,SNRin,n] = noise(EbtoN0,Eb,Fs,Ps_in,BT,OSF,N);
            % Ajout du bruit au signal en bande passante
            r = s+n;
            er = enveloppeRX(s,Fc,t,OSF,N,PB,order);
            ern = enveloppeRX(r,Fc,t,OSF,N,PB,order);

            demMsgFM = demodulationFM(er,kf,Fs,BT);
            demMsgFMnoise = demodulationFM(ern,kf,Fs,BT);

            FOM_SUM = FOM_SUM + fom(demMsgFMnoise,demMsgFM,N,OSF,SNRin);

            bits_rx_FM = interpretationFM(demMsgFMnoise,N,OSF);
            [bits_rx_FSK,bits_rx_FSK_coh] = demodulationFSK(ern,kf,t,N,OSF);

            BER_FM_SUM = BER_FM_SUM + ber(bits_tx,bits_rx_FM,N);
            BER_FSK_SUM = BER_FSK_SUM + ber(bits_tx,bits_rx_FSK,N);
            BER_FSK_COH_SUM=BER_FSK_COH_SUM+ ber(bits_tx,bits_rx_FSK_coh,N);
        end
    
        BER_FM(kfiter,iter)=BER_FM_SUM/repetition;
        BER_FSK(kfiter,iter)=BER_FSK_SUM/repetition;
        BER_FSK_COH(kfiter,iter)=BER_FSK_COH_SUM/repetition;
        FOM(kfiter,iter)=FOM_SUM/repetition;
        iter = iter+1;
        
        if EbtoN0==10 && kf == 500 
            [psdm,f]=welch(m,t,1024,500);
            [psds,f]=welch(s,t,1024,500);
            [psdes,f]=welch(es,t,1024,500);
            [psdr,f] = welch(r,t,1024,500);
            [psdern,f] = welch(ern,t,1024,500);
            figure("Name","PSD Emetteur")
            plot(f,10*log10(psdm));hold on;
            plot(f,10*log10(psds));hold on;
            plot(f,10*log10(psdes));hold off;
            legend("PSD m(t)","PSD s(t)","PSD e_s(t)")
            xlabel("f[Hz]")
            ylabel("PSD[dB]")

            figure("Name","Comparaison des PSD de s(t) et r(t)");
            plot(f,10*log10(psdr),"r");hold on;
            plot(f,10*log10(psds),"b:");hold off;
            ylabel("PSD[dB]")
            xlabel("f[Hz]")
            legend("PSD r(t)","PSD s(t)")
            
            figure('Name',"PSD Recepteur")
            plot(f,10*log10(abs(h)),"color","[0, 0.75, 0.75]");hold on;
            plot(f,10*log10(psdern),'k--');hold off;
            ylabel("PSD[dB]")
            xlabel("f[Hz]")
            legend("H(f) FPB","PSD e_r(t)")
        end   
    end  
kfiter = kfiter +1;
end

figure("Name","Impact de k_f sur le BER");
plot(scale,(BER_FM(1,:)));hold on;
plot(scale,(BER_FM(2,:)));hold on;
plot(scale,(BER_FM(3,:)));hold on;
legend("k_f=400Hz","k_f=500Hz","k_f=600Hz");
xlabel("Eb/N0 [dB]")
ylabel("BER")

figure("Name","Comparaison des BER");
plot(scale,10*log10(BER_FM(2,:)),'r');
hold on;
plot(scale,10*log10(BER_FSK(2,:)),'b');
hold on;
plot(scale,10*log10(BER_FSK_COH(2,:)),'m');
xlabel("Eb/N0 [dB]");
ylabel("BER [dB]");
hold off;
legend("FM","FSK non-coherent","FSK coherent");

figure("Name","Impact de k_f sur la FOM");
plot(scale,FOM(1,:));hold on;
plot(scale,FOM(2,:));hold on;
plot(scale,FOM(3,:));hold on;
xlabel("Eb/N0 [dB]")
ylabel("FoM [dB]")
legend("k_f=400Hz","k_f=500Hz","k_f=600Hz");
hold off;

%===== FONCTIONS =====%

function er = enveloppeRX(r,Fc,t,OSF,N,PB,order)

    r1 = r.*cos(2*pi*Fc*t);
    r2 = r.*sin(2*pi*Fc*t);
    er1=conv(PB,r1);
    er2=conv(PB,r2);
        % NB : produit de convolution discret => le nombre de samples vaut N +
        % ordre du filtre
    er1cut= er1(order/2+1:OSF*N+order/2);
    er2cut= er2(order/2+1:OSF*N+order/2);
    er=2*(er1cut+er2cut*1j);
        % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
        % multiplier par 2

end

function demMsg = demodulationFM(er,kf,Fs,BT)
    % Demodulation analogique
    % Discriminateur de fréquence
    
    Ac = 1; % Amplitude du signal en bande de base
    a=1;    % Facteur de normalisation -> pente = 1/(4*pi*a)
    
    s1=a*([er(2:end)-er(1:end-1) 0]*Fs+1j*pi*BT*er(1:end));
    s2=-a*([er(2:end)-er(1:end-1) 0]*Fs-1j*pi*BT*er(1:end));

    demMsg = (abs(s1)-abs(s2));    
    
    %NB : le BER a une pente positive si demMsg est divise
    %     par 4*pi*kf*a*Ac comme indique dans les slides
    %     => A discuter 
end

function bits_rx = interpretationFM(demMsg,N,OSF)
    
    % Recherche du bit majoritaire sur les segments OSF

    bits_rx = zeros(1,N);
    for i = 1:N
        SUM=0;
        for j=1:OSF
            SUM=SUM+demMsg(OSF*(i-1)+j);
        end
        if (SUM/OSF) >= 0.5
            bits_rx(i) = 1;
        end
    end

end

function [bits_rx,bits_rx_coh] = demodulationFSK(er,kf,t,N,OSF)
	% Hypotheses
        s0 = exp(-1j*2*pi*kf*t);     % Bits 0
        s1 = exp(1j*2*pi*kf*t);      % Bits 1 
        
    % Correlation de l'enveloppe complexe avec les hypotheses
    % le long des segments OSF
        es0 = er.*conj(s0);
        es1 = er.*conj(s1);

        bits_rx = zeros(1, N);      % FSK non-coherent
        for k = 1 : N 
                sum_0 = abs(sum(es0(((k-1)*OSF)+1:OSF*k)));
                sum_1 = abs(sum(es1(((k-1)*OSF)+1:OSF*k)));
                if sum_0 < sum_1
                   bits_rx(k) = 1;
                end
        end
        
        bits_rx_coh = zeros(1,N);   % FSK coherent
        for k = 1 : N 
                sum_0 = real(sum(es0(((k-1)*OSF)+1:OSF*k)));
                sum_1 = real(sum(es1(((k-1)*OSF)+1:OSF*k)));
                if sum_0 < sum_1
                   bits_rx_coh(k) = 1;
                end
        end
end

function [N0,Pn_in,SNRin,n] = noise(EbtoN0,Eb,Fs,Ps_in,BT,OSF,N)
    % PSD du bruit
        N0=Eb/(10^(EbtoN0/10));
	% Puissance du bruit Gaussien
        Pn_in = N0*Fs/2;
	%   Generation du signal de bruit
        n=sqrt(Pn_in)*randn(1,OSF*N);
    % Calcul du SNR a l'entree en dB
        SNRin = 10*log10(Ps_in/(N0*BT));
end

function FOM = fom(demMsgNoise,demMsg,N,OSF,SNRin)
    % Calcul du bruit pur
        pure_noise = demMsgNoise-demMsg;
    % Puissance du signal bruite en sortie
        Ps_out = trapz(demMsgNoise.^2)/(N*OSF);
    % Puissance du bruit pur
        % NB : demodulation non-lineraire => utilisation du bruit pur
        Pn_out = trapz(pure_noise.^2)/(N*OSF);
    % Calcul du SNR en sortie en dB
        SNRout = 10*log10(Ps_out/Pn_out);
    % Calcul de la FOM en dB
        FOM = SNRout-SNRin;
end

function BER = ber(bits_tx,bits_rx,N)
    % Bit Error Rate
    BER=sum(abs(bits_tx-bits_rx))/N;
end

