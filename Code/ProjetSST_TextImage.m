clear;
clc;
close all;
slCharacterEncoding('UTF-8')

%-----Text to binary-----%

text_tx =  fileread('text.txt');
text_tx_int = uint8(double(text_tx)) ; % 8 bits/character
text_tx_bin = de2bi(text_tx_int,8) ; % decimal to binary
nbchar = size(text_tx_bin,1) ;
bits_text = double(reshape(text_tx_bin',1,nbchar*8)) ;

text_msg = 2*bits_text-1;
%----Image to binary----%

image_tx = 'image.jpg' ; % file.jpg
image_tx_int = rgb2gray(uint8(imread(image_tx,'jpeg'))) ; % gray, 8 bits/pixel
image_tx_bin = (image_tx_int>130) ; % black and white, 1 bit/pixel
[nbl,nbc] = size(image_tx_bin) ;
bits_tx = reshape(image_tx_bin,1,nbl*nbc) ;
bits_img = bits_tx*2-1;

img_msg = 2*bits_img-1;
%------Parametres------%
N=max(length(img_msg),length(text_msg));      
F=500;     % Frequence symbole
T=1/F;     % Periode symbole
Ttot=N/F;  % Periode totale
Fc1=6000;  % Frequence porteuse 1
Fc2=3000;  %
            % ATTENTION au calcul BT => largeur bande 2kHz donc bien
            % espacer
Fs=20000;  % Frequence d'echantillonnage
OSF = Fs*T;% Oversampling factor
ts=1/Fs;   % Periode d'echantillonnage
kf=500;    % Selectivite frequentielle
noiseratio = 50;   %Rapport Eb/N0 en dB
%-----Creation du message avec oversampling-----%
       %On a Nb echantillons par bits

m_text = zeros(1,OSF*N);
for i=1:length(text_msg)                      
    for t=1:OSF                
        m_text(OSF*(i-1)+t)=text_msg(i);
    end
end

m_img = zeros(1,OSF*N);
for i=1:length(img_msg)                   
    for t=1:OSF                
        m_img(OSF*(i-1)+t)=img_msg(i);
    end
end
%-----Generation de l'enveloppe complexe-----%

es_text = exp(1j*2*pi*kf*[0 cumsum(m_text(1:end-1))]/Fs);
        % NB : on divise par Fs a cause de l'integrale
es_img = exp(1j*2*pi*kf*[0 cumsum(m_img(1:end-1))]/Fs);

%-----Signal Radio Frequence-----%

t=0:ts:Ttot-ts;
s_text = real(es_text).*cos(2*pi*Fc1*t)+imag(es_text).*sin(2*pi*Fc1*t);
s_img = real(es_img).*cos(2*pi*Fc2*t)+imag(es_img).*sin(2*pi*Fc2*t);
s = s_text+s_img;

[psds,f]=welch(s,t,1024,500);
figure('Name',"PSD Signal")

plot(f,10*log10(psds))

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
                                           
                                           
                                           
%======================= TEXTE ====================================== % 

    %===== SIGNAL BRUITE ====%                                          
    %----Convolution du filtre----%

    r1 = r.*cos(2*pi*Fc1*t);
    r2 = r.*sin(2*pi*Fc1*t);
    er1=conv(PB,r1);
    er2=conv(PB,r2);
        % NB : produit de convolution discret => le nombre de samples vaut N +
        % ordre du filtre
    er1cut= er1(order/2+1:OSF*N+order/2);
    er2cut= er2(order/2+1:OSF*N+order/2);
    ern_text=2*(er1cut+er2cut*j);
        % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
        % multiplier par 2

    %---- Demodulation FM ----%

    Ac = 1; % Amplitude du signal R*(t) en bande de base

    a=1/(Ac*4*pi*kf);

        % Pour faire la derivee on fait la variation avec :
    s1_text=a*([(ern_text(2:end)-ern_text(1:end-1)) 0]*Fs+j*pi*BT*ern_text(1:end));
    s2_text=-a*([(ern_text(2:end)-ern_text(1:end-1)) 0]*Fs-j*pi*BT*ern_text(1:end));

    mrn = abs(s1_text)-abs(s2_text);    %message receive noise

    %===== SIGNAL SANS BRUIT ====%                                          

    r1 = s.*cos(2*pi*Fc1*t);
    r2 = s.*sin(2*pi*Fc1*t);
    er1=conv(PB,r1);
    er2=conv(PB,r2);
    er1cut= er1(order/2+1:OSF*N+order/2);
    er2cut= er2(order/2+1:OSF*N+order/2);
    er_text=2*(er1cut+er2cut*j);
    
    s1_text=a*([(er_text(2:end)-er_text(1:end-1)) 0]*Fs+j*pi*BT*er_text(1:end));
    s2_text=-a*([(er_text(2:end)-er_text(1:end-1)) 0]*Fs-j*pi*BT*er_text(1:end));

    mr = abs(s1_text)-abs(s2_text);    %message receive noise


    %===== Calcul du SNR =====%

    SNRin = 10*log10(Ps_in/(N0*BT));

    pure_noise = mrn-mr;    %bruit pure
    Ps_out = sum(mrn.^2)/(N*OSF);
    Pn_out = sum(pure_noise.^2)/(N*OSF);

    SNRout = 10*log10(Ps_out /Pn_out);

        % Facteur de mérite 
    FoM_text = SNRout-SNRin

        %Pour ameliorer le facteur de merite il faut filtrer le signal de 
        %facon a enlever les parasites qui apparaissent dans la PSD du 
        %signal recu

        %NB : il faut faire un plot Fom en fonction de Eb/N0 et pour chaque
        %point il faut faire une moyenne sur 100 FoM afin d'avoir la courbe la
        %plus précise
        
    %==== Demodulation FSK ====%

    s0 = exp(-1j*2*pi*kf*t);     
    s1 = exp(1j*2*pi*kf*t);      

    es0 = ern_text.*conj(s0);
    es1 = ern_text.*conj(s1);

    %Correlation
    mf2 = zeros(1, N);
    for k = 1 : N 
            sum_0 = abs(sum(es0(((k-1)*OSF)+1:OSF*k)));
            sum_1 = abs(sum(es1(((k-1)*OSF)+1:OSF*k)));
            if sum_0 < sum_1
                mf2(k) = 1;
            end
    end


    text_rx_bin = reshape(mf2(1:length(text_msg)),8,nbchar)' ; % nbchar : number of characters
    text_rx_int = bi2de(text_rx_bin)' ; % binary to decimal
    text_rx = char(text_rx_int)  % decimal to text

%======================= IMAGE ====================================== % 

    %===== SIGNAL BRUITE ====%                                          
    %----Convolution du filtre----%

    r1 = r.*cos(2*pi*Fc2*t);
    r2 = r.*sin(2*pi*Fc2*t);
    er1=conv(PB,r1);
    er2=conv(PB,r2);
        % NB : produit de convolution discret => le nombre de samples vaut N +
        % ordre du filtre
    er1cut= er1(order/2+1:OSF*N+order/2);
    er2cut= er2(order/2+1:OSF*N+order/2);
    ern_img=2*(er1cut+er2cut*j);
        % Faire les calculs avec les sinus pour demontrer le fait qu'il faille
        % multiplier par 2

     xlabel("f[Hz]")

    %---- Demodulation FM ----%

    Ac = 1; % Amplitude du signal R*(t) en bande de base

    a=1/(Ac*4*pi*kf);

        % Pour faire la derivee on fait la variation avec :
    s1_text=a*([(ern_img(2:end)-ern_img(1:end-1)) 0]*Fs+j*pi*BT*ern_img(1:end));
    s2_text=-a*([(ern_img(2:end)-ern_img(1:end-1)) 0]*Fs-j*pi*BT*ern_img(1:end));

    mrn = abs(s1_text)-abs(s2_text);    %message receive noise

    %===== SIGNAL SANS BRUIT ====%                                          
    %----Convolution du filtre----%

    r1 = s.*cos(2*pi*Fc1*t);
    r2 = s.*sin(2*pi*Fc1*t);
    er1=conv(PB,r1);
    er2=conv(PB,r2);
    er1cut= er1(order/2+1:OSF*N+order/2);
    er2cut= er2(order/2+1:OSF*N+order/2);
    er_img=2*(er1cut+er2cut*j);
    
    s1_img=a*([(er_img(2:end)-er_img(1:end-1)) 0]*Fs+j*pi*BT*er_img(1:end));
    s2_img=-a*([(er_img(2:end)-er_img(1:end-1)) 0]*Fs-j*pi*BT*er_img(1:end));

    mr = abs(s1_img)-abs(s2_img);    %message receive noise

    %===== Calcul du SNR =====%

    SNRin = 10*log10(Ps_in/(N0*BT));

    pure_noise = mrn-mr;    %bruit pure
    Ps_out = sum(mrn.^2)/(N*OSF);
    Pn_out = sum(pure_noise.^2)/(N*OSF);

    SNRout = 10*log10(Ps_out /Pn_out);

        % Facteur de mérite 
    FoM_img = SNRout-SNRin

        %Pour ameliorer le facteur de merite il faut filtrer le signal de 
        %facon a enlever les parasites qui apparaissent dans la PSD du 
        %signal recu

        %NB : il faut faire un plot Fom en fonction de Eb/N0 et pour chaque
        %point il faut faire une moyenne sur 100 FoM afin d'avoir la courbe la
        %plus précise
        
    %==== Demodulation FSK ====%

    s0 = exp(-1j*2*pi*kf*t);     
    s1 = exp(1j*2*pi*kf*t);      

    es0 = ern_img.*conj(s0);
    es1 = ern_img.*conj(s1);

    %Correlation
    mf2 = zeros(1, N);
    for k = 1 : N 
            sum_0 = abs(sum(es0(((k-1)*OSF)+1:OSF*k)));
            sum_1 = abs(sum(es1(((k-1)*OSF)+1:OSF*k)));
            if sum_0 < sum_1
                mf2(k) = 1;
            end
    end
    
    figure;
    subplot(2,1,1);
    image_tx = image_tx_bin*255;
    colormap gray ; image(image_tx) ; % plot result
    subplot(2,1,2);
    image_rx_bin = reshape(mf2,nbl,nbc) ; % nbl : number of lines, nbc : number of columns
    image_rx = image_rx_bin*255 ; % white : 255, black : 0
    colormap gray ; image(image_rx) ; % plot result
