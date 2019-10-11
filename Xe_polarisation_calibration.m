%% Xe polarisation calibration program
clear all

%% Xe and Water spectra
for a = 1:2
    
Xp = 'D:\Uni stuff\MSci work\Joe and Jake\OPcellcoil_Xe_optimisation_2\6';
Hp = 'D:\Uni stuff\MSci work\Joe and Jake\cellcoil_H_1.39A\1';
%These are the paths for the Xe and H files, change as necessary

    if a == 1
        addpath(Xp) %selects Xe data
    elseif a == 2
        addpath(Hp) %selects water data
    end
    
    Hdata = dlmread('data.csv'); %reads Prospa data
    Hfreq = Hdata(:,1); Hmag = Hdata(:,2);
    
    
    if size(Hdata) == [8192,1] %Data saved as single vector
        Hfreq = 1:8192;
        Hmag = Hdata(:,1);
    elseif size(Hdata) == [8192,2] %Data saved as xy data (magnitude 
        %spectrum)
        Hfreq = Hdata(:,1); Hmag = Hdata(:,2);
    elseif size(Hdata) == [8192,3] %Data from automatic saving in Prospa 
        %(FID)
        
        %Define parameters of the scan (can be found in acqu.txt)
        Hbandwidth = 20; %Bandwidth (kHz)
        
        Hpoints = 8192; %No. of points
        Hinterval = Hbandwidth/Hpoints; %Freq. interval (kHz)
        T_H = Hpoints/Hbandwidth; %Sample time (ms)
        L_H = length(Hdata); %Length of signal
        t_H = ((0:L_H-1)/Hbandwidth)'; %Time vector (calculated) (ms)
        Htime = (Hdata(:,1)); %Time vector (from data) (ms)
        
        %Create frequency array
        Hfreq=(linspace(-(Hbandwidth/2)+Hinterval, (Hbandwidth/2), 8192))';
        
        %Process the raw FID data
        %convert columns 2 and 3 into complex form
        c_H=complex(Hdata(:,2),Hdata(:,3));
        
        %Perform FFT, shift low freq components to centre of k-space, flip 
        %spectrum in y-axis
        Hspectrum=fftshift(fft(flipud(c_H)));
        
        %Take the absolute value of the data ie. sqrt((real(spectrum)).^2+
        %(imag(spectrum)).^2); and divide by bandwidth to get amplitude in 
        %uV/kHz
        Hmag=(abs(Hspectrum))./Hbandwidth;
        
        % %Plot full magnitude spectrum
        % figure; plot(freq,mag)
        % title('Magnitude spectrum')
        % xlabel('Frequency(kHz)'); ylabel('Amplitude (\muV/kHz)')
        % xlim([min(freq) max(freq)])
    end
        
    Hmag2= Hmag(3900:4300);
    Hfreq2 = Hfreq(3900:4300);
    [val4, idx4] = max(Hmag2);
    Hsize = 150;
    mini = idx4-Hsize;
    maxi = idx4+Hsize;
    
    Hx = Hfreq2(mini:maxi);
    Hy = Hmag2(mini:maxi);
    [val2, idx2] = max(Hy);
    %Data used for calculating PEAK
    
    Hx_err_1 = Hfreq(1000:2000); Hy_err_1 = Hmag(1000:2000);
    Hy_noise=Hy_err_1(900:1000); %Assign noise regions in spectrum
    std_noise=std(Hy_noise); %noise spread in data
    max_signal=max(Hy_err_1);% Maximum peak height in data
    m=mean(Hy_noise);max_signal=max_signal-m; %background correct
    SNR=max_signal/std_noise; % signal to noise ratio
    err=1/SNR;
    
    S_H_1 = trapz(Hx,Hy); %1H magnitude spectrum peak integral (uV)
    S_H_1_V = S_H_1*10^-6;
    S_H_max=max(Hy);
    %H_fit = fit(Hx,Hy,'gauss3');
    %plot(H_fit,Hx,Hy)
    
    n = 0;% Calculate fwhm; find two points where height is half
    
    while Hy(n+1) < S_H_max/2
        n = n + 1;
    end
    
    j = idx2; % account for noise in spectrum dont need to try next few
    %points
    while Hy(j+1) > S_H_max/2
        j = j + 1;
    end
    
    fwhm_H = Hx(j) - Hx(n);
    nu_H_1 = fwhm_H*1000;          %FWHM of magnitude spectrum (Hz)
    t2s_H_1 = (sqrt(3)/pi)*(1/nu_H_1);      %T2* relaxation time (s)
    t2s_H_ms_1 = 1000*t2s_H_1;              %T2* relaxation time (ms)
    
    if a == 1 %saves values for Xe 
        nu_Xe = nu_H_1;
        T2_Xe = t2s_H_1;
        S_Xe = S_H_1_V;
        error_Xe = err;
    elseif a == 2 %saves values for H
        nu_H = nu_H_1;
        T2_H = t2s_H_1;
        S_H = S_H_1_V;
        error_H = err;
    end
    
end
%% Calibration

alphaXe = 75.*pi./180; %flip angle for Xe
alphaH = 77.6.*pi./180; %flip angle for H
P_H_therm = 2.67e-9; %thermal polarisation of water
N_Xe = 1.76e19; %number of Xe atoms
N_H = 6.55e22; %number of H atoms
beta129 = 0.264; %isotopic abundance of 129Xe
gammaXe = 7.441; %gyromganetic ratio of Xe
gammaH = 26.75; %gyromganetic ratio of H
tp_Xe = 350e-6; %pulse length of Xe
tp_H = 100e-6; %pulse length of H

exp_Xe = exp(-tp_Xe./T2_Xe);
exp_H = exp(-tp_H./T2_H);
%calculates exponential terms for both Xe and H

P_Xe = (S_Xe.*P_H_therm.*N_H.*gammaH.*(sin(alphaH)).*exp_H)./...
    (S_H.*N_Xe.*beta129.*gammaXe.*(sin(alphaXe)).*exp_Xe);
%calculates polarisation of Xe

P_Xe_percent = P_Xe.*100; %converts to a percentage

disp(['Xe polarisation = ',num2str(P_Xe_percent),' +/- ',...
    num2str(error_H.*P_Xe_percent),'%']);

