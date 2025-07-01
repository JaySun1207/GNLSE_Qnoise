% task1:   Analytically derive NLSE (keep only 2nd order dispersion)
% task2/3: Critical power (only 2nd order, with & without raman, no self-steepening)
% task4:   Soliton self-freq shift (red-shift) 


% Solves normalized NLS equation with split-step method
% Modifed Raman response function (2006 OL paper) 

clear; tic; 
close all; 

% When lam0 changes, the 0 freq crossing (if there is one) also changes ...
% 



%---Specify input parameters
T0   = 50; % 50;        % pulse width in fs (1e-15 sec)
lam0 = 810;       % in nm
s = lam0/(2*pi*2.998e2*T0);  % self-steepening; 3e2 b/c T0(fs), lam0(nm), 3e8(m/s); s seems to be gamma_1 ? 

chirp0 = -290/(T0^2); %0;     % input pulse chirp (default value; fs^2)
sbeta2 = -1;    % sign of beta_2 (+1 or -1); beta2 normalized to 1
delta3 = 0.0;     % delta3= 0.05;   % beta3/(6*T0*abs(beta2)) 
delta4 = 0.0;   % beta4/(24*T0^2*abs(beta2))
%delta5=0; delta6=0; delta7=0; delta8=0;

distance = 20; % 1;   % normalized to dispersion length
N = 2;          % soliton order
mshape = 0;     % 0 for sech, m>0 for super-Gaussian



raman_on           = 1;    % 0 or 1 
self_steepening_on = 1;    % 0 or 1 
temporal_quantum_noise_on = 0; 
frequency_quantum_noise_on = 1; 
power_adjustment   = sqrt(1); % sqrt(3.5);   

number_of_trials    = 10; % 30; 


%---set silumation parameters
nt       = 2^12;                % N in Agrwal. 
Tmax     = 20; % 100;          % Tm in Agrwal; FFT points and window size; t/T0
step_num = 1000; % 500;          % round(40*distance*N^2); 
zstep    = 100;                 % # of z steps
ng       = step_num/zstep;      % store data every 'ng' steps to 3D 
deltaz   = distance/step_num;   % step size in z
dtau     = (2*Tmax)/nt;         % \delta t in Agrwal; step size in tau 
domega   = (pi/Tmax);           % \delta \omega in Agrwal; 2pi/Tmax % 
% domega   = (1/Tmax);

%-------tau and omega arrays
% When nt = 2048, (-nt/2:nt/2-1 ) is -1024 to +1023
tau = (-nt/2:nt/2-1)*dtau;                 % time array; *T0_SI to get SI 
omega = fftshift((-nt/2:nt/2-1)*domega );  % freq array; 'fftshifted'; /T0_SI to get SI 

% Defination of Raman response function h(t) 
fR = 0.18; % 0.245;  % Raman fraction 
fb = 0.21;   % Boson peak fraction 
tau1 =  12.2/T0; 
tau2 =  32/T0; 
tau_b = 96/T0;
h = (1-fb)*(tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-tau/tau2).*sin(tau/tau1)...
      + fb*((2*tau_b-tau)/tau_b^2).*exp(-tau/tau_b);
h(1:nt/2) = 0;      % causality

% figure() 
% plot( h ); 

beta2_SI = 1.183e-26; % unit: s^2/m 
gamma_SI = 0.11;   % 1e-3;  % unit: 1/(W*m)
A_soliton = sqrt(beta2_SI/gamma_SI) * (1/(T0 * 1e-15));  % tau_sech is T0

L_nl = 1/(gamma_val* (A_soliton)^2);  % nonlinear length / disp length 

% tau is t/T0; dtau = dt * (1/T0) 
center_freq = 2.998e8 / (lam0 * 1e-9);   % Center freq (in Hz) = c/lam0; exact value 

% Temporal quantum noise 
hbar_omega0_by_Dt     = 6.626e-34 * center_freq / (dtau * T0*1e-15 ); 
zeta_factor  = sqrt( hbar_omega0_by_Dt );    % Normalization factor zeta; 



% *** Notation: Omega = omega0 + omega; (omega in [ -omega_max, +omega_max])
% hbar_Omega_by_Domega = 6.626e-34 * abs(center_freq)/(domega/(T0 * 1e-15)); 
hbar_Omega_by_Domega = 6.626e-34 * abs(omega/(T0 * 1e-15 *2*pi) + center_freq)/(pi * (domega/(T0 * 1e-15))); 
zeta_factor_omega = sqrt( 2*pi*hbar_Omega_by_Domega );    % Normalization factor zeta; 


% Plot the range of freqs. 
linear_freq = omega/(T0 * 1e-15 *2*pi) + center_freq; 
wavelength =  2.99792e8*1e9 ./ linear_freq; 
figure() 
plot(omega, linear_freq )
xlabel('Angular Freq (\nu-\nu_0)T_0'); ylabel('Linear Freq (Hz)');
title('Angular freq \nu vs. Linear freq (Hz)'); 
figure() 
plot( fftshift(omega), fftshift(wavelength) )
xlabel('Angular Freq (\nu-\nu_0)T_0'); ylabel('Wavelength (nm)');
title('Angular freq \nu vs. Wavelength (nm)'); 





%----store dispersive phase shifts to speedup code
dispersion = exp(1i*deltaz*(0.5*sbeta2*omega.^2 +...
    delta3*omega.^3 + delta4*omega.^4));
hhz = 1i*N^2*deltaz;    % similar to 'h' in Agrwal





uu_end              = zeros(number_of_trials, nt ); 
spect_end           = zeros(number_of_trials, nt ); 
A_tilde_spect_end   = zeros(number_of_trials, nt ); 
noisy_uu_input_all  = zeros(number_of_trials, nt ); 


% Parsevals Thm?
% uu = sech(tau); 
% figure() 
% plot( tau, abs(uu).^2); 
% trapz( tau, abs(uu).^2 );           % = 2 
% sum(abs(uu.^2))
% 
% figure() 
% plot( omega, abs(ifft(uu)).^2 ); 
% trapz( omega, abs(ifft(uu)).^2 );   % = 3.0642e-04 = 2 / (2^11 * pi ) 
% sum(abs(ifft(uu).^2))
% 
% trapz(  omega/pi , abs(ifft(uu)).^2 );   % = 3.0642e-04 = 2 / (2^11 * pi ) 


ifft_norm_factor = nt*sqrt( 2*pi*dtau*T0*1e-15/(domega/(T0*1e-15)) ); 



for trial_index = 1:number_of_trials
    if( mod(trial_index, 10) == 0  )
        disp( trial_index )
    end 

    %---Input Field profile
    if mshape == 0
        uu_ini = power_adjustment * sech(tau).*exp(-1i*chirp0*tau.^2);     %soliton; power scaling here
    else
        uu_ini = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape)); %super-Gaussian
    end
    uu = uu_ini; 

    
    % If want to see profile 
    % figure() 
    % plot( abs(uu) ); 
    % hold on 
    % plot( real(uu) ); 
    % plot( imag(uu) ); 
    % hold off 
    
    %-----Input temporal noise 
    if temporal_quantum_noise_on == 1 
        gaussian_noise = sqrt( hbar_omega0_by_Dt/2 )/A_soliton * randn(size(uu) ,like=1i); 
        uu = uu + gaussian_noise; 
        noisy_uu_input_all(trial_index, :) = uu; 
    end 

    
    % These 2 are the same 
    % uu_replica1 = fft( ifft(uu)); 
    % uu_replica2 = fft( ifftshift(fftshift(ifft(uu)))); 
    

    %-----Input freq noise 
    if frequency_quantum_noise_on == 1 
        freq_noise = (sqrt( pi * hbar_Omega_by_Domega )/A_soliton) .* randn( size(uu) ,like=1i); 
        % actual_fourier_spectrum = fft(uu) * (dtau*T0*1e-15); 
        actual_fourier_spectrum = ifft(uu) * (dtau*T0*1e-15) * nt; 
        actual_fourier_spectrum = ifftshift(fftshift(actual_fourier_spectrum) + fftshift(freq_noise)); 
        % uu = fft(actual_fourier_spectrum / (dtau*T0*1e-15) ) ; 
        uu = fft(actual_fourier_spectrum ) * (domega/(T0*1e-15)) /(2*pi) ; 
        noisy_uu_input_all(trial_index, :) = uu;         
    end 
    
    if( trial_index == 1 )
        figure() 
        plot( tau, abs(uu_ini),  tau, abs(uu) ); 
        xlabel('Tau (t/T_0)'); ylabel('uu absolute val ');
        title('uu initial with and without noise'); 
        legend('without noise','with noise')
    end 
    
    %-----Input 3D data   
    uu_3d(1,:)=uu.*conj(uu);      % Intensity in time domain
    temp = fftshift(ifft(uu));    % use ifft to do FT 
    spect_3d(1,:)=abs(temp).^2;   % spectrum au unit
    

    % figure() 
    % plot(omega/T0, abs(2*Tmax*T0*1e-15* ifft(uu)) ); 
    % figure()
    % plot(omega/T0, sqrt( pi * hbar_omega_by_Domega )/A_soliton )

    
    % scheme: 1/2N -> D -> 1/2N; first half step nonlinear
    temp = uu.*exp(abs(uu).^2.*hhz/2);    % note hhz/2; This is N_hat*h/2 

    %*********[ Beginning of MAIN Loop]***********
    for n = 1:step_num 
        
        % This step is D_hat*h
        f_temp = ifft(temp).*dispersion;    % Transform to Fourier to apply dispersion; D(omega) 
        uu   = fft(f_temp);                 % Inverse Transform to time domain 
        sst1 = deriv1(abs(uu).^2,dtau); 
        sst2 = deriv1(uu,dtau) ; 
        sst  = 1i*(1-fR)*s*(sst1 + conj(uu).*sst2) ; 

        % Eq. 2.3.36,  Calculate the convolution of 'R' and '|A|^2'
        P      = uu.*conj(uu);   % ??? include A_soliton ??? 
        convl  = (nt*dtau)*fft(ifft(h).*ifft(P));    % Time domain convolution <--> Freq domain multplication
        convl  = fftshift(convl);
        sst3   = deriv1(convl,dtau) ; % off
        sstnew = 1i*s*fR*(sst3 + (conj(uu).*convl)./(P+eps).*sst2);  % s <--> gamma1; gamma absorbed in disp length

        % This part does N_hat*h; N_hat includes raman and self-steepening 
        if( self_steepening_on == 1 && raman_on == 1 )
            temp = uu.*exp(((1-fR)*(abs(uu).^2)+ sst + fR*convl + sstnew).*hhz);   %original
        elseif( self_steepening_on == 1 && raman_on == 0 )
            fR = 0;  % need to move this 
            temp = uu.*exp(((1-fR)*(abs(uu).^2)+ sst + fR*convl + sstnew).*hhz);  % sst only; no raman
        elseif( self_steepening_on == 0 && raman_on == 1 )
            temp = uu.*exp(((1-fR)*(abs(uu).^2)+ sst*0 + fR*convl + sstnew*0).*hhz);  % Raman only; no sst 
        else
            temp = uu.*exp((abs(uu).^2) .*hhz);  % no sst ; no Raman 
        end
        % temp   = uu.*exp((  (1)*(abs(uu).^2) + sst + sstnew).*hhz); 
        spect  = fftshift(ifft(temp));

        % take sample every ng steps 
        if mod(n,ng) == 0
             uu_3d((n/ng + 1),:)     =  temp.*conj(temp);    % temporal Intensity
             spect_3d((n/ng + 1),:)  = spect.*conj(spect);   % spectral Intensity
        end
    end
    %***************[ End of MAIN Loop ]**************

    temp = temp.*exp(-1*abs(temp).^2.*hhz/2); % note hhz/2; This is -1*N_hat*h/2 
    spect  = fftshift(ifft(temp)); 
    uu_end(trial_index, :)    = temp; 
    spect_end(trial_index, :) = spect; 
    % A_tilde_spect_end(trial_index, :) = spect ./ zeta_factor_omega;    

    %----Plot output pulse shape and spectrum 
    % colormap(1-jet)
    % z = linspace(0,distance,zstep+1); % normalised distance
    % subplot(1,2,1);
    %     uu_3d = 10*log10(uu_3d);      % convert to dB units
    %     pmax = max(max(uu_3d));       % maximum value
    %     pcolor(tau,z,uu_3d);          % 3D plot
    %     caxis([pmax-50, pmax]); xlim([-50,50]); shading interp;
    %     set(gca,'FontSize', 12); %title('(a)')
    %     xlabel('Time (t/T_0)'); ylabel('Distance (z/L_D)')
    % 
    % subplot(1,2,2)
    %     spect_3d=10*log10(spect_3d);    % convert to dB units
    %     % spect_3d=spect_3d;
    %     freq = fftshift(omega/(2*pi));  % freq array
    %     pmax = max(max(spect_3d));      % max value for scaling plot
    %     pcolor(freq,z,spect_3d);
    %     caxis([pmax-50, pmax]); xlim([-4,4]); shading interp;
    %     set(gca,'FontSize', 12); %title('(a)')
    %     xlabel('(\nu-\nu_0)T_0'); ylabel('Distance (z/L_D)')
 
end 



% A^2 in units of A_soliton 
% uu_temporal_intensity_mean = mean( abs(uu_end).^2 , 1);
% uu_temporal_intensity_SD = std( abs(uu_end).^2 , 1); 
% noisy_uu_input_intensity_mean = mean( abs(noisy_uu_input_all).^2, 1); 
% noisy_uu_input_intensity_SD = std( abs(noisy_uu_input_all).^2, 1); 
% figure() 
% plot( tau, uu_temporal_intensity_SD, tau, noisy_uu_input_intensity_SD)
% xlabel('Time (t/T_0)'); ylabel('A_{soliton} (sqrt(Watt))');
% title('Standard Deviation of I/O intensities - A_{soliton} unit'); 
% figure() 
% plot( tau, uu_temporal_intensity_mean, tau, noisy_uu_input_intensity_mean)
% xlabel('Time (t/T_0)'); ylabel('A_{soliton} (sqrt(Watt))');
% title('mean of I/O intensities - A_{soliton} unit'); 



% ----- Calculate n(t) initial and final 
n_t_initial = (dtau*T0*1e-15*(abs(noisy_uu_input_all * A_soliton).^2) )/(6.626e-34 * center_freq);  
n_t_final   = (dtau*T0*1e-15*(abs(uu_end * A_soliton ).^2) )/(6.626e-34 * center_freq); 

n_t_final_SD        =  std( n_t_final , 1); 
n_t_final_mean      = mean( n_t_final , 1); 
n_t_initial_SD   =  std( n_t_initial , 1); 
n_t_initial_mean = mean( n_t_initial , 1); 

figure() 
plot( tau, n_t_final_SD .^ 2, tau, n_t_initial_SD .^2)
xlabel('Time (t/T_0)'); ylabel('photon count n(t)');
title('Variance of I/O intensities; temporal; photon unit'); 
legend('Output Intensity var.','Input Intensity var.')
figure() 
plot( tau, n_t_final_mean, tau, n_t_initial_mean)
xlabel('Time (t/T_0)'); ylabel('photon count n(t)');
title('Mean of I/O intensities; temporal; photon unit');
legend('Output Intensity mean','Input Intensity mean')
figure() 
plot( tau, n_t_initial_SD .^2, tau, n_t_initial_mean )  % 2 pi??
xlabel('Time (t/T_0)'); ylabel('photon count n(t)');
title('Input Time Domain mean vs. Input Time Domain variance'); 
legend('Input Intensity var.','Input Intensity mean')


% ----- Calculate n(omega) initial and final 
A_omega_initial = dtau*T0*1e-15*nt*abs(ifft(noisy_uu_input_all * A_soliton, nt, 2)); 
n_omega_initial = (A_omega_initial ./ zeta_factor_omega).^2;
n_omega_initial_mean = mean( n_omega_initial, 1); 
n_omega_initial_SD   =  std( n_omega_initial, 1); 

A_omega_final = dtau*T0*1e-15*nt*abs(ifft(uu_end * A_soliton, nt, 2)); 
n_omega_final = (A_omega_final ./ zeta_factor_omega).^2;
n_omega_final_mean = mean( n_omega_final, 1); 
n_omega_final_SD   =  std( n_omega_final, 1); 


% 
% figure() 
% plot(fftshift(omega), fftshift(n_t_initial_SD.^2), fftshift(omega), fftshift(n_t_final_SD.^2) )
% xlabel('Freq (\nu-\nu_0)T_0'); ylabel('photon count n(\nu) ');
% title('Variance of I/O intensities; spectral; photon unit'); 
% legend('Input Intensity var.','Output Intensity variance')
% figure() 
% plot(fftshift(omega), fftshift(n_t_initial_mean), fftshift(omega), fftshift(n_t_final_mean) )
% xlabel('Freq (\nu-\nu_0)T_0'); ylabel('photon count n(\nu) ');
% title('Mean of I/O intensities; spectral; photon unit'); 
% legend('Input Intensity mean','Output Intensity mean')
figure() 
plot(fftshift(omega), fftshift(n_omega_initial.^2), fftshift(omega), fftshift(n_omega_final_mean) )
xlabel('Freq (\nu-\nu_0)T_0'); ylabel('photon count n(\nu) ');
% xlabel('Freq ( (\omega - \omega_0)T_0 )'); ylabel('photon count n(\omega) ');
title('Input Freq Domain mean vs. Input Freq Domain variance'); 
legend('Input Intensity var.','Input Intensity mean')



figure() 
plot(fftshift(omega), log10(fftshift(n_omega_initial_SD .^2)), ... 
    fftshift(omega), log10(fftshift(n_omega_initial_mean)), ...
    fftshift(omega), log10(fftshift(n_omega_final_SD.^2)), ...
    fftshift(omega), log10(fftshift(n_omega_final_mean)) )
xlabel('Freq (\nu-\nu_0)T_0'); 
ylabel('photon count log10(n(\nu)) ');
title('Variance & Mean of I/O intensities; spectral; photon unit'); 
legend('Input Intensity variance', 'Input Intensity mean', ... 
    'Output Intensity variance', 'Output Intensity mean')


% zeta_factor = 1 ; % overwrite for now
% A_tilde_temp_end  = uu_end / zeta_factor;   % zeta factor ~ 1e-2
% n_temporal_mean = mean( abs(A_tilde_temp_end).^2, 1); 
% n_temporal_SD   =  std( abs(A_tilde_temp_end).^2, 1); 
% n_spectral_mean = mean( abs(A_tilde_spect_end).^2, 1); 
% n_spectral_SD   =  std( abs(A_tilde_spect_end).^2, 1); 
% phi_spectral_mean = mean( angle(A_tilde_spect_end).^2, 1); 
% phi_spectral_SD   =  std( angle(A_tilde_spect_end).^2, 1); 
% figure() 
% plot(tau, abs( sqrt(uu_3d(1,:))/zeta_factor ).^2 , 'DisplayName', "Initial profile" ); 
% hold on
% errorbar(tau(1:10:end), n_temporal_mean(1:10:end), n_temporal_SD(1:10:end)); 
% xlabel('Time (t/T_0)'); ylabel('photon count n(t, t+\Deltat) ');  % squared
% legend()
% hold off 
% If want to plot every trial to see the spread 
% for trial_index = 1:number_of_trials
%     plot(tau, abs(uu_end(trial_index, :))); 
% end 
% errorbar(tau(1:20:end), uu_temporal_mean(1:20:end), uu_temporal_SD(1:20:end)); 
% xlabel('Time (t/T_0)'); ylabel('abs( uu )'); 
% figure() 
% plot(omega, abs( fftshift(sqrt(spect_3d(1,:)) ./ zeta_factor_omega) ).^2, 'DisplayName', "Initial profile"); 
% hold on
% errorbar(omega(1:10:end), fftshift(n_spectral_mean(1:10:end)), fftshift(n_spectral_SD(1:10:end))); 
% xlabel('(\nu-\nu_0)T_0'); ylabel('photon count n(\nu) ');  % squared
% legend()
% hold off 
% figure() 
% % plot(omega, abs( fftshift(sqrt(spect_3d(1,:)) ./ zeta_factor_omega) ).^2, 'DisplayName', "Initial profile"); 
% % hold on
% errorbar(omega(1:10:end), fftshift(phi_spectral_mean(1:10:end)), fftshift(phi_spectral_SD(1:10:end))); 
% xlabel('(\nu-\nu_0)T_0'); ylabel('phi (\nu, \nu+\Delta\nu) ');  % squared
% legend()


toc
