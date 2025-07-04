% task1:   Analytically derive NLSE (keep only 2nd order dispersion)
% task2/3: Critical power (only 2nd order, with & without raman, no self-steepening)
% task4:   Soliton self-freq shift (red-shift) 

% Solves normalized NLS equation with split-step method
% Modifed Raman response function (2006 OL paper) 

clear; tic; 
% close all; 

%---Specify input parameters
T0   = 28.36; % 50;        % pulse width in fs
disp("------------------------------------------")
disp("FWHM Dt = " + num2str(T0 * 1.763) + " (fs)")
lam0 = 835;       % in nm


distance = 2.24; % 1;  % normalized to dispersion length


N = 8.5; %sqrt(25); % sqrt(20); %1;          % soliton order
mshape = 0;     % 0 for sech, m>0 for super-Gaussian
raman_on           = 1;    % 0 or 1 
self_steepening_on = 1;    % 0 or 1 
% temporal_quantum_noise_on = 1; 
% frequency_quantum_noise_on = 0; 
raman_noise_on = 1;   % to be implemented 
power_adjustment   = sqrt(1); %sqrt(2.0); % sqrt(3.5);   



%--- self-steepening
if self_steepening_on == 1
    s = 0.0112; % 0.56fs = 0.0112T0 
    % s = lam0/(2*pi*3e2*T0); % =1/(omega0*T0)=0.0086  --> (norm. by T0) 
else 
    s = 0; 
end 




% Frequency noise only added within range 400 to 1600 nm ?? 

chirp0 =  0; % -290/(T0^2);      % input pulse chirp (default value)
chirp_3rd_order = 0; %-500/(T0^3);  % from paper 
sbeta2 = -1;        % sign of beta_2 (+/-1); beta2 normalized to 1
% delta3 = 0.0;     % delta3= 0.05;   % beta3/(6*T0*abs(beta2)) 
% delta4 = 0.0;     % beta4/(24*T0^2*abs(beta2))
% delta
%delta5=0; delta6=0; delta7=0; delta8=0;


% ********* Dispersion param for Fused Silica: 
% @ 800nm:   beta2 = +35 fs^2/mm =  3.5e-26 s^2/m
% @ 1500nm:  beta2 = -26 fs^2/mm = -2.6e-26 s^2/m
% We are using PCF, beta values above 
% beta2_val = 3.5e-26; % unit: s^2/m 

beta2_SI = 1.183e-26; % unit: s^2/m 
gamma_SI = 0.11;   % 1e-3;  % unit: 1/(W*m)
A_soliton = sqrt(beta2_SI/gamma_SI) * (1/(T0 * 1e-15));  % tau_sech is T0

% ********* Dispersion param for PCF; note unit 
beta_2 = -11.830;      % ps^2/km
beta_3 =  8.1038e-2;   % ps^3/km
beta_4 = -9.5205e-5;   % ps^4/km 
beta_5 =  2.0737e-7;   % ps^5/km 
beta_6 = -5.3943e-10;  % ps^6/km 
beta_7 =  1.3486e-12;  % ps^7/km 
beta_8 = -2.5495e-15;  % ps^8/km 
beta_9 =  3.0524e-18;  % ps^9/km 
beta_10 = -1.7140e-21; % ps^10/km 
% delta3 = 0; 
% delta4 = 0; 
% delta5 = 0; 
% delta6 = 0; 
% delta7 = 0; 
delta3 = beta_3/(factorial(3)*(T0*1e-3)^1*abs(beta_2));  % beta3_SI/(3!* (T0_SI)*abs(beta2_SI))  
delta4 = beta_4/(factorial(4)*(T0*1e-3)^2*abs(beta_2));  % beta4_SI/(4!* (T0_SI)^2 *abs(beta2_SI))  
delta5 = beta_5/(factorial(5)*(T0*1e-3)^3*abs(beta_2));
delta6 = beta_6/(factorial(6)*(T0*1e-3)^4*abs(beta_2));
delta7 = beta_7/(factorial(7)*(T0*1e-3)^5*abs(beta_2));



% Calculate N^2 from params given in paper 
% N^2 = T0^2 * P0(10kW) * gamma / beta2
N_paper_calc = sqrt( (28.4e-15)^2 * 1e4 * gamma_SI / beta2_SI ); 

peak_power_SI = (N*A_soliton)^2; 
% peak_power_SI = (power_adjustment*A_soliton)^2; 

L_dispersion = (T0*1e-15)^2 / beta2_SI;  % 
L_NL = 1/(gamma_SI*peak_power_SI); 
disp("A_soliton = " + num2str(A_soliton) + " (sqrt Watt)")
disp("Peak power SI = " + num2str(peak_power_SI) + " (Watt)")
disp("Dispersion Length = " + num2str(L_dispersion) + " (m)")
disp("Nonlinear  Length = " + num2str(L_NL) + " (m)")

% 10000/640 ok 
% 
%---set silumation parameters
nt       = 10000; %2^12;                % N in Agrwal. 
Tmax     = 200; %100;                 % Tm in Agrwal; FFT points and window size
step_num = 8000; % 1000; % 500;                 % round(40*distance*N^2); 
zstep    = 100;                 % # of z steps
ng       = step_num/zstep;      % store data every 'ng' steps to 3D 
deltaz   = distance/step_num;   % step size in z
dtau     = (2*Tmax)/nt;         % \delta t in Agrwal; step size in tau 
domega   = (pi/Tmax);           % \delta \omega in Agrwal 

%-------tau and omega arrays
% When nt = 2048, (-nt/2:nt/2-1 ) is -1024 to +1023
tau = (-nt/2:nt/2-1)*dtau;                 % temporal grid
omega = fftshift((-nt/2:nt/2-1)*domega );  % freq array

%-----Input Field profile
if mshape == 0
    uu_ini = N * sech(tau).*exp(-1i* (chirp0*tau.^2 + chirp_3rd_order*tau.^3 ) );     %soliton;
    % uu_ini = power_adjustment * sech(tau) .^ (1 + 1i) ;    
else
    uu_ini = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape)); %super-Gaussian
end
uu = uu_ini; 


%-----Defination of Raman response function h(t) 
fR = 0.18; %0.18; %0.245;  % fraction Raman
fb = 0.21;   % fraction Boson peak 
tau1 =  12.2/T0; 
tau2 =  32/T0; 
tau_b = 96/T0;
h = (1-fb)*(tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-tau/tau2).*sin(tau/tau1)...
      + fb*((2*tau_b-tau)/tau_b^2).*exp(-tau/tau_b);
h(1:nt/2) = 0;      % causality


%-----Input 3D data   
uu_sq_3d(1,:) = uu.*conj(uu);      % Intensity in time domain
temp = fftshift(ifft(uu));    % use ifft to do FT 
spect_sq_3d(1,:)=abs(temp).^2;   % spectrum au unit



% ** N * sum(abs(ifft(A))^2) = sum(abs(A)^2)*dt/(h*f_0)
% ** Photon number at each frequency bin = N*abs(ifft(A))^2
% ***** Assume total number of photons is correct when calculating from
% time domain 


%-----Calculate pulse energy 
center_freq = 2.998e8 / (lam0 * 1e-9); % in Hz 
% pulse_total_energy = 1e9* (T0 * 1e-15 * dtau )  * trapz( (abs(uu_ini)*A_soliton * N ).^2 ); 
pulse_total_energy = (T0*1e-15*dtau) * trapz( abs(uu_ini).^2 ) * A_soliton^2 ; % in J
disp("Total pulse energy = " + num2str(1e9 * pulse_total_energy) + " nJ")



% ----- Parseval's theorem in matlab ifft/fft ; sum and trapz yields same result 
uu_sq_init_integral = trapz( abs(uu).^2 )
uu_ifft = fftshift(ifft( uu )); 
spect_sq_init_integral = trapz(abs(uu_ifft).^2)*nt

% ----- check conservation of photon number before and after FT
dt_SI =  dtau * T0*1e-15; 
domega_SI = domega / (T0*1e-15); 
total_photon_number_init_temp   = trapz( abs(uu).^2 ) * (dt_SI)              * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_init_spec   = trapz(abs(uu_ifft).^2)*nt * (dt_SI)        * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_init_spec_2 = trapz(abs(uu_ifft).^2)* 2*pi / (domega_SI) * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_init_spec_3 = trapz(abs(uu_ifft).^2)* nt^2 * dt_SI^2 * domega_SI / (2*pi) * A_soliton^2 / (6.626e-34 * center_freq)
% ---> 3 equivalent ways to calculate total number of photons in pulse 

% total_initial_photon_spectral = trapz(init_spectrum_sq) 
% trapz((init_spectrum_sq * A_soliton.^2 * dtau*T0*1e-15)/(6.626e-34 * center_freq) ) * domega/(2*pi*T0*1e-15)



%----store dispersive phase shifts to speedup code
dispersion = exp(1i*deltaz* ...
    (0.5*sbeta2*omega.^2 +...
    delta3*omega.^3 + ...
    delta4*omega.^4 + ...
    delta5*omega.^5 + ...
    delta6*omega.^6 + ...
    delta7*omega.^7 )); 

% hhz = 1i*N^2*deltaz;    % similar to 'h' in Agrwal; original 
hhz = 1i*deltaz;    % Modified 07/02/2025: N^2 absorbed into uu 

  
%*********[ Beginning of MAIN Loop]***********
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear

temp = uu.*exp(abs(uu).^2.*hhz/2);    % note hhz/2; This is N_hat * h/2 

for n = 1:step_num
    f_temp = ifft(temp).*dispersion;    % Transform to Fourier to apply dispersion; D(omega) 
    uu   = fft(f_temp);                 % Inverse Transform to time domain 

    sst1 = deriv1(abs(uu).^2,dtau); % off
    sst2 = deriv1(uu,dtau) ; % off
    sst  = 1i*(1-fR)*s*(sst1 + conj(uu).*sst2) ;   % off
    
    % Eq. 2.3.36,  Calculate the convolution of 'R' and '|A|^2'
    P      = uu.*conj(uu);
    convl  = (nt*dtau)*fft(ifft(h).*ifft(P));    % Time domain convolution <--> Freq domain multplication
    convl  = fftshift(convl);
    sst3   = deriv1(convl,dtau) ; % off
    sstnew = 1i*s*fR*(sst3 + (conj(uu).*convl)./(P+eps).*sst2);  %???

    % This part does N_hat * h; 
    % N_hat slightly different (raman and self-steepening added); 

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

    spect  = fftshift(ifft(temp));

    % take sample every ng steps 
    if mod(n,ng) == 0
         uu_sq_3d((n/ng + 1),:)    =  temp.*conj(temp);  % temporal power
         spect_sq_3d((n/ng + 1),:) = spect.*conj(spect); % spectral power
    end
end
%***************[ End of MAIN Loop ]**************

temp = temp.*exp(-1*abs(temp).^2.*hhz/2); % note hhz/2; This is -1*N_hat*h/2 
spect  = fftshift(ifft(temp)); 


% n_2 = 2.6e-20; % m^2/W
% A_eff = 2 * (1e-6)^2; 
% gamma_cal = 2*pi*n_2/(lam0*1e-9 * A_eff )   ~ 0.1
 
%----Plot output pulse shape and spectrum 
figure() 
% colormap(1-jet)
colormap(jet)
z = linspace(0,distance,zstep+1); % normalised distance





subplot(1,2,1);
% uu_sq_3d = 10*log10(uu_sq_3d);      % convert to dB units
% pmax = max(max(uu_sq_3d));       % maximum value
% pcolor(tau,z,uu_sq_3d);          % 3D plot
pcolor(tau,z, 10*log10(uu_sq_3d) );          % 3D plot
pmax = max(max( 10*log10(uu_sq_3d) ));       % maximum value
% clim([pmax-100, pmax]); 
clim([pmax-40, pmax]); 
xlim([-50,50]); shading interp;
set(gca,'FontSize', 12); %title('(a)')
xlabel('Time (t/T_0)'); ylabel('Distance (z/L_D)')
colorbar

subplot(1,2,2)
% spect_sq_3d=10*log10(spect_sq_3d);    % convert to dB units
% spect_sq_3d=spect_sq_3d;
freq = fftshift(omega/(2*pi));  % freq array; --> linear freq ? 
% pmax = max(max(spect_sq_3d));      % max value for scaling plot
pmax = max(max( 10*log10(spect_sq_3d) ));       % maximum value
% pcolor(freq,z,spect_sq_3d);
pcolor(freq,z, 10*log10(spect_sq_3d) );
hold on 
% clim([pmax-500, pmax]); 
clim([pmax-40, pmax]); 
xlim([-4,4]); shading interp;
set(gca,'FontSize', 12); %title('(a)')
xlabel('(\nu-\nu_0)T_0'); ylabel('Distance (z/L_D)')
colorbar

% At t=0 calculate average freq; in normalized units 
ini_index = 2; 
mean_freq_ini = trapz(freq, freq .* spect_sq_3d(ini_index ,:) ) / trapz( spect_sq_3d(ini_index  ,:) );

% At the end, calculate average freq  
mean_freq_end = trapz(freq, freq .* spect_sq_3d(end-1,:) )/  trapz(freq, spect_sq_3d(end-1,:) );
% Then plot this point on spectrum ... 
plot(mean_freq_end, z(end-1), "kx", markersize=15 )


% Original figure: Amplitude^2 in log unit
% subplot(1,2,1);
% uu_sq_3d = 10*log10(uu_sq_3d);      % convert to dB units
% pmax = max(max(uu_sq_3d));       % maximum value
% pcolor(tau,z,uu_sq_3d);          % 3D plot
% clim([pmax-50, pmax]); 
% xlim([-50,50]); shading interp;
% set(gca,'FontSize', 12); %title('(a)')
% xlabel('Time (t/T_0)'); ylabel('Distance (z/L_D)')
% subplot(1,2,2)
% spect_sq_3d=10*log10(spect_sq_3d);    % convert to dB units
% % spect_sq_3d=spect_sq_3d;
% freq = fftshift(omega/(2*pi));  % freq array
% pmax = max(max(spect_sq_3d));      % max value for scaling plot
% pcolor(freq,z,spect_sq_3d);
% hold on 
% clim([pmax-50, pmax]); 
% xlim([-4,4]); shading interp;
% set(gca,'FontSize', 12); %title('(a)')
% xlabel('(\nu-\nu_0)T_0'); ylabel('Distance (z/L_D)')


% ------- Initial/Final temporal signal as function of time 
% figure() 
% plot(tau*T0 , uu_sq_3d(2,:), tau, uu_sq_3d(100,:) );
% plot(tau*T0 , uu_sq_3d(1,:), tau*T0, uu_sq_3d(end,:));
% xlabel('time (fs)'); ylabel(' uu unit ');
% xlim([-10,10]);
% ------- Initial/Final spectrum as function of (diff) frequency 
% figure() 
% plot(freq , 10*log10(init_spectrum),freq, 10*log10(final_spectrum) );
% xlabel('freq. (a.u.)'); ylabel(' dB unit ');
% pmax = max(max( 10*log10(spect_sq_3d ) ));
% ylim([pmax-30,pmax]);


init_spectrum_sq = spect_sq_3d(1,:); 
final_spectrum_sq = spect_sq_3d(end,:); 



diff_freq_hz = freq ./ (T0 * 1e-15); 
linear_freq_hz = center_freq + diff_freq_hz; % A fraction of the freq < 0 

% convert to wavelength (nm) 
wavelength_nm = 1e9 * 2.998e8 ./ linear_freq_hz; % only keep positive freq
[wavelength_nm_sorted, sorted_idx] = sort(wavelength_nm);
zero_crossing_idx = find(wavelength_nm_sorted > 0, 1 ); 

% y_sorted = y(idx);
init_spectrum_wvl_sorted = init_spectrum_sq(sorted_idx); 
final_spectrum_wvl_sorted = final_spectrum_sq(sorted_idx); 
wavelength_nm_sorted = wavelength_nm_sorted(zero_crossing_idx:end); 
init_spectrum_wvl_sorted = init_spectrum_wvl_sorted(zero_crossing_idx:end); 
final_spectrum_wvl_sorted = final_spectrum_wvl_sorted(zero_crossing_idx:end); 

% Initial/Final spectrum (in uu unit) as function of wvl 
% figure() 
% plot( wavelength_nm_sorted , 10*log10(init_spectrum_wvl_sorted), wavelength_nm_sorted, 10*log10(final_spectrum_wvl_sorted) );
% xlabel('wvl (nm)'); ylabel(' ifft(uu)^2 in dB ');
% xlim([400,1400]);
% pmax = max(max( 10*log10(init_spectrum_wvl_sorted) ));
% ylim([pmax-60,pmax]);


% Initial/Final spectrum (in photon unit) as function of wvl 
eta = A_soliton * 1; 
init_spectrum_wvl_sorted_photon_unit = ((init_spectrum_wvl_sorted * eta).^2 * dtau*T0*1e-15)/(6.626e-34 * center_freq); 
final_spectrum_wvl_sorted_photon_unit = ((final_spectrum_wvl_sorted * eta).^2 * dtau*T0*1e-15)/(6.626e-34 * center_freq); 
figure() 
plot( wavelength_nm_sorted , log10(init_spectrum_wvl_sorted_photon_unit), wavelength_nm_sorted, log10(final_spectrum_wvl_sorted_photon_unit) );
xlabel('wvl (nm)'); ylabel('Photon unit: log10( n(\lambda)) ');
xlim([400,1400]);
pmax = max(max( log10(init_spectrum_wvl_sorted_photon_unit) ));
% ylim([pmax-15,pmax]);
ylim([-15, 0]);



figure() 
plot( tau, 10*log10(uu_sq_3d(end, :)))


% ------- calculate total number of photoncs: initial and final 

total_photon_number_final_temp   = trapz( uu_sq_3d(end, :) ) * (dt_SI)              * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_final_spec   = trapz( spect_sq_3d(end, :) )* nt * (dt_SI)        * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_final_spec_2 = trapz( spect_sq_3d(end, :) )* 2*pi / (domega_SI) * A_soliton^2 / (6.626e-34 * center_freq)
total_photon_number_final_spec_3 = trapz( spect_sq_3d(end, :) )* nt^2 * dt_SI^2 * domega_SI / (2*pi) * A_soliton^2 / (6.626e-34 * center_freq)



% ------ Plot Spectogram
% temp = fftshift(omega)/(2*pi);  %store frequecnies
% sample=4; gate_w=1;
% spxf = xfrog(gate_w,sample,tau,uu);
% figure(9);
% pcolor(tau(1:sample:nt),temp(1:sample:nt), 10*log10(spxf+1e-5)');
% shading interp;

toc
  
  