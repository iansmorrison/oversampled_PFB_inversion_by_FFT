function compare(samples,start,len,fname_compare,compare_offset)
% Compares waveform in "samples" with the original - in both the time and
% frequency domains.
% Ian Morrison
% 26-8-16

fid_in = fopen(fname_compare);

% Read original input to PFB
Vstream = single(fread(fid_in, len, 'single'));

if feof(fid_in)
    error('Error - hit end of input file!');
end;

Vdat = reshape(Vstream, 1, []);

% Look at original waveform
figure;
subplot(211); plot((1:len),real(Vdat(1,1:len))); box on; grid on;
title('Input Real'); 
subplot(212); plot((1:len),imag(Vdat(1,1:len))); box on; grid on;
title('Input Imag'); xlabel('time');

% Look at reconstituted waveform
figure;
subplot(211); plot((1:len),real(samples(1:len))); box on; grid on; title('Reconstituted Output Real'); 
subplot(212); plot((1:len),imag(samples(1:len))); box on; grid on; title('Reconstituted Output Imag'); xlabel('time');

figure;
subplot(211); plot((1:len),10.0*log10(abs(real(samples(1:len)))+1e-12)); box on; grid on; title('Reconstituted Output Real - Log scale');
axis([1 len -100 10]);
subplot(212); plot((1:len),10.0*log10(abs(imag(samples(1:len)))+1e-12)); box on; grid on; title('Reconstituted Output Imag - Log scale'); xlabel('time');
axis([1 len -100 10]);


% Time domain comparison with original input - good for integer sample
% delays when those delays have been sync'd out
if (1)
    centre_Vdat = start;
    centre_samples = centre_Vdat + compare_offset;
 
    plot_range = 25;
    figure;
    subplot(211); plot((-plot_range+1:plot_range),real(samples(centre_samples-plot_range+1:centre_samples+plot_range)), (-plot_range+1:plot_range),real(Vdat(1,centre_Vdat-plot_range+1:centre_Vdat+plot_range))); box on; grid on; title('Reconstituted Output vs Input Real'); 
    subplot(212); plot((-plot_range+1:plot_range),imag(samples(centre_samples-plot_range+1:centre_samples+plot_range)), (-plot_range+1:plot_range),imag(Vdat(1,centre_Vdat-plot_range+1:centre_Vdat+plot_range))); box on; grid on; title('Reconstituted Output vs Input Imag'); xlabel('time');

    % calculate total RMS error over Nsamp samples
    Nsamp = 200;
    Realerr = 0;
    Imagerr = 0;

    for j = centre_Vdat-(Nsamp/2):centre_Vdat+(Nsamp/2-1),
        Realerr = Realerr + (real(Vdat(1,j)) - real(samples(j+compare_offset)))^2;
        Imagerr = Imagerr + (imag(Vdat(1,j)) - imag(samples(j+compare_offset)))^2;
    end;

    RMSerr_real = (Realerr/Nsamp)^0.5
    RMSerr_imag = (Imagerr/Nsamp)^0.5
            
end;


% Frequency domain comparison with original input - good for fractional
% sample delays to avoid the need for interpolation of the input time series
if (1)
    Vdat_shift = Vdat(1,1-compare_offset:len);
    VDAT = fft(Vdat_shift);
%     figure;
%     subplot(211); plot((1:len+compare_offset),abs(VDAT)); box on; grid on; title('Input Mag'); 
%     subplot(212); plot((1:len+compare_offset),angle(VDAT)); box on; grid on; title('Input Phase'); xlabel('time');
    
    samples_shift = samples(1:len+compare_offset);
    SAMPLES = fft(samples_shift);
%     figure;
%     subplot(211); plot((1:len+compare_offset),abs(SAMPLES)); box on; grid on; title('Output Mag'); 
%     subplot(212); plot((1:len+compare_offset),angle(SAMPLES)); box on; grid on; title('Output Phase'); xlabel('time');
    
    % cross-power spectrum of output and a similar length input
    CP = SAMPLES.*conj(VDAT);
    
    figure;
    subplot(211); plot((1:len+compare_offset),abs(CP)); box on; grid on; title('Cross-Power Mag'); 
    subplot(212); plot((1:len+compare_offset),angle(CP)); box on; grid on; title('Cross-Power Phase'); xlabel('time');
    
    % average cross-power magnitude over whole CP spectrum
    average_cp_mag = sum(abs(CP))/length(CP)
   
  end;

return
end

