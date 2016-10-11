function fine_chan_proc(chan,Nin,OS_Nu,OS_De,input_offset,fname_in,fname_out,equalise_ripple)
%
% Reads one block of length Nin from a fine channel of the PFB, performs a
% forward FFT and discards the oversampled portions (transition bands).
% Optionally applies pass-band equalisation (de-ripple).
% Optionally applies a phase shift of value that depends on the channel number.

savefile = fname_out;

fid = fopen(fname_in);
 
% Shift starting point for reading file (8 bytes per complex sample)
fseek(fid, input_offset*8, 'cof');

% Read stream of complex voltages, forming a single column
Vstream = single(fread(fid, 2*Nin, 'single'));

if feof(fid)
    error('Error - hit end of input file!');
end;

% Parse real and imag components
Vstream = reshape(Vstream, 2, []);
Vstream = complex(Vstream(1,:), Vstream(2,:));
    
Vdat = reshape(Vstream, 1, []);


% Optional phase shift - 8 channels only at present
phase_shift = 1.0;
if (0)
    switch chan
        case 1
            phase_shift = 0.; % channel 1 must be 0: can't rotate upper and lower sidebands by different amounts with a single complex multiplier
        case 2
            phase_shift = 1i;
        case 3
            phase_shift = 0.5 + (sqrt(3.0)/2.0)*1i;
        case 4
            phase_shift = sqrt(3.0)/2.0 + 0.5i;
        case 5
            phase_shift = 1;
        case 6
            phase_shift = sqrt(3.0)/2.0 - 0.5i;
        case 7
            phase_shift = 0.5 - (sqrt(3.0)/2.0)*1i;
        case 8
            phase_shift = -1i;
    end;
end;


% Forward FFT
F1 = fftshift(fft(Vdat(1,:).*phase_shift, Nin));
F1 = reshape(F1, Nin, 1);


% Keep only the pass-band
discard = (1.0 - (OS_De/OS_Nu))/2.0;
F1 = F1(round(discard*Nin)+1:round((1.0-discard)*Nin));

 
% Optionally equalise PFB pass-band ripple
if(equalise_ripple)
    % load PFB prototype filter transfer function
    load('TF_points.mat');

    % use just the baseband passband section of transfer function
    % - apply to both halves of channel
    passband_len = (Nin/2)*OS_De/OS_Nu;
    for ii = 1:passband_len,
        F1(ii) = F1(ii)/abs(H0(passband_len-ii+2));
        F1(passband_len+ii) = F1(passband_len+ii)/abs(H0(ii+1));
    end;
end;

save(savefile,'F1');

fclose(fid);
 
return
end
 

