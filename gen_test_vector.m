function gen_test_vector(Wave_type,impulse_offset,impulse_width,Nout,Nblocks,f_sample_out,period,noise,fname)
% Generates a file containing one continuous series of real samples
% (32 bit floating point I and Q), representing the desired test waveform.
% The output filename is passed as an input parameter.
%
% The output consists of Nblocks blocks, each of length Nout.
% 
% Waveform options are:
%    1. an M-sample impulse starting at a specified location within each
%       output block
%    2. one polarisation of a simulated pulsar signal with specified period
%       (and SNR)
%
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% I. Morrison      27-Jun-2016  Original version
% 
% ----------------------------------------------------------------------
 
% Local parameters
nbins = 2^10; % Number of bins within a pulse period
shift = 0; % performs an fftshift before the inverse FFT (0 if using FPB)
Tin = 2.0/abs(f_sample_out)*1E-6; % Time spacing between input data elements
Nin = Nout/2; % Number of data elements in input time series
Pmul = 0.5; % Power multiplication factor for all but the DC channel

% Vector of relative times
trel = (0:Nin-1)*Tin;
 
% Open file for writing
fid = fopen(fname, 'w');

%===============
 
for ii = 1:Nblocks,
    % Print loop number
    fprintf('Loop # %i of %i\n', ii, Nblocks);
    
    % Time vector
    if ii == 1,
        tt = trel;
    else
        tt = ttclip(end) + Tin + trel;
    end;
    
    f = mod(tt, period)/period; 
    tindex = transpose(round(f*nbins+0.5));
    index = unique(tindex);
    
  
    if(Wave_type == 1)  % Impulse case
        
        % Initialize data vector
        z = zeros(Nout,1,'single');
        
        % Introduce an impulse in z
        z = zeros(Nout,1);
        z(impulse_offset:impulse_offset+impulse_width-1,1) = 1;
        
    else  % Pulsar case
    
        % Calculate phase-dependent Stokes parameters and coherency matrix
        % using the rotating vector model
        [~, J] = rotvecmod(nbins,noise);
        
        % Initialize data vector
        z = zeros(Nin,2,'single');
        
        % Loop through groups of data that share the same phase. Random data 
        % in each group are generated from the same coherency matrix

        for jj = 1:length(index),
            %Get coherency matrix for this pulsar phase
            Jcoh = [J(index(jj),1), J(index(jj),3); ...
                    J(index(jj),2), J(index(jj),4)];

            % Indices of elements with a given phase
            iphase = find(tindex == index(jj));
            nL = length(iphase);

            %Generate two randomly-phased, unit-length phasors  
            %z0 = exp(complex(0,1)*2*pi()*rand(nL,npol));
            z0 = sqrt(0.5)*[complex(randn(nL,1),randn(nL,1)), ...
                            complex(randn(nL,1),randn(nL,1))];

            %Generate covariant vectors via Cholesky decomposition
            zjj = z0*chol(Jcoh, 'upper');

            % Concatenate with data from other phases
            z(iphase, :) = zjj;
        end;
        
        f1a = fft(z(:,1), Nin);

        f1 = [real(f1a(1)); f1a(2:Nin)*Pmul; ...
              imag(f1a(1)); flipud(conj(f1a(2:Nin)))*Pmul];
    
        % Inverse FFT
        % Optionally include an fftshift before the inverse FFT, as needed
        if shift == 1,
            f1 = fftshift(f1);
        end;
        z = ifft(f1, Nout);
    end;
 
    ttclip = tt(1:Nin);
    
    dat = reshape(transpose(z),Nout,1);

    %Write vector to file
    fwrite(fid, dat, 'single');
end;
 
fclose(fid);

return;
 
end




function [S, J, p] = rotvecmod(N, noise, showplot)
% Rotating vector model for pulsar emission
 
if ~exist('N','var'),
    N = 1024;
end;
 
esig = 5. ; % emission half angle (polar angle, degrees)
epeak = 0. ; % emission peak angle (polar angle, degrees)
flin = 0.3; % linear polarized fraction amplitude
 
zeta = 30.; % observing angle (degrees) relative to rotation axis
alpha = 40.; % magnetic axis (degrees) relative to rotation axis
 
pmin = -180.;
pmax = 180.;
 
% Angle of rotation: p=0 for aligned dipole. 
% This is equivalent to pulsar longitude or phase
p = transpose(linspace(pmin, pmax, N)); 
 
% Polarization angle w.r.t projected rotation axis from observing direction
%psi = atand(sind(alpha)*sind(p)./(sind(zeta)*cosd(alpha) - ...
%    sind(alpha)*cosd(zeta)*cosd(p)));
psi = atan2d(sind(alpha)*sind(p),  ...
    (sind(zeta)*cosd(alpha) - sind(alpha)*cosd(zeta)*cosd(p)));
 
% Polar observation angle in magnetic axis reference frame
cosO = cosd(p)*sind(zeta)*sind(alpha) + cosd(alpha)*cosd(zeta);
tanO = sqrt(1./(cosO.^2)-1);
 
% Polar emission angle in magnetic axis reference frame
thetaE = atand(1.5*(sqrt(1+(8/9)*tanO.^2) - 1)./tanO);
%thetaE = atand(1.5*(-sqrt(1+(8/9)*tanO.^2) - 1)./tanO);
 
% Intensity (model-based assumption)
S0 = (1./sqrt(2*pi()*esig^2))*exp(-(thetaE-epeak).^2/(2.*esig^2));
S0 = S0/max(S0); %normalize max to 1
 
% Linear polarization fraction (model-based assumption)
L = flin*S0.*cosd(thetaE);
 
% Other Stokes parameters
S1 = L.*cosd(2*psi);
S2 = L.*sind(2*psi);
S3 = -(1-flin)*S1; % Fake circular polarization to avoid zero signal
%S3 = single(zeros(N,1)); % Zero circular polarization component
 
% Add noise, typically such that max(S/N) = 1
S0 = S0 + noise;
 
% Normalize Stokes 4-vector so that S0 = 1. 
factor = max(S0);
S0 = S0/factor;
S1 = S1/factor;
S2 = S2/factor;
S3 = S3/factor;
 
% Create Coherency matrix
Jxx = 0.5*(S0 + S1);
Jyy = 0.5*(S0 - S1);
Jxy = 0.5*(S2 + 1i*S3);
Jyx = 0.5*(S2 - 1i*S3);
 
% Plot results, if requested. Useful for debugging.
if exist('showplot','var'),
    clf();
 
    subplot(2,2,1);
    plot(p, transpose([S0, S1, S2, S3])); 
    legend('S0', 'S1', 'S2', 'S3');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,3);
    plot(S0, transpose([S1, S2]));
    hleg1 = legend('S1', 'S2');
    set(hleg1,'Location','NorthWest')
    axis([0, 2, -Inf, Inf]);
    xlabel('S0','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S1 or S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,2);
    plot(p, transpose([Jxx, Jyy, real(Jxy), imag(Jxy)]));
    legend('Jxx', 'Jyy', 'Real(Jxy)', 'Imag(Jxy)');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,4);
    plot(S1, S2, 'b');
    xlabel('S1','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
end;
 
S = [S0, S1, S2, S3];
J = [Jxx, Jyx, Jxy, Jyy];

return;

end
