function samples = invert(Nchan,OS_Nu,OS_De,Nin,fname_in)
% Combines mutiple sub-channel pass-band chunks (from an oversampled PFB
% that have had their transition bands discarded) into a single contiguous
% block, then inverse FFTs.
% Ian Morrison
% 26-8-16
 
% load and concatenate the chunks
chan = 1;
load(strcat(fname_in,int2str(chan),'.mat'));
FFFF = F1((length(F1)/2)+1:length(F1)); % upper half is first part of FFFF
for chan = 2 : Nchan,
    load(strcat(fname_in,int2str(chan),'.mat'));
    FFFF = [FFFF; F1];
end;
chan = 1;
load(strcat(fname_in,int2str(chan),'.mat'));
FFFF = [FFFF; F1(1:(length(F1)/2))]; % lower half is last part of FFFF

len = length(FFFF);

save('N_channels','FFFF');

% figure;
% subplot(211); plot((1:len),abs(FFFF)); box on; grid on; title('FFFF Mag'); 
% subplot(212); plot((1:len),angle(FFFF)); box on; grid on; title('FFFF Phase'); xlabel('time');

% back transform
z1 = (ifft((FFFF), len))./(OS_Nu/OS_De);  % re-scale by OS factor

% figure;
% subplot(211); plot((1:len),real(z1(1:len))); box on; grid on; title('z1 Real'); 
% subplot(212); plot((1:len),imag(z1(1:len))); box on; grid on; title('z1 Imag'); xlabel('time');
% 
% figure;
% subplot(211); plot((1:len),10.0*log10(abs(real(z1(1:len)))+1e-12)); box on; grid on; title('z1 Real - Log scale');
% axis([1 len -100 10]);
% subplot(212); plot((1:len),10.0*log10(abs(imag(z1(1:len)))+1e-12)); box on; grid on; title('z1 Imag - Log scale'); xlabel('time');
% axis([1 len -100 10]);

samples = z1;

return
end

