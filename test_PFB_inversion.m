%% Model for evaluating inversion by FFT of an oversampled PFB
%  Supports:
%     - variable number of PFB channels
%     - variable PFB oversampling factor
%     - variable number of PFB taps per channel
%     - variable forward FFT length
%     - variable test data block length
%     - continuous reconstruction over multiple blocks
%       (with variable number of blocks and overlap fraction)
%     - choice of impulse or simulated pulsar test waveforms
%     - optional PFB response equalisation (de-ripple)
% 
%  Ian Morrison, Swinburne Centre for Astrophysics and Supercomputing
%  October 2016
%

fprintf('\nTest of OS-PFB Inversion via FFT\n');


%% GLOBAL PARAMETERS

% Number of PFB output channels - power of 2, min OS_Nu, max 256
N = 32;

% PFB oversampling factor
OS_Nu = 32;  % numerator - should be a sub-multiple of N
OS_De = 27;  % denominator

% Width of PFB channel passband in MHz = spacing of PFB output channels
fine_chan_passband = 0.003617;

% Length of forward FFT to process fine channels
ffft_length = 2^10;

% Length of test vector blocks (spacing of impusles)
block_length = N*ffft_length;

% Number of blocks
Nblocks = 50;

 
%% GENERATE TEST VECTOR (input to PFB)

test_vector_filename = 'test_vec.dump';

Wave_type = 0;  % 0 for pulsar, 1 for impulse
impulse_offset = block_length/4;  % location of impulse within each block
impulse_width = 1;  % number of samples width of impusle
f_sample_out = N*fine_chan_passband;  % sample rate in MHz
period = 0.01;  % simulated pulsar period in seconds
noise = 1;  % sets SNR of simulated pulsar signal

% function gen_test_vector(Wave_type,impulse_offset,impulse_width,Nout,Nblocks,noise,f_sample_out,period,fname)
fprintf('\nGenerating test vector...\n');
gen_test_vector(Wave_type,impulse_offset,impulse_width,block_length,Nblocks,f_sample_out,period,noise,test_vector_filename);


%% DESIGN PFB PROTOTYPE FILTER
% function design_PFB(Nchan,OS_Nu,OS_De,Ntaps,ffft_len,display)
taps_per_chan = 16;
Ntaps = taps_per_chan*N;

display = 1;    % 1 to display filter design plot, 0 otherwise

fprintf('\nDesigning PFB prototype filter...\n');
if (display == 1)
    fprintf('\nPress any key to continue...\n');
end;
design_PFB(N,OS_Nu,OS_De,Ntaps,ffft_length,display);
  

%% PFB Channelize - all blocks
% function PFBchannelizer(Nchan,OS_Nu,OS_De,Nin,Nblocks,fname_in,fname_out)
% minimum Nin is (block_length/OS_factor) - can be longer
fprintf('\nChannelizing...\n');
PFB_channelizer(N,OS_Nu,OS_De,OS_De*block_length/OS_Nu,Nblocks,test_vector_filename,'fine_channel_');

 
%% PROCESS EACH FINE CHANNEL
initial_input_offset = 128;  % number of samples to drop at the start of the PFB output data, to ensure impulse within window
input_offset = initial_input_offset;
equalise_ripple = 1;  % 1 to equalise PFB ripple, 0 to not
overlap_fraction = 0.125;
overlap = overlap_fraction*ffft_length;
Num_overlapped_blocks = Nblocks + floor(Nblocks*overlap_fraction);
usable_length = N*(OS_De*ffft_length/OS_Nu)*(1 - overlap_fraction);
usable_start = N*overlap*OS_De/OS_Nu/2 + 1;
total_usable = Num_overlapped_blocks*usable_length;
output_samples = complex(zeros(1,total_usable));
sample_index = 1;
% Process each block
fprintf('\nProcessing blocks...\n');
for block = 1:Num_overlapped_blocks,
    fprintf('\nBlock %d\n', block);   
    
    % process each channel of the current block
    fprintf('Processing each channel...\n');
    for chan = 1:N,
        % function fine_chan_proc(chan,Nin,OS_Nu,OS_De,input_offset,fname_in,fname_out,equalise_ripple)
        fine_chan_proc(chan,ffft_length,OS_Nu,OS_De,input_offset,strcat('fine_channel_',int2str(chan),'.dump'),strcat('chunk_',int2str(chan),'.mat'),equalise_ripple);
    end;
    
    % combine chunks and back-transform
    % function invert(Nchan,OS_Nu,OS_De,Nin,fname_in)
    fprintf('Combining channels and back transforming...\n');
    samples = invert(N,OS_Nu,OS_De,block_length,'chunk_');
    output_samples(1,sample_index:sample_index+usable_length-1) = samples(usable_start:usable_start+usable_length-1);
    sample_index = sample_index + usable_length;
    input_offset = input_offset + ffft_length - overlap;
end;


%% COMPARE RECONSTRUCTED WAVEFORM WITH ORIGINAL
% function compare(samples,start,length,fname_compare,compare_offset)
compare_offset = Ntaps/2 + 1 - (OS_De*N/OS_Nu)*initial_input_offset - (usable_start - 1);
compare(output_samples,total_usable/2,total_usable,test_vector_filename,compare_offset); % compare 200 samples around the middle of the output sequence


%% SAVE RECONSTRUCTED WAVEFORM TO FILE

reconstructed_filename = 'reconstructed.dump';
fid = fopen(reconstructed_filename, 'w');
fwrite(fid, transpose(output_samples), 'single');
fclose(fid);


%% Finished
fprintf('\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
