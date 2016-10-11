function design_PFB(Nchan,Num,Den,Ntaps,ffft_len,display)

    % Oversampling Factor
    OS = Num/Den;

    % Filter specs for the prototype filter
    % Cut-off frequency 
    Fp = 1./Nchan;
    % Stop-band frequency
    Fs = 1.*(2*OS-1)/Nchan;

    % Filter Transfer Function specs
    Ap = 0.01;
    As = 60;

    % Design filter
    Hf = fdesign.lowpass('N,Fp,Fst',Ntaps,Fp,Fs); 
    H_Obj_0 = design(Hf,'firls','Wstop',15,'systemobject',true);
    h = H_Obj_0.Numerator;

    % Save impulse response h, and other parameters
    save Prototype_FIR.mat h Nchan Fp Fs Ap As;
    
    % Save a sampled version of the Transfer Function for later equalisation
    % - length should be Nchan times the half-channel width (where width is FFTlength/OS_factor)
    % e.g. 64 channels, ffft_len = 1024: 28,672 is 448*64, which gives 448 points per half-channel, 896 per channel
    [H0,W] = freqz (h, 1, ffft_len*Den*Nchan/(2*Num));
    save TF_points.mat H0 W;

    % Optionally display design
    if(display==1)
        [H0,W] = freqz (h, 1, Ntaps*Nchan);

        %Rescaling the frequency axis
        W = W/pi;

        figure; 
        subplot(3,1,1)
        plot (W, abs(H0));
        axis ([0 3.5*Fp -0.15 1.15]); 
        title('Transfer Function of the Prototype Filter')
        grid on; box on; 

        subplot(3,1,2)
        hold on;
        plot (W, 20*log10(abs(H0)));
        plot([0 Fp], [-0.5*Ap -0.5*Ap],'k-.','LineWidth',1);
        plot([0 Fp], [ 0.5*Ap  0.5*Ap],'k-.','LineWidth',1);
        hold off;
        axis ([0 1.5*Fp -4.5*Ap 4.5*Ap]); 
        title ('Passband')
        grid on; box on;

        subplot (3,1,3);
        hold on;
        plot (W, 20*log10(abs(H0)));
        plot([Fs 1], [-As -As],'r-','LineWidth',1);
        hold off;
        axis ([0 1 -(As+10) 3]); 
        title ('Stopband')
        grid on; box on;
        
        pause;
    end;
    
    close all;
    
return
end
