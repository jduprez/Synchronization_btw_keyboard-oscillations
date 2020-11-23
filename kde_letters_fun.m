% This function performs Kernel Density Estimation and was written by Mike X Cohen

function [kde] = kde_letters_fun(vecIKI, vecfft, fwhm)
%% create kernel

% fwhm = 1; % in Hz
gwin = exp(-4*log(2)*linspace(-3,3,101).^2 / fwhm^2);

% frequencies vector (arbitrary upper limit and resolution)
hz = linspace(0,100,1000); % Hz
emptydelta = zeros(size(hz));

% show an example
% plot(hz(1:length(gwin)),gwin)

% convolution parameters
nConv = length(gwin) + length(hz) - 1;
halfk = (length(gwin)-1)/2;
gausX = fft(gwin,nConv);

%%

% initialize kernel density estimation
kde = zeros(size(vecIKI,1),length(hz));

for wordi = 1:size(vecIKI,1)
    
    % IKIs in Hz
    wordIKIs = 1000./diff(find(vecfft{wordi}));
    nIKIs = length(wordIKIs);
    
    % find closest frequency in Hz
    closehz = dsearchn(hz',wordIKIs');
    
    % per letter
    for letti=1:nIKIs
        
        % delta function
        lethz = emptydelta;
        lethz(closehz(letti)) = 1;
        
        % kernel density via convolution with Gaussian
        tkde = ifft(fft(lethz,nConv).*gausX);
        tkde = tkde(halfk:end-halfk-1);
        
        % then add to this word's KDE (and normalize by number of letters)
        kde(wordi,:) = kde(wordi,:) + tkde/nIKIs;
    end
end

%%

