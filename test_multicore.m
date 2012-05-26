%
% Test out the 12 cores
%
clear all; close all;

X = peaks(10028);

ncores = feature('numcores');
ntrials = 20;

mtx_power_times = zeros(ncores,ntrials);
fft_times = zeros(ncores, ntrials);

for i=1:ncores
    
    fftw('planner','exhaustive')
    tic;fft(X,[],1);toc
    
    for j=1:ntrials
        
        maxNumCompThreads(i);

        tic;
        %X^2;
        mtx_power_times(i,j) = toc;
        
        tic
        fft(X,[],1);
        fft_times(i,j) = toc;
        
    end
end


subplot(1,2,1);
plot(mtx_power_times,'x-')
title('mtx power time vs number of cores');

subplot(1,2,2);
plot(fft_times,'x-');
title('fftn time vs num of cores');
