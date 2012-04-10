%
% Parallel timing plot
%

ncores = feature('numcores');
ntrials = 20;

recon_times = zeros(ncores,ntrials);

for aye = 1:ncores
    maxNumCompThreads(enn);
    
    for jay=1:ntrials
        tic;
        attempt_CS_framelets_3d;
        recon_times(aye, jay) = toc;
    end
end

save parall_data.mat

%plot(recon_times,'x-');
