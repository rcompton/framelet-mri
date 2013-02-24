%ugh
load plt_excibig.mat
errors_both_perbreg = errors_per_breg;
errors_both_periter = errors;
errorsr_both_perbreg = res_error_per_breg;
errorsr_both_periter = res_errors;
biters = iters;

load plt_nu_only.mat
errors_tvonly_perbreg = errors_per_breg;
errors_tvonly_periter = errors;
errorsr_tvonly_perbreg = res_error_per_breg;
errorsr_tvonly_periter = res_errors;
tviters = iters;

load plt_exci_only.mat
errors_fonly_perbreg = errors_per_breg;
errors_fonly_periter = errors;
errorsr_fonly_perbreg = res_error_per_breg;
errorsr_fonly_periter = res_errors;
foiters = iters;

load plt_excism.mat
errors_bothsm_perbreg = errors_per_breg;
errors_bothsm_periter = errors;
errorsr_bothsm_perbreg = res_error_per_breg;
errorsr_bothsm_periter = res_errors;
biterssm = iters;

%%

figure()
hold on;
semilogy(errors_both_perbreg,'bs-');
semilogy(errors_fonly_perbreg,'g');
semilogy(errors_tvonly_perbreg,'r');
semilogy(errors_bothsm_perbreg,'bx-')
xlabel('Bregman iterations');
ylabel('error: log(||\rho - \rho_{exact}||)');
%%
figure()
semilogy(errorsr_both_perbreg,'b');
hold on;
semilogy(errorsr_fonly_perbreg,'g');
semilogy(errorsr_tvonly_perbreg,'r');
xlabel('Bregman iterations');
ylabel('residual error: log(||A\rho - s||)');
%%

figure()
semilogy(cumsum(biters), errors_both_periter,'bs-');
hold on;
semilogy(cumsum(foiters), errors_fonly_periter,'g');
semilogy(cumsum(tviters), errors_tvonly_periter,'r');
xlabel('matrix multiplications');
ylabel('error: log(||\rho - \rho_{exact}||)');
%%

figure()
semilogy(cumsum(biters), errorsr_both_periter(2:end),'bs-');
hold on;
semilogy(cumsum(foiters), errorsr_fonly_periter(2:end),'g');
semilogy(cumsum(tviters), errorsr_tvonly_periter(2:end),'r');
xlabel('matrix multiplications');
ylabel('residual error: log(||A\rho - s||)');



%% Brain Image

load brain_hybrid_128_again.mat
errors_breg_hyb128 = errors_per_breg;
errorsr_breg_hyb128 = errorsr_per_breg;
errors_iter_hyb128 = errors;
errorsr_iter_hyb128 = errorsr;
hybiters = iters;

load brain_fonly_128_again.mat
errors_breg_fonly128 = errors_per_breg;
errorsr_breg_fonly128 = errorsr_per_breg;
errors_iter_fonly128 = errors;
errorsr_iter_fonly128 = errorsr;
fiters = iters;

load brain_tvonly_128_again.mat
errors_breg_tvonly128 = errors_per_breg;
errorsr_breg_tvonly128 = errorsr_per_breg;
errors_iter_tvonly128 = errors;
errorsr_iter_tvonly128 = errorsr;
tviters = iters;

%%

figure()
semilogy(errors_breg_hyb128,'b');
hold on;
semilogy(errors_breg_fonly128,'g');
semilogy(errors_breg_tvonly128,'r');
xlabel('Bregman iterations');
ylabel('error: log(||\rho - \rho_{exact}||)');

figure()
semilogy(errorsr_breg_hyb128,'b');
hold on;
semilogy(errorsr_breg_fonly128,'g');
semilogy(errorsr_breg_tvonly128,'r');
xlabel('Bregman iterations');
ylabel('residual error: log(||A\rho - s||)');


figure()
semilogy(cumsum(hybiters), errors_iter_hyb128,'b');
hold on;
semilogy(cumsum(fiters), errors_iter_fonly128,'g');
semilogy(cumsum(tviters), errors_iter_tvonly128,'r');

xlabel('matrix multiplications');
ylabel('error: log(||\rho - \rho_{exact}||)');


figure()
semilogy(cumsum(hybiters), errorsr_iter_hyb128,'b');
hold on;
semilogy(cumsum(fiters), errorsr_iter_fonly128,'g');
semilogy(cumsum(tviters), errorsr_iter_tvonly128,'r');
xlabel('matrix multiplications');
ylabel('residual error: log(||A\rho - s||)');


