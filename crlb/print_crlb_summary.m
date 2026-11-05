function print_crlb_summary(R)
% Print the default-SNR summary.

fprintf('\nRESULTS AT THE DEFAULT SNR\n');
fprintf('UE PEB   = %.6g m\n', R.uePEB_m);
fprintf('UE VEB   = %.6g m/s\n', R.ueVEB_mps);
fprintf('eta LB   = %.6g ppm\n', R.etaStdPpm);
fprintf('delta LB = %.6g ns\n', R.deltaStdNs);
fprintf('FIM rank = %d/%d, condition number = %.6g\n', ...
    R.infoEquivalent.rank, R.infoEquivalent.dimension, ...
    R.infoEquivalent.conditionNumber);
end
