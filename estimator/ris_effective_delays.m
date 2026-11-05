function a = ris_effective_delays(theta, times, cfg, idx)
% Compute a_lk = tau_lk - (nu_lk-eta)T_l - delta_c.

s = ris_unpack_state(theta, idx);
L = numel(times);
a = zeros(L, idx.M);

for ell = 1:L
    T = times(ell);
    if strcmpi(cfg.geometryModel, 'time-varying')
        u = s.u0 + s.uDot*T;
        p = s.p0 + s.pDot*T;
    else
        u = s.u0;
        p = s.p0;
    end

    d0 = u - cfg.q;
    r0 = norm(d0);
    check_distance(r0, 'UE-BS');
    tau = r0/cfg.c;
    nu = (s.uDot.'*d0)/(cfg.c*r0);
    a(ell, 1) = tau - (nu-s.eta)*T - s.deltaC;

    for k = 1:idx.K
        dBR = cfg.q - p(:, k);
        dRU = p(:, k) - u;
        rBR = norm(dBR);
        rRU = norm(dRU);
        check_distance(rBR, sprintf('BS-RIS%d', k));
        check_distance(rRU, sprintf('RIS%d-UE', k));

        tau = (rBR+rRU)/cfg.c;
        nu = (-(s.pDot(:, k).')*dBR/rBR ...
            + (s.pDot(:, k)-s.uDot).'*dRU/rRU)/cfg.c;
        a(ell, k+1) = tau - (nu-s.eta)*T - s.deltaC;
    end
end
end


function check_distance(d, label)
if d < 1e-8
    error('Degenerate distance: %s.', label);
end
end
