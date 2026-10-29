function [hCycles, dAmin_all, dAmax_all] = compute_zuecco_stats(W, R, time, startMonth, minLen)
    % W = "wetness" = PET
    % R = "Response" = AET
    % Same cycle-detection logic as compute_zuecco , but returns the
    % per-year vectors instead of collapsing to mean/min/max, so the
    % caller can take mean, median, or whatever summary stat it wants.

    % Zuecco, G., Penna, D., Borga, M. & van Meerveld, H. J. 
    % A versatile index to characterize hysteresis between hydrological 
    % variables at the runoff event timescale. Hydrological Processes 30, 
    % 1449–1466 (2016).

    yrs  = year(time);
    mons = month(time);

    % --- Hydrological year assignment ---
    hydroYear = yrs;
    hydroYear(mons < startMonth) = hydroYear(mons < startMonth) - 1;

    uYears = unique(hydroYear);
    uYears(isnan(uYears)) = [];

    if length(uYears) < minLen
        hCycles = nan; dAmin_all = nan; dAmax_all = nan;
        return
    end

    hCycles   = nan(length(uYears),1);
    dAmin_all = nan(length(uYears),1);
    dAmax_all = nan(length(uYears),1);

    for i = 1:length(uYears)

        idx = hydroYear == uYears(i);

        if sum(idx) < 200
            continue
        end

        Wy = W(idx);
        Ry = R(idx);

        if range(Wy) < 0.05
            continue
        end

        % --- Normalise within year ---
        u = (Wy - min(Wy)) ./ (max(Wy) - min(Wy));
        v = (Ry - min(Ry)) ./ (max(Ry) - min(Ry));

        % --- Smooth ---
        u = smoothdata(u,'movmean',30);
        v = smoothdata(v,'movmean',30);

        % --- Turning point = last seasonal max ---
        ip = find(u == max(u), 1, 'last');

        if ip <= 10 || ip >= length(u)-10
            continue
        end

        % --- Split limbs ---
        u_r = u(1:ip);   v_r = v(1:ip);   % wetting
        u_f = u(ip:end); v_f = v(ip:end); % drying

        if length(u_r) < 20 || length(u_f) < 20
            continue
        end

        % --- Sort limbs by wetness ---
        [u_r, ir] = sort(u_r);
        v_r = v_r(ir);

        [u_f, if_] = sort(u_f);
        v_f = v_f(if_);

        % --- Make wetness unique by averaging duplicates ---
        [u_r_unique, ~, ir2] = unique(u_r);
        v_r_unique = accumarray(ir2, v_r, [], @mean);

        [u_f_unique, ~, if2] = unique(u_f);
        v_f_unique = accumarray(if2, v_f, [], @mean);

        u_r = u_r_unique;
        v_r = v_r_unique;

        u_f = u_f_unique;
        v_f = v_f_unique;

        % --- Only use overlapping wetness range ---
        u_min = max(0.15, max(min(u_r), min(u_f)));
        u_max = min(max(u_r), max(u_f));

        if u_max <= u_min
            continue
        end

        u_common = linspace(u_min, u_max, 100);

        vr = interp1(u_r, v_r, u_common, 'linear');
        vf = interp1(u_f, v_f, u_common, 'linear');

        if any(isnan(vr)) || any(isnan(vf))
            continue
        end

        % --- Difference at same wetness ---
        d = vr - vf;

        hCycles(i)   = mean(d);
        dAmin_all(i) = min(d);
        dAmax_all(i) = max(d);
    end
end