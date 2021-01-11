function out = mut_dist_base(rho, grp, delta, add_dist, cmp_dist)
%MUT_DIST Summary of this function goes here
%   Detailed explanation goes here

    %% Getting information.
    m = height(rho);
    n = width(rho);
    total_per_day = sum(rho, 1);
    days = total_per_day > 1e-5;

    %% Compute maximum delta amount.
    cmp_dist = cmp_dist - min(cmp_dist);
    cmp_dist = cmp_dist / sum(cmp_dist);

    available_values = cmp_dist(days) / max(cmp_dist(days)) .* rho(grp,days);

    available_total = sum(available_values);
    delta = min([delta, available_total]);
    if delta < 1e-7
        out = rho;
        return
    end

    available_dist = available_values / available_total;

    %% Compute compensation for addition.
    add_dist = add_dist - min(add_dist);
    add_dist = add_dist / sum(add_dist);
    compensate_add = zeros(m, n);
    compensate_add(:,days) = ((rho(:,days) ./ total_per_day(:,days)) .* add_dist(:,days) * delta);

    %% Compute compensation for substraction.
    compensate_rem = zeros(m, n);
    compensate_rem(:,days) = sum(compensate_add(:,days), 2) * available_dist;

    %% Create the new matrix.
    out = rho - compensate_add + compensate_rem;
    out(grp,:) = out(grp,:) + sum(compensate_add, 1) - sum(compensate_rem, 1);
end

