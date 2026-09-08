%% Step 5: compute the static RMPC robustness constraint X_R^i minus D^i
% For X_R = {x | H*x <= g} and D = E*[uPrevMin,uPrevMax],
% X_R minus D = {x | H*x <= g - max_{d in D}(H*d)}.

clear; clc; close all;
mpt_init;

if ~isfile('robust_sets_step4.mat')
    error(['robust_sets_step4.mat was not found. ', ...
        'Run step4_robust_invariant_set.m first.']);
end
load('robust_sets_step4.mat', ...
    'X', 'Xrobust', 'E', 'uPrevMin', 'uPrevMax');

Xrobust.minHRep();
H = Xrobust.A;
g = Xrobust.b;

% Worst-case displacement of each inequality face caused by d in D.
HdAtLower = H*E*uPrevMin;
HdAtUpper = H*E*uPrevMax;
robustMargin = max(HdAtLower, HdAtUpper);

XRminusD = Polyhedron('A', H, 'b', g - robustMargin);
XRminusD.minHRep();

if XRminusD.isEmptySet()
    error('X_R^i minus D^i is empty: this RMPC configuration is infeasible.');
end

origin = [0; 0];
originContained = XRminusD.contains(origin);
slackAtOrigin = XRminusD.b - XRminusD.A*origin;
originInterior = all(slackAtOrigin > 1e-8);

fprintf('Number of X_R^i minus D^i inequalities: %d\n', size(XRminusD.A,1));
fprintf('Origin belongs to X_R^i minus D^i: %s\n', logicalText(originContained));
fprintf('Origin is strictly inside X_R^i minus D^i: %s\n', ...
    logicalText(originInterior));
fprintf('Minimum inequality slack at the origin: %.6g\n', min(slackAtOrigin));

%% Plot the nested sets
figure('Color', 'w', 'Name', 'Step 5 - Pontryagin difference');
X.plot('color', [0.80 0.80 0.80], 'linewidth', 1.0);
hold on
Xrobust.plot('color', [0.49 0.18 0.56], 'linewidth', 1.4);
XRminusD.plot('color', [0 0.45 0.74], 'linewidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on
xlabel('e_p [m]')
ylabel('e_v [m/s]')
title('Static RMPC robustness constraint X_R^i minus D^i')
legend('State constraint X^i', 'Robust invariant set X_R^i', ...
    'RMPC constraint X_R^i minus D^i', 'Origin', 'Location', 'best')

save('rmpc_constraint_step5.mat', 'XRminusD', 'robustMargin', ...
    'originContained', 'originInterior', 'slackAtOrigin');

function text = logicalText(value)
    if value
        text = 'true';
    else
        text = 'false';
    end
end
