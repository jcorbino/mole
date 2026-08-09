% Main script to run all the tests
%
% Runs every test, prints one PASS/FAIL line each, and ends with a summary.
% Works in both MATLAB and Octave.

clc
clear all;
close all;

% Description, file
mole_tests = { ...
    'Nullity test of Divergence operator', 'test1.m'; ...
    'Nullity test of Gradient operator',   'test2.m'; ...
    'Nullity test of Laplacian operator',  'test3.m'; ...
    'Energy test (Schrödinger equation)',  'test4.m'; ...
    'Accuracy test (Poisson equation)',    'test5.m'; ...
    'BC Ops consistency test',             'test6.m'};

mole_pass = false(size(mole_tests, 1), 1);
mole_log = cell(size(mole_tests, 1), 1);

% The mole_ prefix matters: run() executes each test in this same workspace, and
% the tests use short names like i, k and m, which would clobber anything plainer.
for mole_t = 1 : size(mole_tests, 1)
    try
        mole_log{mole_t} = evalc(sprintf('run(''%s'')', mole_tests{mole_t, 2}));

        % Each test reports by printing. Count it as passed only if it announced
        % PASSED and never announced FAILED -- both halves are needed. test1 to
        % test3 print PASSED even after printing FAILED, so the first half alone
        % would miss them; and test4 sends its FAILED to stderr, which MATLAB's
        % evalc does not capture, so there the absence of PASSED is the only
        % signal. Octave's evalc does capture stderr, and the rule holds either way.
        mole_pass(mole_t) = ~isempty(strfind(mole_log{mole_t}, 'PASSED')) && ...
                             isempty(strfind(mole_log{mole_t}, 'FAILED'));
    catch mole_err
        % A test that raised rather than reported
        mole_log{mole_t} = ['error: ' mole_err.message];
        mole_pass(mole_t) = false;
    end

    if mole_pass(mole_t)
        fprintf('  [PASS]  %s\n', mole_tests{mole_t, 1});
    else
        fprintf('  [FAIL]  %s\n', mole_tests{mole_t, 1});
    end
end

fprintf('\n%d of %d tests passed\n', sum(mole_pass), numel(mole_pass));

if all(mole_pass)
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('TESTS FAILED\n');
    % Echo what the failing tests said, since evalc swallowed it
    for mole_t = 1 : numel(mole_pass)
        if ~mole_pass(mole_t)
            fprintf('\n--- %s (%s) ---\n%s', mole_tests{mole_t, 1}, ...
                    mole_tests{mole_t, 2}, mole_log{mole_t});
        end
    end
end
