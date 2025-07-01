classdef ContamxJr_FuncTests < matlab.unittest.TestCase

    properties(Constant)
        verbose = false;             % Set to true for detailed output
        codePath = '..//functions';
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setupOnce(testCase)
            addpath(ContamxJr_FuncTests.codePath);
            if ContamxJr_FuncTests.verbose
                fprintf('Running TestClassSettup - setupOnce() %s\n', ContamxJr_FuncTests.codePath);
            end
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function test_DensityCalc(testCase)
            relTol = 1e-6;

            if ContamxJr_FuncTests.verbose
                fprintf('\nRunning test_DensityCalc()\n')
            end
            Pambt = 101325.0;
            T = [20., 0., 30.];
            expVal = [1.204097, 1.292261, 1.164378];
            actVal = zeros(length(T));
            for i = 1:length(T)
                actVal(i) = Density(Pambt, T(i));
                if ContamxJr_FuncTests.verbose
                    fprintf('\texpVal= %.7f actVal= %.7f diff= %+.5e eps= %g\n', ...
                        expVal(i), actVal(i), (actVal(i) - expVal(i))/expVal(i), relTol);
                end
                testCase.assertEqual(actVal(i), expVal(i), 'RelTol', relTol);
            end
        end % end testDensityCalc()

        function test_FlowCoefCalc(testCase)
            relTol = 1e-6;
            if ContamxJr_FuncTests.verbose
                fprintf('\nRunning test_FlowCoefCalc()\n')
            end

            Paths = struct( ... 
                ... Paths struct is abbreviated version of that used in
                ... ContamxJr examples.
                'area', {       0.01        0.04          0.01          0.04  },...    % Opening area m2/m2
                'Cd',   {        0.6         0.6           0.6           0.6  },...    % Discharge coefficients
                'expt' ,{        0.5         0.5          0.65          0.65  },...    % Use with PL_LEAK
                'Pref' ,{        4.0         4.0          10.0          10.0  },...    % Use with PL_LEAK
                'ReTr' ,{       100.        100.           30.           30.  },...    % Transition Reynolds number
                'Clam', {8.124330e-6 6.499464e-5  1.315118E-06  7.640431E-06  },...    % Laminar flow coefficient
                'Cturb',{ 0.00848528  0.03394113   0.006007119   0.024028477  },...    % Turbulent flow coefficient
                'elem' ,{  'PL_ORFC'   'PL_ORFC'    'PL_LEAK3'    'PL_LEAK3'  }...
                );

            N = length(Paths);
            coeffs = zeros(N,N);    % [Cturb, Clam]
            for i = 1:N
                [coeffs(i,1),coeffs(i,2)] = setFlowCoef(Paths(i), Paths(i).elem);
                expVal = Paths(i).Cturb;
                actVal = coeffs(i,1);
                testCase.assertEqual(actVal, expVal, 'RelTol', relTol);
                if ContamxJr_FuncTests.verbose
                    fprintf('\t%s\n', Paths(i).elem);
                    fprintf('\tCturb: expVal= %.7e actVal= %.7e diff= %+.5e eps= %g\n', ...
                        expVal, actVal, (actVal - expVal)/expVal, relTol);
                end
                expVal = Paths(i).Clam;
                actVal = coeffs(i,2);
                testCase.assertEqual(actVal, expVal, 'RelTol', relTol);
                if ContamxJr_FuncTests.verbose
                    fprintf('\tClam:  expVal= %.7e actVal= %.7e diff= %+.5e eps= %g\n', ...
                        expVal, actVal, (actVal - expVal)/expVal, relTol);
                end
            end
        end % test_FlowCoefCalc()

    end % end methods(Test)

end % end class ContamxJr_FuncTests