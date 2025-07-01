classdef AirnetTests < matlab.unittest.TestCase
    properties(Constant)
        verbose = false;             % Set to true for detailed output
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setupOnce(testCase)
            %addpath('..\\');
            %addpath('..\\funtions')
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function test_AirnetPl1(testCase)
            if AirnetTests.verbose
                fprintf('\nRunning %s\t%s\n',coder.mfunctionname, pwd());
            end
            relTol = 1e-5;

            mdot_exp = [
                0.00903301,  -0.00903301,  0.00285649,  -0.00285649,...
                0.000903301, -0.000903301, 0.000362448, -0.000362448...
                0.000362178, -0.000362178, 2.03953E-05, -2.03954E-05...
                2.01081E-05, -2.01081E-05, 1.91505E-05, -1.91505E-05];
            [~, mdot, ~] = airnet_pl1();
            NumAct = length(mdot);
            NumExp = length(mdot_exp);
            if AirnetTests.verbose
                fprintf('NumAct= %d, NumExp= %d\n', NumAct, NumExp);
            end
            testCase.verifyEqual(NumAct, NumExp);
            for i=1:NumExp
                testCase.verifyEqual(mdot(i), mdot_exp(i), 'RelTol', relTol);
            end
        end % end test_AirnetPl1()


        function test_AirnetPl2(testCase)
            if AirnetTests.verbose
                fprintf('\nRunning %s\t%s\n',coder.mfunctionname, pwd());
            end
            relTol = [1e-5, 1e-5, 5e-5, 1e-5, 1e-5, 1e-5, 1e-5];

            mdot_exp = [
                5.375717E-04, 6.567484E-04, 6.583718E-04, 6.583881E-04,...
                6.583882E-04, 6.583882E-04, 6.583882E-04
                ];

            mult = logspace(0,6,7);
            Ncases = length(mult);
            testCase.assertLength(mdot_exp, Ncases);
            for i=1:length(mult)
                [~,mdot, ~] = airnet_pl2(0, mult(i) );
                abs_mdot = abs(mdot);
                for j=1:length(abs_mdot)
                    if AirnetTests.verbose
                        relErr = (abs_mdot(j) - mdot_exp(i))/mdot_exp(i);
                        fprintf('Mult= %7d: absMdot= %g, MdotExp= %g, err=%g\n', ...
                            mult(i), abs_mdot(j), mdot_exp(i), relErr);
                    end
                    testCase.verifyEqual( abs_mdot(j), mdot_exp(i), 'RelTol', relTol(i));
                end
            end
        
        end % end test_AirnetPl2


        function test_AirnetPl3(testCase)
            if AirnetTests.verbose
                fprintf('\nRunning %s()\t%s\n', coder.mfunctionname, pwd());
            end
            relTol = 1e-5;

            mdot_inlet  =  6.11010E-02;
            mdot_outlet = -6.11010E-02;
            dPexp = 100.0;
            [~, mdot, dp] = airnet_pl3();
            NumAct = length(mdot);
            if AirnetTests.verbose
                fprintf('NumAct= %d, NumExp= %d\n', NumAct, 20);
            end
            testCase.verifyEqual(NumAct, 20);
            testCase.verifyEqual(mdot(9), mdot_inlet, 'RelTol', relTol);
            testCase.verifyEqual(mdot(15), mdot_outlet, 'RelTol', relTol);

            % Each pathway in series through the flow network should have 
            % the same drop in pressure in summation.
            dPact = zeros(1,7);
            dPact(1)= abs(dp(9))+abs(dp(1))+abs(dp(3))+abs(dp(2))+abs(dp(15));
            dPact(2)= abs(dp(9))+abs(dp(4))+abs(dp(6))+abs(dp(5))+abs(dp(15));
            dPact(3)= abs(dp(9))+abs(dp(7))+abs(dp(6))+abs(dp(8))+abs(dp(15));
            dPact(4)= abs(dp(9))+ sum(abs(dp(10:15)));
            dPact(5)= abs(dp(9))+abs(dp(17))+abs(dp(16))+abs(dp(19))+abs(dp(15));
            dPact(6)= abs(dp(9))+abs(dp(17))+abs(dp(18))+abs(dp(19))+abs(dp(15));
            dPact(7)= abs(dp(9))+abs(dp(17))+abs(dp(20))+abs(dp(19))+abs(dp(15));

            NumAct = length(dPact);
            for i=1:NumAct
                testCase.verifyEqual(dPact(i), dPexp, 'RelTol', relTol);
            end
        end % end test_AirnetPl3()

        function test_AirnetStack1(testCase)
            if AirnetTests.verbose
                fprintf('\nRunning %s()\t%s\n', coder.mfunctionname, pwd());
            end
            relTol = 1e-5;

            stackDpExp = 8.644910;
            [~, ~, dP] = airnet_stack1(0);
            stackDpAct = abs(dP(1)) + abs(dP(2));
            testCase.verifyEqual(stackDpAct, stackDpExp, 'RelTol', relTol);

        end

    end % end methods(Test)

end % end class AirnetTests