classdef Test_NMFManopt < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 13-May-2024 14:24:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 24.1.0.2578822 (R2024a) Update 2 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mladni.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        
        function test_self_test(this)
            this.testObj_.self_test()
        end

        function test_original(this)
            call(this.testObj_);
        end

        function test_voxel_partitions(this)

            import mladni.NMFManopt

            m = 1000;
            mhalf = floor(m/2);
            mqtr = floor(m/4);
            n = 200;

            % X ~ [0.25, 0.5, 0.25, 1]
            X = ones(m, n);
            X(1:mhalf,:) = 0.5;
            X(1:mqtr,:) = 0.25;
            X(mhalf+1:mhalf+mqtr,:) = 0.25;
            X = X + 0.05*randn(m, n);
            X = X / norm(X);
            
            r = 2;
            W0 = NMFManopt.NNDSVD(X, r, 3);
            obj = NMFManopt( ...
                single(X), r, W0, @grassmannfactory, ...
                do_plot=true, do_AD=false, epsilon0=0, ...
                solver_name="stochasticgradient");
            tic
            call(obj);
            toc

            disp(obj)
        end

        function test_pca_voxel_partitions(this)

            import mladni.NMFManopt

            m = 1000;
            mhalf = floor(m/2);
            mqtr = floor(m/4);
            n = 200;

            % X ~ [0.25, 0.5, 0.25, 1]
            X = ones(m, n);
            X(1:mhalf,:) = 0.5;
            X(1:mqtr,:) = 0.25;
            X(mhalf+1:mhalf+mqtr,:) = 0.25;
            X = X + 0.05*randn(m, n);
            X = X / norm(X);
            
            r = 2;
            W0 = NMFManopt.NNDSVD(X, r, 3);
            obj = NMFManopt( ...
                single(X), r, W0, @stiefelfactory, ...
                do_plot=true, do_AD=false, epsilon0=0, ...
                solver_name="stochasticgradient");
            call(obj);

            disp(obj)            
        end

        function test_pca_original(this)
            obj = mladni.PCAManopt();
            call(obj);
        end

        function test_pca_rand_X(this)
            X = mladni.PCAManopt.rand_X();
            disp(X)
            r = 1;
            obj = mladni.PCAManopt(X, r);
            call(obj)
        end

        function test_pca_rand_X_opts(this)
            X = mladni.PCAManopt.rand_X(m=3, n=30, outliers=10);
            disp(X)
            r = 2;
            obj = mladni.PCAManopt(X, r, n_iterations=12, reduction=0.618, do_plot=true);
            call(obj)
        end

        function test_pca_onehot(this)
            n = 10;
            X1 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0] * ones(1, n);
            X1 = X1 + 0.05*rand(size(X1, 1), n);
            X2 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 1] * ones(1, n);
            X2 = X2 + 0.05*rand(size(X2, 1), n);
            X = X1 + X2;
            X = X / norm(X);
            r = 2;
            obj = mladni.PCAManopt(X, r, n_iterations=6, do_plot=false);
            call(obj);
        end
    end
    
    methods (TestClassSetup)
        function setupNMFManopt(this)
            import mladni.*
            this.testObj_ = NMFManopt();
        end
    end
    
    methods (TestMethodSetup)
        function setupNMFManoptTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
