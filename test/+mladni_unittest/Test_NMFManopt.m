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

        function test_X(this)

            import mladni.NMFManopt

            ld = load( ...
                fullfile(getenv("ADNI_HOME"), "NMF_FDG", "baseline_cn", "NumBases24", "X.mat"));
            X = ld.X;  % 228483 x 269
            X = X / norm(X, "fro");
            X(X < 1e-14) = 1e-14;
            
            for r = [2,6,8,10,12,14]
                W0 = NMFManopt.NNDSVD(X, r, 3);
                obj = NMFManopt( ...
                    X, r, W0, @mgrassmannfactory, ...
                    maxiter=1e5, ...
                    stepsize_lambda=1e-3, ...
                    do_force_nonneg=true, do_plot=false, do_AD=false, epsilon0=0, ...
                    solver_name="stochasticgradient");
                tic
                call(obj);
                toc
                disp(obj)
                fqfp = fullfile(getenv("ADNI_HOME"), "NMF_FDG", "baseline_cn", "NumBases"+r, ...
                    strrep(stackstr()+"_maxiter"+obj.maxiter, ".", "p"));
                ensuredir(myfileparts(fqfp));
                save(fqfp+".mat", "obj");

                ifc = mlfourd.ImagingFormatContext2( ...
                    fullfile(getenv("ADNI_HOME"), "VolBin", "mask.nii.gz"));
                ifc.img = double(ifc.img);
                img = ifc.img;
                bin = logical(img);
                for rho = 1:(r-1)
                    ifc.img = cat(4, ifc.img, img);  % 3D -> 4D
                end
                for rho = 1:r
                    img_ = img;
                    img_(bin) = obj.W(:, rho);
                    ifc.img(:, :, :, rho) = img_;
                end
                ifc.fqfp = fqfp;
                ifc.save();
                % ifc.view()

                % reconstruction
                ifc.fileprefix = ifc.fileprefix + "_recon";
                ifc.img = img;
                X_hat = obj.W * obj.H;
                for nu = 1:size(X_hat, 2)
                    img_ = img;
                    img_(bin) = X_hat(:, nu);
                    ifc.img(:,:,:,nu) = img_;
                end
                ifc.save()
                % ifc.view()
            end
        end

        function test_X_3mm(this)

            import mladni.NMFManopt

            ld = load( ...
                fullfile(getenv("ADNI_HOME"), "NMF_FDG", "baseline_cn", "NumBases24", "X_3mm.mat"));
            X = ld.X_3mm;  % 67019 x 269
            % X = X(:, 1:10);
            X = X / norm(X, "fro");
            X(X < 1e-14) = 1e-14;
            
            r = 6;
            W0 = NMFManopt.NNDSVD(X, r, 3);
            obj = NMFManopt( ...
                X, r, W0, @grassmannfactory, ...
                maxiter=1e5, ...
                stepsize_lambda=1, ...
                batchsize=26, ...
                do_force_nonneg=true, do_plot=false, do_AD=false, epsilon0=0, ...
                solver_name="stochasticgradient");
            tic
            call(obj);
            toc
            disp(obj)
            fqfp = fullfile(getenv("ADNI_HOME"), "NMF_FDG", "baseline_cn", "NumBases"+r, ...
                strrep(stackstr()+"_maxiter"+obj.maxiter, ".", "p"));
            ensuredir(myfileparts(fqfp));
            save(fqfp+".mat", "obj");

            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(getenv("ADNI_HOME"), "VolBin", "mask_3mm.nii.gz"));
            ifc.img = double(ifc.img);
            img = ifc.img;
            bin = logical(img);            
            for rho = 1:(r-1)
                ifc.img = cat(4, ifc.img, img);  % 3D -> 4D
            end
            for rho = 1:r
                img_ = img;
                img_(bin) = obj.W(:, rho);
                ifc.img(:, :, :, rho) = img_;
            end
            ifc.fqfp = fqfp;
            ifc.save();
            ifc.view()
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
            X = X + 0.02*randn(m, n);
            % X = X / norm(X, "fro");
            X(X < 1e-16) = 1e-16;
            
            r = 2;
            W0 = NMFManopt.NNDSVD(X, r, 3);
            obj = NMFManopt( ...
                X, r, W0, @grassmannfactory, ...
                batchsize=20, ...
                maxiter=2e4, ...
                stepsize_lambda=1, ...
                do_force_nonneg=true, do_plot=true, do_AD=false, epsilon0=0, ...
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
            X = X + 0.02*randn(m, n);
            % X = X / norm(X, "fro");
            X(X < 1e-16) = 1e-16;
            
            r = 2;
            W0 = NMFManopt.NNDSVD(X, r, 3);
            obj = NMFManopt( ...
                X, r, W0, @stiefelfactory, ...
                do_force_nonneg=false, do_plot=true, do_AD=false, epsilon0=0, ...
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
