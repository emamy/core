classdef FileMatrix < matlab.unittest.TestCase
    % Tests for the data.FileMatrix class
    
    methods(TestClassSetup)
        function setup(~)
            warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym');
        end
    end
    
    methods(TestClassTeardown)
        function teardown(~)
            warning('on','MATLAB:eigs:TooManyRequestedEigsForRealSym');
        end
    end
    
    methods(Test)
        function General(this)
            B = rand(99,100);
            A = data.FileMatrix(99,100,'BlockSize',data.FileMatrix.blockSizeOf(B)/5);
            % col-wise setting (fast)
            for k=1:100
                A(:,k) = B(:,k);
            end
            this.verifyEqual(A(:,:),B,'Column-wise assignment mismatch');
            
            % Row-wise setting
            for k=1:10
                A(:,k) = B(:,k);
            end
            this.verifyEqual(A(:,:),B,'Row-wise assignment mismatch');
            
            this.verifyEqual(A.toMemoryMatrix,B,'toMemoryMatrix mismatch');
            
            % Direct constructor test
            A2 = data.FileMatrix(B);
            this.verifyEqual(A(:,:),B,'Direct constructor test');
            
            % IsEqual test
            this.verifyTrue(A == B,'isEqual test case');
            B2 = B; B2(1,:) = -B2(1,:);
            this.verifyTrue(A ~= B2,'isEqual test case');
            
            % Transpose test
            A.InPlaceTranspose = false;
            At = A';
            this.verifyEqual(At(:,:),B','Transpose test');
            A.InPlaceTranspose = true;
            At = A';
            this.verifyEqual(At(:,:),B','Transpose test 2');
            
            % Bounding box test
            [bm, bM] = Utils.getBoundingBox(B);
            [am, aM] = A.getColBoundingBox;
            this.verifyEqual(bm,am,'Bounding box test Min');
            this.verifyEqual(bM,aM,'Bounding box test Max');
        end
        
        function transposeTimes(this)
            [A1,B1] = this.getTestPair(80,100,5);
            [A2,B2] = this.getTestPair(80,100,5);
            A3 = A1.transposedTimes(A2);
            this.verifyLessThan(norm(B1'*B2 - A3.toMemoryMatrix,'fro'),sqrt(eps));
        end
        
        function SVD(this)
            [A,B] = this.getTestPair(99,100,5);
            % SVD test
            p = sqrt(eps);
            [u,s,v] = svd(B,'econ');
            [U,S,V] = A.getSVD;
            V = V.toMemoryMatrix;
            this.verifyLessThan(norm(U*S*V-B,'fro'),p);
            [U5,S5,V5] = A.getSVD(5);
            V5 = V5.toMemoryMatrix;
            this.verifyLessThan(norm(abs(V)-abs(v'),'fro'),p);
            this.verifyLessThan(norm(abs(U)-abs(u),'fro'),p);
            this.verifyLessThan(norm(diag(S)-diag(s)),p);
            
            this.verifyLessThan(norm(abs(V5)-abs(v(:,1:5)'),'fro'),p);
            this.verifyLessThan(norm(abs(U5)-abs(u(:,1:5)),'fro'),p);
            this.verifyLessThan(norm(diag(S5)-diag(s(1:5,1:5))),p);
            
            % Exclude test
            o = general.Orthonormalizer;
            exclu = o.orthonormalize(rand(size(B,1),5));
            [U,~,~] = A.getSVD(10,exclu);
            this.verifyLessThan(norm(exclu'*U),1e-12);
        end
        
        function SaveLoad(this)
            [A,B] = this.getTestPair(100,10000,1);
            % Load/save
            save Atmp A;
            clear A;
            load Atmp;
            this.verifyTrue(A == B,'Loading FileMatrix');
            rmdir(A.DataDirectory,'s');
            delete Atmp.mat;
            
            % Multiblock save/load
            [A,B] = this.getTestPair(100,10000,10);
            
            % Load/save
            save Atmp A;
            clear A;
            load Atmp;
            this.verifyTrue(A == B,'Loading FileMatrix');
            rmdir(A.DataDirectory,'s');
            delete Atmp.mat;
        end
        
        function SumPower(this)
            [A,B] = this.getTestPair(100,10000,40);
            
            % Sum test
            bs = sum(B,2);
            as = sum(A,2);
            this.verifyEqual(sum(B,1),sum(A,1),'Sum operation');
            this.verifyLessThan(abs((as-bs) ./ bs),sqrt(eps),'Relative comparison');
            
            % Power test
            this.verifyTrue(A.^2 == B.^2,'Power operation');
        end
        
        function ScalarMult(this)
            % @todo need implementation for >1 nBlock matrices and scalar values
            % @todo need separate test for each overridden operator
            
            [A,a] = this.getTestPair(400,400);
            B = 4*A;
            this.verifyEqual(4*a,B);
            C = A*1;
            this.verifyEqual(a,C);
            
            B = 4.*A;
            this.verifyEqual(4*a,B);
            C = A.*1;
            this.verifyEqual(a,C);
            
            [A,a] = this.getTestPair(400,400,4);
            B = 4*A;
            this.verifyTrue(4*a == B);
            C = A*1;
            this.verifyTrue(a == C);
            
            B = 4.*A;
            this.verifyTrue(4*a == B);
            C = A.*1;
            this.verifyTrue(a == C);
        end
        
        function Times_MTimes(this)
            % @todo need implementation for >1 nBlock matrices and scalar values
            % @todo need separate test for each overridden operator
            
            % No blocks
            [A,B] = this.getTestPair(100,1000);
            AB = A.*B;
            this.verifyEqual(AB,B.*B);
            AB = B.*A;
            this.verifyEqual(AB,B.*B);
            
            % With blocks
            [A,B] = this.getTestPair(100,1000,4);
            AB = A.*B;
            this.verifyTrue(AB == B.*B);
            AB = B.*A;
            this.verifyTrue(AB == B.*B);
            
            % Vector Multiplication test
            v = rand(size(A,2),1);
            this.verifyLessThan(norm(A*v-B*v),1e-10);
            v = rand(size(A,2),100);
            this.verifyLessThan(Norm.L2(A*v-B*v),1e-10);
            v = rand(1,size(A,1));
            this.verifyLessThan(norm(v*A-v*B),1e-10);
            v = rand(100,size(A,1));
            this.verifyLessThan(Norm.L2(v*A-v*B),1e-10);
            
            % Left / right multiplication test
            R = rand(1000,200);
            L = rand(2000,50);
            % With no blocks
            [A, Amat] = this.getTestPair(50,1000);
            LA = L*A;
            AR = A*R;
            this.verifyEqual(LA,L*Amat,'Left / right multiplication test with no blocks');
            this.verifyEqual(AR,Amat*R,'Left / right multiplication test with no blocks');
            % With blocks
            [A, Amat] = this.getTestPair(50,1000,4);
            LA = L*A;
            AR = A*R;
            % Left multiplication is exact as whole blocks can be used
            this.verifyTrue(LA == L*Amat,'Left multiplication is exact as whole blocks can be used');
            % Right multiplication requires accumulation of values -> rounding errors
            this.verifyLessThan(norm(AR - Amat*R),1e-11);
        end
        
        function TwoFileMatrix_MTimes(this)
            %% Two FileMatrices
            % No blocks
            [A,B] = this.getTestPair(100,1000);
            [C,D] = this.getTestPair(1000,100);
            AC = A*C;
            this.verifyTrue(AC == B*D);
            CA = C*A;
            this.verifyTrue(CA == D*B);
            
            % With blocks
            [A,B] = this.getTestPair(100,1000,4);
            [C,D] = this.getTestPair(1000,100,5);
            AC = A*C;
            this.verifyLessThan(norm(B*D - AC.toMemoryMatrix,'fro'),1e-10);
            CA = C*A;
            this.verifyLessThan(norm(D*A - CA.toMemoryMatrix,'fro'),1e-10);
        end
        
        function SpeedSVDTransp(this)
            % Test results for 100000x100 matrix:
            % BlockSVD: Computing truncated 50-SVD on 100x1000000 matrix (3 blocks)...
            % Direct time: 2475.53, Transposed time: 756.023, transposed SVD time: 661.171
            %
            % Thus: Auto-transpose for matrices with nBlocks > 1
            A = this.getTestPair(100,1000,5);
            ti = tic;
            [u,s,v] = A.getSVD(5);%#ok
            t(1) = toc(ti);
            
            ti = tic;
            B = A';
            t(2) = toc(ti);
            
            ti = tic;
            [v2,s2,u2] = B.getSVD(5);%#ok
            t(3) = toc(ti);
            t(2) = t(2) + t(3);
            
            fprintf('Direct SVD time: %g, Transposition time: %g, Transposed SVD time: %g\n',t);
        end
        
        function PartialSVD(this)
            % Tests the selective SVD/POD block-wise algorithm
            d = 100;
            A = this.getTestPair(d,1000,3);
            B = this.getTestPair(d,1000,1);
            idx = [1 d-1 floor(rand(1,10)*d)]+1;
            for k = 1:length(idx)
                [u,s] = B.getSVD(100,[],1:idx(k));
                this.verifyEqual(size(u,1),size(s,1),'Size mismatch!');
                [u,s] = A.getSVD(100,[],1:idx(k));
                this.verifyEqual(size(u,1),size(s,1),'Size mismatch!');
            end
        end
    end
    
    methods(Access=private)
        function [A,B] = getTestPair(~,n,m,nb)
            if nargin < 4
                nb = 1;
            end
            B = rand(n,m);
            A = data.FileMatrix(n,m,'BlockSize',data.FileMatrix.blockSizeOf(B)/nb);
            A(:,:) = B;
            % Tedious formulation as use of direct operators inside the defining class of
            % overloads will not work.
            % See http://www.mathworks.com/help/matlab/matlab_oop/indexed-reference-and-assignment.html#br09nsm
            %             A.subsasgn(struct('type',{'()'},'subs',{{':', ':'}}),B);
        end
    end
    
    methods(Static)
        function res = test_FileMatrix
            % Legacy function so that MUnit also starts 
            res = runtests('data.tests.FileMatrix');
            res = all([res(:).Passed]);
        end
    end
    
end

