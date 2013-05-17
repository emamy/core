classdef JKerMorExport < export.JaRMoSExport
% AppExport: Export class for Android ROMSim App
%
%
% @author Daniel Wirtz @date 2011-08-02
%

%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    methods
        function this = JKerMorExport
            this.ModelType = 'JKerMor';
        end
    end
    
    methods(Access=protected)
        
        function [sourcebase, sourcefiles] = typeSpecificExport(this, f, model, folder)
            rm = model;
            
            fprintf(f,'\t<kermor_model>\n');
            fprintf(f,'\t\t<T>%17.17f</T>\n',rm.T);
            fprintf(f,'\t\t<dt>%17.17f</dt>\n',rm.dt);
            
            % Get system dimension from initial value
            mu = [];
            if ~isempty(rm.ParamSamples)
                mu = rm.ParamSamples(:,1);
            end
            dim = size(rm.System.x0.evaluate(mu),1);
            fprintf(f,'\t<dim>%d</dim>\n',dim);
            
            %% Export model data
            this.saveRealMatrix(rm.V,'V.bin',folder);
            this.saveRealMatrix(rm.W,'W.bin',folder);
            %export.AppExport.saveRealMatrix(rm.W,'G.bin',folder);
                        
            % ODE solver type
            stype = 'explicit';
            if isa(rm.ODESolver,'solvers.IImplSolver')
                stype = 'implicit';
            end 
            fprintf(f,'\t<solvertype>%s</solvertype>\n',stype);
            fprintf(f,'\t<outputtodof>TimoOutToDoF</outputtodof>\n');
            sources = {'TimoOutToDoF'};
            
            %% Export system data
            % Parameters
            fprintf(f,'\t<parameters>\n');
            p = rm.FullModel.System.Params;
            if ~isempty(p)
                %pvals = zeros(length(p),2);
                for k=1:length(p)
                    fprintf(f,'\t\t<param name="%s" min="%17.17f" max="%17.17f" label="%s"/>\n',...
                        p(k).Name,p(k).MinVal,p(k).MaxVal,p(k).Name);
                end
                %export.AppExport.saveRealMatrix(pvals, 'paramvalues.bin', folder);
            end
            fprintf(f,'\t</parameters>\n');
            
            % Kernel expansion
            s = rm.System;
            cf = s.f;
            fprintf(f,'\t<corefun type="%s">\n',class(cf));
            if isa(cf,'kernels.KernelExpansion')
                this.saveRealMatrix(cf.Ma,'Ma.bin',folder);
                this.saveRealMatrix(cf.Centers.xi,'xi.bin',folder);
                exportKernel(cf.Kernel,'kernel.bin',folder);
                fprintf(f,'\t\t<statekernel>%s</statekernel>\n',class(cf.Kernel));
                if isa(cf,'kernels.ParamTimeKernelExpansion')
                    this.saveRealVector(cf.Centers.ti,'ti.bin',folder);
                    exportKernel(cf.TimeKernel,'timekernel.bin',folder);
                    fprintf(f,'\t\t<timekernel>%s</timekernel>\n',class(cf.TimeKernel));
                    
                    this.saveRealMatrix(cf.Centers.mui,'mui.bin',folder);
                    exportKernel(cf.ParamKernel,'paramkernel.bin',folder);
                    fprintf(f,'\t\t<paramkernel>%s</paramkernel>\n',class(cf.ParamKernel));
                end
            elseif isa(cf, 'dscomponents.LinearCoreFun')
                this.saveRealMatrix(cf.A,'A.bin',folder);
            elseif isa(cf, 'dscomponents.AffLinCoreFun')
                if isempty(cf.CoeffClass)
                    error('AffLinCoreFuns must have the CoeffClass value set for export.');
                end
                % Set path to IAffineCoefficients class to compile
                fprintf(f,'\t\t<coeffclass>%s</coeffclass>\n',cf.CoeffClass);
                sources{end+1} = cf.CoeffClass;
                this.saveRealMatrix(cf.AffParamMatrix.Matrices,'A.bin',folder);
            else
                error('System function type unknown for export.');
            end
            fprintf(f,'\t</corefun>\n');
            
            % Input
            if rm.System.InputCount > 0 && ~isempty(s.B)
                fprintf(f,'\t<inputconv type="%s">\n',class(s.B));
                if isa(s.B,'dscomponents.LinearInputConv')
                    this.saveRealMatrix(s.B.B,'B.bin',folder);
                elseif isa(s.B,'dscomponents.AffLinInputConv')
                    if isempty(s.B.CoeffClass)
                        error('AffLinInputConv instances must have the CoeffClass value set for export.');
                    end
                    % Set path to IAffineCoefficients class to compile
                    fprintf(f,'\t\t<coeffclass>%s</coeffclass>\n',s.B.CoeffClass);
                    sources{end+1} = s.B.CoeffClass;
                    this.saveRealMatrix(s.B.Matrices,'B.bin',folder);
                end
                fprintf(f,'\t</inputconv>\n');
                sources{end+1} = 'Inputs';
            end
            
            % Mass matrix
            if ~isempty(s.M)
                fprintf(f,'\t<massmatrix type="%s">\n',class(s.M));
                if isa(s.M,'dscomponents.ConstMassMatrix')
                    this.saveRealMatrix(s.M.M,'M.bin',folder);
                elseif isa(s.M,'dscomponents.AffLinMassMatrix')
                    if isempty(s.M.CoeffClass)
                        error('AffLinMassMatrix instances must have the CoeffClass value set for export.');
                    end
                    % Set path to IAffineCoefficients class to compile
                    fprintf(f,'\t\t<coeffclass>%s</coeffclass>\n',s.M.CoeffClass);
                    sources{end+1} = s.M.CoeffClass;
                    this.saveRealMatrix(s.M.Matrices,'M.bin',folder);
                end
                fprintf(f,'\t</massmatrix>\n');
            end
            
            % Output
            if ~isempty(s.C)
                if isa(s.C,'dscomponents.LinearOutputConv')
                    C = s.C.C;
                    if isscalar(C)
                        C = eye(dim);
                    end
                    this.saveRealMatrix(C,'C.bin',folder);
                elseif isa(s.C,'dscomponents.AffLinOutputConv')
                    error('Not yet implemented.');
                end
                fprintf(f,'\t<outputconvtype>%s</outputconvtype>\n',class(s.C));
            end
            
            % Initial value
            if isa(s.x0,'dscomponents.ConstInitialValue')
                this.saveRealVector(s.x0.x0,'x0.bin',folder);
            elseif isa(s.B,'dscomponents.AffineInitialValue')
                error('Not yet implemented.');
            end
            fprintf(f,'\t<initialvaluetype>%s</initialvaluetype>\n',class(s.x0));
            fprintf(f,'\t</kermor_model>\n');
            
            % Assign source files and base here (from settings; here as
            % JaRMoSExport does not "know" the KerMor etc)
            sourcebase = KerMor.App.JKerMorSourceDirectory;
            sourcefiles = sources;
            
            function exportKernel(k, file, folder)
                if isa(k,'kernels.GaussKernel')
                    this.saveRealVector(k.Gamma,file,folder);
                elseif isa(k,'kernels.LinearKernel')
                    % do nothing.
                end
            end
        end
    end
end