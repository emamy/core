function f = approx_f(model_data, triples)
%APPROX_F Summary of this function goes here
%   Detailed explanation goes here

% Support vectors in format [t; x; mu] each column
sv = model_data.suppvect;
% Contains struct with fields ai,b,svidx
data = model_data.approx_data;

dims = size(data,2);

f = zeros(dims,1);

% Compute kernel matrix (only to be computed once for all regressions)
K = model.kernels.combined(sv,triples,dims,...
                    model.kernels.time,...
                    model.kernels.space,...
                    model.kernels.param);
                
for idx = 1:dims
    f(idx,:) = data(idx).ai'*K + b;
end

end

