function model = test_basemodel
%BASE_TESTMODEL Creates the basic test model according to testsettings

ts = testsettings;

model = base_model;

model.times = ts.times;

model.sampling.mode = ts.mode;
model.sampling.samples = ts.samples;

% Set system here, too to have running testmodel external to testKerMor
model.system = test_dynsys;

model.approx.mode = 'scalar_svr';
model.approx.scalar_svr.eps = 10;
model.approx.scalar_svr.C = 100000;
end

