function [m,md] = testKerMor(full)
%testKerMor - Tests the KerMor framework according to testsettings

ts = testsettings;

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nRunning KerMor tests with settings:\n');
disp(ts);

if nargin == 0
    full = false;
end
%full = true;

% Full test on all availabe test models
if full
    modelTest(base_model);
    modelTest(test_model_lin);
    modelTest(test_model_lin_param);
    modelTest(test_model_lin_input);
    modelTest(test_model_lin_param_input);
    modelTest(test_model_nonlin);
    modelTest(test_model_nonlin_param);
    modelTest(test_model_nonlin_input);
    modelTest(test_model_nonlin_param_input);
end

% Default full test model
m = test_model_nonlin_param_input;

fprintf('Checking parameter sampling with ''grid''...');
m.mode = 'grid';
s = gen_param_samples(m);
fprintf('done\n');

fprintf('Checking parameter sampling with ''rand''...');
m.mode = 'rand';
m.samples = ts.samples;
ps = gen_param_samples(m);
fprintf('done\n');

fprintf('Checking snapshot generation...');
md = gen_snapshots(m, ps);
fprintf('done\n');

fprintf('Checking space reduction with ''POD''...\n');
m.reduction.mode = 'POD';

fprintf('Checking POD with ''eps''...');
m.reduction.POD.mode = 'eps';
m.reduction.POD.value = 1e-2;
md = gen_reduced_space(m,md);
fprintf('done\n');
disp(size(md.V));

fprintf('Checking POD with ''sign''...');
m.reduction.POD.mode = 'sign';
m.reduction.POD.value = .3;
md = gen_reduced_space(m,md);
fprintf('done\n');
disp(size(md.V));

fprintf('Checking POD with ''rel''...');
m.reduction.POD.mode = 'rel';
m.reduction.POD.value = .3;
md = gen_reduced_space(m,md);
fprintf('done\n');
disp(size(md.V));

fprintf('Checking POD with ''abs''...');
m.reduction.POD.mode = 'abs';
m.reduction.POD.value = round(ts.testdim/2);
md = gen_reduced_space(m,md);
fprintf('done\n');
disp(size(md.V));

fprintf('Checking scalar SVR ...');
testsvr;
fprintf('done\n');

fprintf('KerMor tests done.\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

    function modelTest(model)
        disp(['Checking model ' model.info.name ' ..........']);
        model.verbose = 2;
        sn = offline_phase(model);
        disp(sn);
    end

    function m = test_model_lin
        m = test_basemodel;
        m.info.name = 'Test Model linear';
        s = base_dynsystem;
        s.f = ts.flin;
        s.x0 = ts.x0;
        m.system = s;
    end

    function m = test_model_lin_param
        m = test_basemodel;
        m.info.name = 'Test Model linear with parameters';
        s = base_dynsystem;
        s.f = ts.flin_p;
        s.params = ts.params;
        s.x0 = ts.x0_p;
        m.system = s;
    end

    function m = test_model_lin_input
        m = test_basemodel;
        m.info.name = 'Test Model linear with inputs';
        s = base_dynsystem;
        s.inputs = ts.inputs;
        s.f = ts.flin;
        s.x0 = ts.x0;
        m.system = s;
    end

    function m = test_model_lin_param_input
        m = test_basemodel;
        m.info.name = 'Test Model linear with params and input';
        s = base_dynsystem;
        s.f = ts.flin_p;
        s.inputs = ts.inputs;
        s.params = ts.params;
        s.x0 = ts.x0_p;
        m.system = s;
    end

    function m = test_model_nonlin
        m = test_basemodel;
        m.info.name = 'Test Model nonlinear';
        s = base_dynsystem;
        s.f = ts.fnlin;
        s.x0 = ts.x0;
        m.system = s;
    end

    function m = test_model_nonlin_param
        m = test_basemodel;
        m.info.name = 'Test Model nonlinear with parameters';
        s = base_dynsystem;
        s.f = ts.fnlin_p;
        s.params = ts.params;
        s.x0 = ts.x0_p;
        m.system = s;
    end

    function m = test_model_nonlin_input
        m = test_basemodel;
        m.info.name = 'Test Model nonlinear with inputs';
        s = base_dynsystem;
        s.inputs = ts.inputs;
        s.f = ts.fnlin;
        s.x0 = ts.x0;
        m.system = s;
    end

    function m = test_model_nonlin_param_input
        m = test_basemodel;
        m.info.name = 'Test Model nonlinear with params and inputs';
        s = base_dynsystem;
        s.f = ts.fnlin_p;
        s.inputs = ts.inputs;
        s.params = ts.params;
        s.x0 = ts.x0_p;
        m.system = s;
    end
end