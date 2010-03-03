function odefun = gen_odefun(system, mu, inputidx)
    % System without inputs
    if isempty(inputidx) || isempty(system.inputs)
        u = @(dummy)0;
    else
        % generates the ode function for given parameter and input function
        u = system.inputs{inputidx};
    end
    odefun = @(t,x)(system.f(x,t,mu) + system.B(t,mu)*u(t));
end

