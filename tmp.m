function tmp

    mu(1)=pi/2; % rotation angle speed
    mu(2)=-0.02; % angle offset/difference
    t = 0:.1:40;
    
    [to,yo] = ode45(@odefun,t,[0;1]);
    
    plot(yo(:,1),yo(:,2),'r');

    function y = odefun(t,x)
            a = mu(1);
            b = a+mu(2);
            %b = a+mu(2)+sin(t)/2;
            A = [cos(a) -sin(b); 
                 sin(a) cos(b)];
            y = A*x;
    end

end