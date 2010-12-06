% step = .025;
% matsize = 300;
% intvl = step:step:1;
%
% n = length(intvl);
% by = zeros(2,n);
% for cnt = 1:n
%     somesparsemat = sparse(+(rand(matsize) < intvl(cnt)));
%     somefullmat = full(somesparsemat);
%     in = whos('some*mat');
%     by(1,cnt) = in(1).bytes;
%     by(2,cnt) = in(2).bytes;
% end
% plot(intvl,by(1,:),'r',intvl,by(2,:),'b');

% dx = .01;
% x = -5:dx:5;
% b = .5;
% y = exp(-b*x.^2);
% y2 = [0 diff(y)]/dx;
% y3 = [0 diff(y2)]/dx;
% f1 = -2*b*x.*exp(-b*x.^2);
% f2 = (4*b^2*x.^2-2*b).*exp(-b*x.^2);
% plot(x,0,'r',x,y,'g',x,y2,'b',x,y3,'m',x,f2,'m--',x,f1,'b--');

% dt = .01;
% % b = 18.5;
% % ep = 0.01;
% b = 3;
%
% x0 = sqrt(b/2);
% x = -1:dt:x0+2;
% fx = 2*x/b.*exp(-x.^2/b);
% %fx = abs(2/b - 4*x.^2/b^2);
% %fx2 = fx .* exp(-x.^2/b);
% plot(x,fx,'r');%,x,fx2,'b');
% ep = 0.001;
%
% x0 = sqrt(b*log(1/ep));
% x1 = x0-2:dt:x0+10;
% fx = exp(-x1.^2/b);
% dx = abs(gradient(fx,dt,dt));
% plot(x1,ep,'b',x1,fx,'r',x1,ep*exp(2*sqrt(log(1/ep)/b)),'b--',x1,dx,'g');
% axis tight;
% max(dx)
% pause;
%
% [x1,x2] = ndgrid(-5:dt:5);
% fx = exp(-x1.^2-x2.^2);
% [dx1,dx2] = gradient(fx,dt,dt);
% quiver(dx1,dx2);
% hold on;
% %contour(sqrt(x1.^2+x2.^2) -sqrt(.5),0);
% contour(sqrt(dx1.^2 + dx2.^2)-sqrt(.5),[0,0]);
% hold off;
% max(max(sqrt(dx1.^2 + dx2.^2)))
%
% figure;
% surf(sqrt(dx1.^2 + dx2.^2));

% [x1,x2,x3] = ndgrid(-5:dt:5);
% fx = exp(-x1.^2-x2.^2-x3.^2);
% [dx1,dx2,dx3] = gradient(fx,dt,dt,dt);
%quiver(dx1,dx2,dx3);
%hold on;
%contour(sqrt(x1.^2+x2.^2) -sqrt(.5),0);
%contour(sqrt(dx1.^2 + dx2.^2 + dx3.^2),[0,0]);
%hold off;
% max(max(max(sqrt(dx1.^2 + dx2.^2 + dx3.^2))))

function tmp

%% Presettings

dt = .2;
area = [-20 20];
gamma = 15;
% epsilon-threshold range
epsi = .005;
% rnd = false;
rnd = true;
% Initial current state point
xp = 3;
yp = 2;

% Positions
if rnd
    kNum = 50;
    xP = zeros(1,kNum);
    yP = xP;
    for cnt = 1:kNum
        xP(cnt) = area(1) + rand * (area(2)-area(1));
        yP(cnt) = area(1) + rand * (area(2)-area(1));
    end
else
    kNum = 4;
    xP = [0 .2 .2 .2 8];
    yP = [-2.5 .2 1 .5 5];
    xP = xP(1:kNum);
    yP = yP(1:kNum);
end
% Coefficients
if rnd
    %     c = (rand(1,kNum)-.5)*5;
    c = (rand(1,kNum))*5;
else
    c = [1 2 3 2 15];
    c = c(1:kNum);
end

%% Automatic
k = kernels.GaussKernel(gamma);

%r = area(2)-area(1);
range = area(1):dt:area(2);
[X,Y] = meshgrid(range,range);

% Plot data
vals = zeros(size(X,1),size(X,2),kNum);
for cnt = 1:kNum
    arg = sqrt((X-xP(cnt)).^2+(Y-yP(cnt)).^2);
    vals(:,:,cnt) = c(cnt)*k.evaluateScalar(arg);
end
Z = sum(vals,3);

% Distance matrix
D = getDistanceMat(xP,yP);
% Set large distances to zero

nullrad = sqrt(-gamma*log(epsi./abs(c)));
nullrad = repmat(nullrad,kNum,1);
F = k.evaluateScalar(D);
%W = D./nullrad;
W = 1-D./nullrad;
W(W<0)=epsi;
if ~rnd
    nullrad
    D
    W
    F
end
%W = ones(size(D));

% Angle matrix
%A = getAngleMat(xP,yP)

%% Postsettings
lipFun = @k.getImprovedLocalSecantLipschitz;

% Precomputations
% globEst = [Inf,-Inf];
% globEff = [Inf,-Inf];
% coarsest = -Inf;
% sharpest = Inf;
% for cnt = 1:numel(X)
%    x = X(cnt); y = Y(cnt);
%    [estMax, di, lip] = getEstGradient(x,y);
%    if globEst(1) > estMax
%       globEst(1) = estMax;
%    end
%    if globEst(2) < estMax
%       globEst(2) = estMax;
%    end
%    [effMax, xm, ym, dist] = getEffGradient(x,y);
%    effMax = abs(effMax);
%    if globEff(1) > effMax
%       globEff(1) = effMax;
%    end
%    if globEff(2) < effMax
%       globEff(2) = effMax;
%    end
%    rel = estMax/effMax;
%    if rel > coarsest
%        coarsest = rel;
%    end
%    if rel < sharpest
%        sharpest = rel;
%    end
%    if mod(cnt,200) == 0
%        fprintf('%d%%..',round(cnt/numel(X)*100));
%    end
% end
% fprintf('\n');
% coarsest = round(coarsest*100);
% sharpest = round(sharpest*100);

% Plotting
h = figure;
cm = datacursormode(h);
cm.Enable = 'on';
cm.SnapToDataVertex = 'on';
cm.DisplayStyle = 'window';
cm.Updatefcn = @pointSelected;


subplot(1,2,1);
%mesh(X,Y,Z);
surf(X,Y,Z);

shading interp
%grid off;
hidden off;
%lighting gouraud;
colormap hot;
axis equal;
hold on;
h = contour3(X,Y,Z,20,'black');

lineobj = cell.empty(1,0);

oldpos = [Inf Inf];
updatePlot(xp,yp);

    function updatePlot(xp,yp)
        
        try
            %% Altes zeug
            [estMax, di, lip] = getEstGradient(xp,yp);
            [effMax, xm, ym, dist, Zsecant] = getEffGradient(xp,yp);
            
            %% Session Jan
            
            xidiff = [xP - xp; yP-yp];
            Umat = xidiff'*xidiff;
            em = 0;
            for idx=1:length(xP)
                sameside = Umat(idx,:) > 0;
                [estMaxLoc, didummy, lipdummy] = getEstGradientJan(xp,yp,sameside);
                em(end+1) = abs(estMaxLoc);
            end
            em
            max(em)
            estMax
            effMax
            
            %% Altes zeug
            
            [dummy, yi] = min(abs(range - xp));
            [dummy, xi] = min(abs(range - yp));
            %         yi = find(range == xp,1);
            %         xi = find(range == yp,1);
            
            % Gradient vectors - real
            for idx=1:length(lineobj)
                lin = lineobj{idx};
                if ~isempty(lin) && ishandle(lin)
                    delete(lin);
                end
            end
            lineobj = {};
            
            
            
            
            subplot(1,2,1);
            
            % Local secant lipschitz estimations
            %         for idx = 1:kNum
            %             z1 = c(idx)*k.evaluateScalar(di(idx));
            %             z2 = z1 + c(idx)*lip(idx)*di(idx);
            %             lineobj{end+1} = plot3([xp xP(idx)],[yp yP(idx)],[z1 z2],'LineWidth',2);%#ok
            %         end
            
            % Plot single cone maxima
            centers = [xP; yP];
            fxi = c*k.evaluate(centers,centers);
            % Plot cone with current maximum secant
            [v,idx] = max(abs(fxi - Z(xi,yi)) ./ sum(sqrt([xP-xp;yP-yp].^2)));
            lineobj{end+1} = plot3(xP(idx),yP(idx),fxi(idx),'red.','MarkerSize',15);
            idx = setdiff(1:length(xP),idx);
            % Plot other peaks black
            lineobj{end+1} = plot3(xP(idx),yP(idx),fxi(idx),'black.','MarkerSize',15);
            
            % Max effective secant grad
            hlp = Z(xi,yi) + effMax * dist;
            lineobj{end+1} = plot3([xp X(xm,ym)],[yp Y(xm,ym)], [Z(xi,yi) hlp],'g','LineWidth',2);
            
            % Max computed secant grad
            hlp = Z(xi,yi) + sign(Z(xm,ym)-Z(xi,yi)) * estMax * dist;
            lineobj{end+1} = plot3([xp X(xm,ym)],[yp Y(xm,ym)], [Z(xi,yi) hlp],'r','LineWidth',2);
            
            % Estimate function secants roughly
            y = [xp; yp];
            inner = di < k.x0;
            theta = k.xR ./ di(inner);
            t1 = repmat(theta,2,1);
            yvec = repmat(y,1,size(centers(:,inner),2));
            centers(:,inner) = (1-t1).*centers(:,inner) + t1.*yvec;
            fxi = c * k.evaluate(centers,centers)';
            xisec = (fxi - k.evaluateScalar(norm(y))) ./  di;
            
            hlp = c .* lip;
            posi = max(hlp,0);
            [maxp, idp] = max(xisec);
            neg = 0;
            idm = 1;
            %neg = min(hlp,0);
            %[maxm, idm] = min(neg);
            
            %         posidx = hlp >=0;
            %         negidx = hlp < 0;
            %         posnum = sum(posidx);
            %         negnum = sum(negidx);
            %         e1 = zeros(posnum,negnum);
            %         e2 = zeros(posnum,negnum);
            %         for pidx = 1:max(1,posnum)
            %             for nidx = 1:max(1,negnum)
            %                 e1(pidx,nidx) = sum( (posi .* F(pidx,:) + abs(neg) .* F(nidx,:)));
            %                 e2(pidx,nidx) = sum( (posi .* W(pidx,:) + abs(neg) .* W(nidx,:)));
            %             end
            %         end
            %          if all(e1(:) < abs(effMax))
            %             e1/abs(effMax)
            %             disp(round(e1/abs(effMax)*100));
            %             fprintf('e1 all too low: ');
            %             fprintf('%2.2f%%, ',round(e1/abs(effMax)*100));
            %             fprintf('\n');
            %          end
            %         if all(e2(:) < abs(effMax))
            %disp(round(e2/abs(effMax)*100));
            %             fprintf('e2 all too low: ');
            %             fprintf('%2.2f%%, ',round(e2/abs(effMax)*100));
            %             fprintf('\n');
            %         end
            
            % Experiment 1
            %experi = maxp + abs(maxm);
            experi = sum( (posi .* F(idp,:) + abs(neg) .* F(idm,:)));
            
            % Plot Exp1
            hlp = Z(xi,yi) + sign(Z(xm,ym)-Z(xi,yi)) * experi * dist;
            lineobj{end+1} = plot3([xp X(xm,ym)],[yp Y(xm,ym)], [Z(xi,yi) hlp],'m','LineWidth',2);
            
            % Experiment 2
            %acoeff = 1-getAngle([xp; yp],xP,yP)/pi;
            experi2 = sum( (posi .* W(idp,:) + abs(neg) .* W(idm,:)));
            
            % Plot Exp2
            hlp = Z(xi,yi) + sign(Z(xm,ym)-Z(xi,yi)) * experi2 * dist;
            lineobj{end+1} = plot3([xp X(xm,ym)],[yp Y(xm,ym)], [Z(xi,yi) hlp],'c.-','LineWidth',2);
            
            
            
            %         title(sprintf(['Secant estimation: %f, effective secant: %f (ratio: %3d%%)\n'...
            %             'Global (estimated/effective): sharpest %d%%, coarsest %d%%'],estMax,abs(effMax),round(estMax/abs(effMax)*100),sharpest,coarsest));
            title(sprintf(['Secant estimation: %f, effective secant: %f (ratio: %3d%%)\n'...
                'Experimental 1: %f (ratio:%d%%), Experimental 2: %f (ratio:%d%%)'],estMax,abs(effMax),round(estMax/abs(effMax)*100),...
                experi,round(experi/abs(effMax)*100),...
                experi2,round(experi2/abs(effMax)*100)));
            
            subplot(1,2,2);
            hold off;
            % Max real secant
            Zs = abs(Zsecant);
            surf(X,Y,Zs)
            shading interp;
            hidden off;
            colormap hot;
            axis tight;
            hold on;
            h = contour3(X,Y,Zs,20,'black');
            [dummy, maxi] = max(Zs(:));
            plot3(X(maxi),Y(maxi),Zs(maxi),'black*','MarkerSize',15);
            
        catch ME
            printErrorStack(ME);
        end
    end

    function txt = pointSelected(hObject, eventdata)
        info = getCursorInfo(cm);
        pos = info.Position;
        if ~isequal(oldpos,pos(1:2))
            oldpos = pos(1:2);
            updatePlot(pos(1),pos(2));
        end
        txt = num2str(pos);
    end

    function [estMax, di, lip] = getEstGradient(xp,yp)
        % Estimations
        di = sqrt((xP-xp).^2 + (yP-yp).^2);
        lip = lipFun(di,Inf,0,[]);
        % estimated secant gradient
        estMax = abs(c) * lip';
    end

    function [estMax, di, lip] = getEstGradientJan(xp,yp,sameidx)
        % Estimations
        di = sqrt((xP-xp).^2 + (yP-yp).^2);
        liptest = lipFun(di,Inf,0,[]);
        di = di(sameidx);
        lip = lipFun(di,Inf,0,[]);
        % estimated secant gradient
        estMax = abs(c(sameidx)) * lip';
    end

    function [effMax, xm, ym, dist, sgrad] = getEffGradient(xp,yp)
        
        [dummy, yi] = min(abs(range - xp));
        [dummy, xi] = min(abs(range - yp));
        %yi = find(range == xp),1);
        %xi = find(range == yp),1);
        
        % Max real secant
        fdiff = Z - Z(xi,yi);
        xdiff = sqrt((X-xp).^2+(Y-yp).^2);
        sgrad = fdiff ./ (xdiff+eps);
        %sgrad(xi,yi) = k.evaluateD1(sqrt(xp^2+yp^2));
        [val, id] = max(abs(sgrad(:)));
        [xm,ym] = ind2sub(size(sgrad),id);
        
        % Distance
        dist = sqrt((X(xi,yi)-X(xm,ym)).^2 + (Y(xi,yi)-Y(xm,ym)).^2);
        % effective secant gradient
        effMax = (Z(xm,ym)-Z(xi,yi)) / dist;
    end

    function D = getDistanceMat(xP,yP)
        %[mx,my] = meshgrid(xP,yP);
        %D = sqrt((mx-mx').^2+(my-my').^2);
        vec = [xP;yP];
        nsq = sum(vec.^2,1);
        n = size(vec,2);
        D = sqrt((ones(n,1)*nsq)' + ones(n,1)*nsq - 2*(vec'*vec));
    end

    function A = getAngleMat(xP,yP)
        vec = [xP; yP];
        sp = vec' * vec;
        no = sqrt(sum(vec.^2));
        no = no'*no;
        A = acos(sp ./ no);
    end

    function A = getAngle(v,xP,yP)
        vec = [xP; yP];
        sp = vec' * v;
        no1 = sqrt(sum(vec.^2));
        no2 = sqrt(sum(v.^2));
        no = no1'*no2;
        A = acos(sp ./ no);
    end

end








