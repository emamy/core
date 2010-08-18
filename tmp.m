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







