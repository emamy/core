function EpsLoss(ep)
% Plots the epsilon-insensitive loss function for a given epsilon.

if nargin == 0
    ep = .1;
end

x = -5*ep:ep/10:5*ep;
lx = max(0, abs(x)-ep);
h = figure;
plot(x,zeros(size(x)),'k--',x,lx,'b',[0 0+eps],[-ep/3 max(lx)],'black--');
annotation(h,'textbox',[0.4258 0.113 0.0461 0.07664],...
    'String',{'-\epsilon'},...
    'FontWeight','bold',...
    'FontSize',20,'LineStyle','none');

% Create textbox
annotation(h,'textbox',[0.5858 0.1229 0.0394 0.06729],...
    'String',{'\epsilon'},...
    'FontWeight','bold',...
    'FontSize',20,'LineStyle','none');
axis([x(1) x(end) -ep/3 max(lx)]);
end