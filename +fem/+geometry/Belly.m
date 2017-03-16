function [geo27, rx, ry] = Belly(xgrid, varargin)
% Code to generate a belly-shaped geometry with various options.
if ~isnumeric(xgrid)
    error('At least an xgrid numeric value has to be specified.');
end
i = inputParser;
i.KeepUnmatched = true;
i.addParamValue('Radius',1);
i.addParamValue('InnerRadius',.5);
i.addParamValue('Gamma',2);
i.addParamValue('Center',[]);
i.addParamValue('Layers',1);
i.addParamValue('MinZ',[]);
i.addParamValue('Half',false);
i.parse(varargin{:});
opt = i.Results;
if isscalar(xgrid)
    xgrid = linspace(0,xgrid,4);
end
% Make sure its a row vector
xgrid = reshape(xgrid,1,[]);
% refine the points as we have quadratic elements
xgrid = [xgrid; xgrid(1:end-1)+diff(xgrid)/2 0];
xgrid = xgrid(1:end-1);

if isa(opt.Radius,'function_handle')
    fx = opt.Radius(xgrid);
    if size(fx,1) == 1
        fx = [fx; fx];
    end
else
    if isscalar(opt.Gamma)
        opt.Gamma = [opt.Gamma opt.Gamma];
    end
    if isscalar(opt.Radius)
        opt.Radius = [opt.Radius opt.Radius];
    end
    k = kernels.GaussKernel;
    k.Gamma = opt.Gamma(1);
    if ~isempty(opt.Center)
        cent = opt.Center;
    else
        cent = (xgrid(:,end)-xgrid(:,1))/2;
    end
    fx(1,:) = k.evaluate(xgrid,cent)'*opt.Radius(1);
    k.Gamma = opt.Gamma(2);
    fx(2,:) = k.evaluate(xgrid,cent)'*opt.Radius(2);
end
if isscalar(opt.InnerRadius)
    opt.InnerRadius = [1 1]*opt.InnerRadius;
end

posfun = @(omega,offset)[sin(omega+offset); cos(omega+offset)];

qpart = getQuarterPart(opt.Layers);
baseelems = zeros(size(qpart,1),size(qpart,2),4);
baseelems(:,:,1) = qpart;
nlayerelems = size(baseelems,1)/2;
for k=1:nlayerelems
    pos = (1:2)+2*(k-1);
    baseelems(pos,:,2) = [-1 0; 0 1]*baseelems(pos,[3 2 1 6 5 4 9 8 7],1);
end
if ~opt.Half
    baseelems(:,:,3) = -baseelems(:,9:-1:1,1);
    baseelems(:,:,4) = -baseelems(:,9:-1:1,2);
end

nelems = floor(length(xgrid)/2);
npp = 27*nelems*nlayerelems;
nodes = zeros(3,npp*4);
rx = opt.InnerRadius(1) + fx(1,:);
ry = opt.InnerRadius(2) + fx(2,:);
for p = 1:nelems
    elempos = (p-1)*2 + (1:3);
    rx_ = ones(9,1) * rx(elempos);
    ry_ = ones(9,1) * ry(elempos);
    z = ones(9,1) * xgrid(elempos);
    for k=1:nlayerelems
        elem = baseelems((1:2)+2*(k-1),:,:);
        for rotpart = 1:4
            nodepos = (((p-1) * nlayerelems + (k-1))*4 + (rotpart-1))*27 + (1:27);
            nodes(1,nodepos) = bsxfun(@times,repmat(elem(1,:,rotpart),1,3),rx_(:)');
            nodes(2,nodepos) = bsxfun(@times,repmat(elem(2,:,rotpart),1,3),ry_(:)');
            nodes(3,nodepos) = z(:);
        end
    end
end

% Filter nodes to unique ones
[nodes, ~, elems] = unique(nodes','rows','stable');

% For the "not so unique but same" nodes use pdist
D = pdist(nodes);
[D, idx] = sort(D,'ascend');
same = D < eps*max(D);
[i,j] = find(tril(ones(size(nodes,1)),-1));
i = i(idx(same)); j = j(idx(same));
% Kick out double nodes
nodes(j,:) = [];
% Correct element indices after node removal
% - use i indices (instead of to-remove-j)
for k=1:length(i)
    elems(elems == j(k)) = i(k);
end
% - Get effective number of nodes
effidx = unique(elems(:));
% - Set mapping indices at old positions
invidx(effidx) = 1:length(effidx);
% - Use those!
elems = invidx(elems);

% Get the right shape
nodes = nodes';
elems = reshape(elems,27,[])';

% This little hack allows to restrict the minimum z node
% position to a certain value, so as if the muscle belly is
% "lying on a table"
%
% This only works on the "downside" faces, so if the threshold
% is larger than the intermediate node position of the side
% faces, the muscle geometry will appear "dented".
if ~isempty(opt.MinZ)
%     error('Fixme');
    hlp_g27 = fem.geometry.Cube27Node;
    centerline_facenodeidx = hlp_g27.MasterFaces(3,:);
    for elem = [3:4:4*nelems 4:4:4*nelems]
        centerline_nodes = elems(elem,centerline_facenodeidx);
        toScale = nodes(2,centerline_nodes) < opt.MinZ;
        factor = opt.MinZ./nodes(2,centerline_nodes(toScale));
        nodes(2,centerline_nodes(toScale)) = opt.MinZ;
        % Also scale inner node positions
        inner_nodes = elems(elem,[4:6 13:15 22:24]);
        nodes(2,inner_nodes(toScale)) = nodes(2,inner_nodes(toScale)).*factor;
    end
end

geo27 = fem.geometry.Cube27Node(nodes,elems);
geo27.swapYZ;

    function elems = getQuarterPart(elem_layers)
        if nargin < 1
            elem_layers = 9;
        end
        if isscalar(elem_layers)
            elem_layers = linspace(0,1,elem_layers+1);
            elem_layers = elem_layers(2:end);
        end
        elem_layers = elem_layers / max(elem_layers);
        num_elem_layers = length(elem_layers);
        
        % Augment layer count if odd
        if mod(num_elem_layers,2) == 1
            elem_layers(end+1) = elem_layers(end)*1.1;
        end
        
        % Create innermost element
        o = posfun(linspace(0,pi/2,5),0);
        c = posfun(pi/4,0)/2;
        elems(1:2,:) = [0 .5 1 0 c(1) o(1,4) o(1,1:3)
            0 0  0 .5 c(2) o(2,4) o(2,1:3)]*elem_layers(1);
        
        if num_elem_layers > 1
            o = elem_layers(2)*posfun(linspace(0,pi/2,5),0);
            elems(3:4,:) = [elems(1:2,7:9) (o(:,1:3)+elems(1:2,7:9))/2 o(:,1:3)];
            elems(5:6,:) = [elems(1:2,[9 6 3]) (o(:,3:5)+elems(1:2,[9 6 3]))/2 o(:,3:5)];
            % Reorder
            elems(5:6,:) = elems(5:6,[3 6 9 2 5 8 1 4 7]);
        end
        
        % number of blocks needed (two elem layers per block)
        num_block_layers = ceil(num_elem_layers/2)-1;
        
        omega = pi/2;
        blocks_per_loop = 1;
        for bn=1:num_block_layers
            layerpos = (2:4)+2*(bn-1);
            for b = 1:blocks_per_loop
                off = (b-1)*omega;
                new = getblock(elem_layers(layerpos),omega,off);
                if bn == num_block_layers && mod(num_elem_layers,2) ~= 0
                    % Cut out the outer layer if not wanted
                    new = new(1:6,:);
                end
                elems(end+1:end+size(new,1),:) = new;
            end
            omega = omega / 2;
            blocks_per_loop = blocks_per_loop*2;
        end
    end

    function elems = getblock(radii,omega,offset)
        i = radii(1)*posfun(linspace(0,omega,5),offset);
        m = radii(2)*posfun(linspace(0,omega,13),offset);
        m = m(:,[1 3 5:9 11 13]);
        o = radii(3)*posfun(linspace(0,omega,9),offset);
        
        % E1
        elems(1:2,:) = [i(:,1:3) (m(:,1:3)+i(:,1:3))/2 m(:,1:3)];
        % E3
        elems(5:6,:) = [i(:,3:5) (m(:,7:9)+i(:,3:5))/2 m(:,7:9)];
        center = (radii(1)+radii(2))/2*posfun(omega/2,offset);
        % E2 (split part)
        elems(3:4,:) = [elems(5:6,[1 4 7]) elems(1:2,6) center m(:,[6 3:5])];
        
        % E4-E7 top layer
        elems(7:8,:) = [elems(1:2,7:9) (o(:,1:3)+elems(1:2,7:9))/2 o(:,1:3)];
        elems(9:10,:) = [elems(3:4,7:9) (o(:,3:5)+elems(3:4,7:9))/2 o(:,3:5)];
        elems(11:12,:) = [elems(3:4,[9 6 3]) (o(:,5:7)+elems(3:4,[9 6 3]))/2 o(:,5:7)];
        elems(13:14,:) = [elems(5:6,7:9) (o(:,7:9)+elems(5:6,7:9))/2 o(:,7:9)];
        
        % Reorder for x-axis increase first
        %elems(5:6,:) = elems(5:6,[3 6 9 2 5 8 1 4 7]);
        %elems(11:12,:) = elems(11:12,[3 6 9 2 5 8 1 4 7]);
        %elems(13:14,:) = elems(13:14,[3 6 9 2 5 8 1 4 7]);
    end
end

