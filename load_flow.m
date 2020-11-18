function [ node_voltage,line_currents ] = load_flow( branches,root_node,base_voltage ,spot_loads,DG_loads )
%load_flow backward/forward loadflow of radial distribution systems
%   Inputs : branches,root_node,base_voltage ,spot_loads[,DG_loads]
%   Outputs: node_voltage,line_currents

assert(is_radial(branches),...
    'Please check your branches and make sure it represents a radial system.')
%% Init Data
if (nargin == 6) % combine new load with base load
    for i = 1:size(DG_loads,1)
        not_in_list = true;
        for j = 1:size(spot_loads,1)
            % check if there is already a load in the node the DG is added
            if (spot_loads(j,1) == DG_loads (i,1))
                spot_loads(j,2) = spot_loads(j,2) + DG_loads(i,2);
                not_in_list = false;
                break;
            end
        end
        % add new load to end of spot_loads if there is not any load in the DG node
        if (not_in_list)
            spot_loads(size(spot_loads,1)+1,:) = [DG_loads(i,1) DG_loads(i,2)];
        end
    end
end
branches = change_root_node(branches, root_node); % root node check
nodes = unique(branches(:,1:2)); % extract nodes
path = path_finder(branches,nodes);
Nx = size(nodes,1);
node_voltage = [nodes, base_voltage * ones(Nx,1)]; % Init volatage
node_voltage_prev = [nodes, zeros(Nx,1)]; % Could be any value other than node_voltage
tolerance = 1e-2; % termination criteria
%% FF LoadFlow
while any(any (abs(node_voltage(:,2) - node_voltage_prev(:,2)) > tolerance))
    % Set the last computed voltage to current value
    node_voltage_prev = node_voltage;
    load_currents = load_current(node_voltage,spot_loads);
    % Since this is a radial network, end node of each line is unique,
    % so a line can be represented by end node of it.
    line_currents = branch_current(path,load_currents);
    for i = 1:Nx % compute new voltages
        V1 = node_voltage(i,2); % strat node voltage
        V1n = node_voltage(i,1); % strat node voltage number
        for j = 1: Nx-1
            if (V1n == branches(j,1))
                V2n = branches(j,2); % the end node voltage number
                Z = branches(j,3); % the branch impedance
                for x = 1:Nx
                    if (V2n == line_currents(x,1))
                        I = line_currents(x,2); % the line current
                        V2 = V1 - Z * I; % end node voltage
                        for k = 1: Nx
                            if (node_voltage(k,1) == V2n)
                                node_voltage(k,2) = V2; % update computed voltage
                                break;
                            end
                        end
                        break;
                    end
                end
            end
        end
    end
end
end

%% sub functions:
function y = is_radial(branches)
%is_radial checks whether a grid defined by Branch is radial
%   The idea is checking if there exist two lines ending in one node.
%   if there are more than one line ending in one node the the unique
%   function will get rid of it! and number of unique ending are less than
%   all of them.
y=length(unique(branches(:,2))) == length(branches(:,2));
end

function [ Branch ] = change_root_node( Branch,root_node )
%change_root_node changes the direction of branches to make the input node
% be the root of the tree
Bx = size (Branch,1);
try
    for i = 1:Bx
        if (Branch(i,2) == root_node)
            NewRoot = Branch(i,1);
            Branch(i,2) = NewRoot;
            Branch(i,1) = root_node;
            Branch = change_root_node_recursive(Branch,NewRoot,root_node);
        elseif (Branch(i,1) == root_node)
            NewRoot = Branch(i,2);
            Branch = change_root_node_recursive(Branch,NewRoot,root_node);
        end
    end
catch ME
    if strcmp('MATLAB:lang:StackOverflow',ME.identifier)
        msg = ['It is likely provided configuration is not radial.', ...
            '\nPlease check your branches and', ...
            ' make sure it represents a radial system.\n',...'
            'If the error persists contact me at', ...
            ' amin.gholizad@gmail.com.'];
        causeException = MException('MATLAB:Amin:StackOverflow',msg);
        ME = addCause(ME,causeException);
    end
    rethrow(ME)
end
end
function [ Branch ] = change_root_node_recursive( Branch,root_node,PreviuseRoot )
%change_root_node_recursive Internal part of change_root_node function
Bx = size (Branch,1);
for i = 1:Bx
    if (Branch(i,2) == root_node) && (Branch(i,1) ~= PreviuseRoot)
        NewRoot = Branch(i,1);
        Branch(i,2) = NewRoot;
        Branch(i,1) = root_node;
        Branch = change_root_node_recursive(Branch,NewRoot,root_node);
    elseif (Branch(i,1) == root_node)
        NewRoot = Branch(i,2);
        Branch = change_root_node_recursive(Branch,NewRoot,root_node);
    end
end
end

function [ path ] = path_finder( Branch,nodes )
%path_finder will find all sub nodes of a node
Bx = size(Branch,1);
path = zeros(Bx+2);
path(:,1) = [0; nodes];
path(1,:) = [0, nodes.'];
for i = 2:Bx+2
    path(i,i) = 1;
    A = path(i,1);
    j = 1;
    while (j <= Bx)
        if (A == Branch(j,2))
            B = Branch(j,1);
            for k = 2:Bx+2
                if(path(1,k) == B)
                    path(k,i) = 1;
                    A = B;
                    j = 0;
                    break;
                end
            end
        end
        j = j + 1;
    end
end
end

function [ load_current ] = load_current( NodeVoltage,SpotLoadNode )
%load_current computes loads current
Current = 0;
Lx = size(SpotLoadNode,1);
load_current = zeros(Lx,2);
for i = 1:Lx
    for j = 1:size(NodeVoltage,1)
        if (NodeVoltage(j,1) == SpotLoadNode(i,1))
            Current = conj(SpotLoadNode(i,2)./NodeVoltage(j,2));
            break;
        end
    end
    load_current(i,:) = [SpotLoadNode(i,1) Current];
end
end

function [ branch_current ] = branch_current( path,load_current )
%branch_current copmutes Current of each line
Lx = size (load_current,1);
Px = size(path,1);
Current = 0;
branch_current = zeros(Px-1,2); % create the output
for i = 2:Px
    for j = 2:Px
        if path(i,j)
            for k = 1:Lx
                if path(1,j) == load_current(k,1)
                    Current = Current + load_current(k,2);
                    break;
                end
            end
        end
    end
    branch_current(i-1,:) = [path(i,1) Current];
    Current = 0;
end
end