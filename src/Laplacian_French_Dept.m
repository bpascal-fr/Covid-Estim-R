% Create the graph discrete gradient for French departments, in which 
% departments sharing a terrestrial border are considered connected.
%
% From
% - Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P.,
% Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P.,
% & Garnier, N. (2020). Spatial and temporal regularization to estimate
% COVID-19 reproduction number R(t): Promoting piecewise smoothness via
% convex optimization. PlosOne, 15(8), e0237901


function G_User = Laplacian_French_Dept(Departments)

    % Inputs:  - Departments: structure containing information about the selected departments
    %                 Names: list of French departments names
    %                 Number: official INSEE numerotation of French departments
    %                 Indices: index of each department in the storage from 1 to D
    %
    % Ouputs:  - G_User: graph discrete gradient operator encoding the spatial connectivity structure

    % Load the information of French departments sharing a terrestrial border
    load('data/French_Contiguous_Departments.mat','terrestrial')

    

    % Select all departments but Corsica
    matrice_no_Corsica = terrestrial;   % Corsica, departments 2A and 2B, corresponds to on lines 29 and 30.
    num_depts          = [1:28,31:96]; % INSEE Numbers of French metropolitan departments:
    % - from 1 to 19, the column index is the department INSEE number,
    % - from 20 to 94, the column index is the department INSEE number plus one.
    terrestrial         = matrice_no_Corsica(num_depts,num_depts); % Connectivity between the 94 continental departments.


    % Initialize the graph discrete gradient operator
    G = zeros(length(find(terrestrial==-1)),size(terrestrial,2));

    % Fill it with -1 and 1 to implement pairwise differences
    kk = 1;
    for ll = 1:size(terrestrial,1)
        ind2 = find(terrestrial(ll,:)==-1);
        for jj = 1:length(ind2)
            G(kk,terrestrial(ll,:)==1) = 1;
            G(kk,ind2(jj)) = -1;
            kk = kk+1;
        end
    end

    % Patch for Corsica
    [E,C]          = size(G) ;
    Gtmp           = G ;
    G              = zeros(E+2,C+2) ;
    G(1:E,1:C)     = Gtmp ;
    G(E+1,C+1)     = 1 ;
    G(E+1,C+2)     = -1 ;
    G(E+2,C+1)     = -1 ;
    G(E+2,C+2)     = 1 ;

    % Add some zeros for French overseas territories
    G_All          = zeros(E+2,104);
    G_All(:,1:C+2) = G;

    if nargin == 0
        
        % All French departments including overseas territories
        G_User = G_All;

    else

        % Selected French departments with possible French overseas territories
        G_User = G_All(:,Departments.Indices);
        
    end

    % Remove the edges which are one-sided
    valid_edges = [];
    for edge    = 1:size(G_User,1)
        if sum(abs(G_User(edge,:)) > 0) < 2
            G_User(edge,:) = 0;
        else
            valid_edges    = [valid_edges, edge];
        end
    end
    if isempty(valid_edges)
        G_User             = 0; % no edge at all between the territories
        disp('here')
    else
        G_User             = G_User(valid_edges,:); % store only valid edges
    end
end
