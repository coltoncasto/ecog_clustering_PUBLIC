function [IDX_new,C_new,SUMD_new,D_new] = reorderClusters(X,C_template,IDX,C_current,SUMD,D,dist)

    k = size(C_current,1);
    if strcmp(dist,'sqeuclidean')
        dist = 'squaredeuclidean';
    end
    dist = 'correlation'; % template matching

    % --- MAKE PAIRINGS --- %

    % calculate mean of cluster (in case using kmedoids)
    C_current_means = zeros(k,size(X,2));
    for kk=1:k
        C_current_means(kk,:) = mean(X(IDX==kk,:),1);
    end

    % calculate distances
    D_between = pdist2(C_template,C_current_means,dist);
    
    % algorithm for finding optimal pairs
    assign_new = zeros(1,k);
    for kk=1:k
        % find minimum of distance matrix (rows = template, cols = current)
        [col_mins, new_idxs] = min(D_between);
        [~, old_idx] = min(col_mins);

        % make assignment 
        assign_new(old_idx) = new_idxs(old_idx);

        % set row and column of distance matrix to inf
        D_between(:,old_idx) = inf;
        D_between(new_idxs(old_idx),:) = inf;

    end

    % --- REORDER VALUES ---

    % relabel clusters in manually specified order
    IDX_new = IDX;
    % C_new = zeros(size(C_current));
    % SUMD_new = zeros(size(SUMD));
    % D_new = zeros(size(D));
    for kk=1:k
        eval(strcat("idx",num2str(kk),"=find(IDX_new==kk);")); 
    end
    for kk=1:k
        eval(strcat("IDX_new(idx",num2str(kk),")=",num2str(assign_new(kk)),";")); 
        % C_new(kk,:) = C_current(assign_new(kk),:);
        % SUMD_new(kk) = SUMD(assign_new(kk));
        % D_new(:,kk) = D(:,assign_new(kk));
    end

    % relabel other things
    C_new = C_current(assign_new,:);
    SUMD_new = SUMD(assign_new);
    D_new = D(:,assign_new);
    
 
end