function [CLUSTER_PATH,SAVE_PATH]=initialize(varargin)
    % Initializes paths.
        
    p = inputParser();
    addRequired(p,'saveName');
    parse(p, varargin{:});
    ops = p.Results;
    
    % paths
    SAVE_PATH = [pwd filesep 'output' filesep ops.saveName filesep];
    CLUSTER_PATH = [pwd filesep];
    addpath(genpath([CLUSTER_PATH 'utils']));
        
end