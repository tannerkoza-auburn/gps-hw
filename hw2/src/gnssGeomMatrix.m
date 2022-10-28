function  G = gnssGeomMatrix(UnitVecs)
% DESCRIPTION: This function produces the geometry/observation matrix 
% used in the least squares estimation of global position using
% satellite pseudoranges.
% PARAMS:
%       UnitVecs: mx3 matrix of unit vectors to satellite(s) positions
% OUTPUT:
%       G: GNSS geometry matrix
% AUTHOR: Tanner Koza, M.E. (Master of Engineering) Candidate

%% Initialization

    numMeas = length(UnitVecs);
    
    uhat_x = UnitVecs(:,1);
    uhat_y = UnitVecs(:,2);
    uhat_z = UnitVecs(:,3);

%% Geometry Matrix Population

    G = [-uhat_x -uhat_y -uhat_z ones(numMeas,1) zeros(numMeas, 4);
         zeros(numMeas,4) -uhat_x -uhat_y -uhat_z ones(numMeas,1)];

end