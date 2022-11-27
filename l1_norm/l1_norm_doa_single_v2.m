function [time,weights,selectedIndex, recSigmaSquare,recoveredValuesSD] = ...
    l1_norm_doa_single_v2 (phi, y, eta)

% phi : the measurement matrix (size: N*M)
% phi{i} = [phi(i,1),phi(i,2).....phi(i,M)], where phi(i,n) is N*1 column
% vector
% y : the noised observed measurement vector (size N*1)



% version: May, 2018
% Hua Bai

maxIteration = 1000;
testMode = true;
lambda = 0;
if iscell(y)
    L = length(y); % number of snapshots
else
    L = 1;
%     y = {y};
%     phi = {phi};
end

for i = 1 : L
    [Row(i), Col(i)] = size(phi{i});
end
if sum(abs(diff(Col))) ~= 0 
    error('The basis matrix should have same columns.');
else
    phiColNum = Col(1);
end

% initial guess for the noise
sigmaSquare = 0;
sigmaSquareCpp = 0;
for i = 1 : L
    sigmaSquare = sigmaSquare + var(y{i})*0.01/L;
    sigmaSquareCpp = sigmaSquareCpp + (y{i}- mean(y{i}))'*(y{i}- mean(y{i}))/(length(y{i}-1))*0.1/L;
end

% calculate the phiInnerProd and phiYprod
for i = 1 : L
   phiInnerPord(:,i) = sum(phi{i}.*phi{i})'; 
   phiYProd(:,i) = phi{i}'*y{i};
end

% only for initial likelihood 
ml = phiYProd.^2./phiInnerPord; 
ml_sum = sum(ml,2);

% choose the first basis
while 1
    [maxML, selectedIndex] = max(ml_sum);
     mlPlot(1) = maxML;
    % equation (39) at multitask compressive sensing
    alpha = L./...
        sum((phiYProd(selectedIndex,:).^2./phiInnerPord(selectedIndex,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex,:));
    gamma = 1/alpha;
    if alpha > 0
        break;
    else
        ml_sum(selectedIndex) = 0;
    end
end


% calculate S,Q and G matrix
% only the selected phi are included 
for i = 1 : L
    phiSelected{i} = phi{i}(:,selectedIndex);
    covMatrix{i} = inv(phiSelected{i}' * phiSelected{i} / sigmaSquare + alpha);
    mu{i} = covMatrix{i} * phiSelected{i}' * y{i} / sigmaSquare;
    phiBPhiProd = phi{i}' * phiSelected{i} / sigmaSquare;
    phiBYProd = phi{i}' * y{i} / sigmaSquare;
    S(:,i) = phiInnerPord(:,i)/sigmaSquare - covMatrix{i}*(phiBPhiProd).^2;
    Q(:,i) = phiBYProd - covMatrix{i}*phiBPhiProd*phiSelected{i}'*y{i}/sigmaSquare;
end
deleted = [];
tic
for i = 2 : maxIteration
    s = S; q = Q; 
    alphaRow = repmat(alpha,1,L);
    % update the s matrix, q matrix and g matrix
    % s( and q,g) matix, the row stands for the index of slected basis column,
    % the column stands for which snapshot
    s(selectedIndex,:) = alphaRow.*S(selectedIndex,:)./(alphaRow - S(selectedIndex,:));
    q(selectedIndex,:) = alphaRow.*Q(selectedIndex,:)./(alphaRow - S(selectedIndex,:));
     
    % update lambda
    lambda = (length(selectedIndex) - 1) / sum(1/(2*alpha));
    % search the next candidate basis according to the marginal likelihood
    twoDeltaML = -inf * ones(phiColNum, L);
    theta = sum((q.^2 - s),2);
    thetaLg0Index = find(theta > lambda);
    thetaSm0Index = find(theta < lambda);
    % single snapshot
%     newGamma = sum((q.^2 - s)./(s.^2),2)/L;
    A = lambda + s - q.^2;
    B = 2*lambda.*s + s.^2;
    C = lambda.*s.^2;
    discriminant = B.^2 - 4.*A.*C;
    
    newAlpha = (-B - sqrt(discriminant) ) ./ (2*A);
    newGamma = 1./newAlpha;
    % the reestimate index: theta>0 and the basis is selected (in the selected index) 
    [reestimateIndex,~,alphaIndex] = intersect(thetaLg0Index,selectedIndex);
    % the change of likelihood with reestimate
    if ~isempty(reestimateIndex)
        newAlphaRow = repmat(newAlpha(reestimateIndex),1,L);
        oldAlphaRow = repmat(alpha(alphaIndex),1,L);
        deltaAlpha = 1./newAlphaRow - 1./oldAlphaRow;
        qReestimate = q(reestimateIndex,:);
        sReestimate = s(reestimateIndex,:);
        twoDeltaML(reestimateIndex,:) = ...
            (qReestimate.^2)./(newAlphaRow + sReestimate) + log(newAlphaRow ./ (newAlphaRow + sReestimate)) -...
            (qReestimate.^2)./(oldAlphaRow + sReestimate) - log(oldAlphaRow ./ (oldAlphaRow + sReestimate))...
            -  lambda./ newAlphaRow +  lambda./ oldAlphaRow;
    end
    
    % the add index: theta>0 and the basis is NOT selected ( NOT in the selected index) 
    addIndex = setdiff(thetaLg0Index,selectedIndex);
    if ~isempty(addIndex)
       qAdd = q(addIndex,:);
       sAdd = s(addIndex,:);
       twoDeltaML(addIndex,:) = (qAdd.^2-sAdd)./sAdd...
                              +log(sAdd./(qAdd.^2)) - lambda./newAlpha(addIndex);
       which = intersect(deleted,addIndex);
       twoDeltaML(which) = -inf;
    end
    
    % the delete index: theta<=0 and the basis is selected ( in the selected index) 
    [deleteIndex,~,alphaIndex] = intersect(thetaSm0Index,selectedIndex);
    if ~isempty(deleteIndex)
        if length(selectedIndex) == 1
            twoDeltaML(deleteIndex,:) = -inf;
        else
            qDel = q(deleteIndex,:);
            sDel = s(deleteIndex,:);
            twoDeltaML(deleteIndex,:) = -qDel.^2./(sDel + alpha(alphaIndex))...
                -log(alpha(alphaIndex)./(alpha(alphaIndex)+sDel)) + lambda./alpha(alphaIndex);
        end
    end
    
    [maxML(i),maxDeltaMLIdx] = max(sum(real(twoDeltaML),2));
    mlPlot(i) = mlPlot(i-1) + maxML(i);
    if (i > 2) && (abs(maxML(i) - maxML(i-1)) < abs(maxML(i) - maxML(1)) * eta)
        break;
    end
    
    % update the alpha, mu and covariance and (S, Q)
    alphaIndex = find(selectedIndex == maxDeltaMLIdx);
    if ismember(maxDeltaMLIdx, reestimateIndex)
        if testMode
            disp(['Iteration ', num2str(i), ': Reestimate']);
        end
        deltaAlpha = newAlpha(maxDeltaMLIdx) - alpha(alphaIndex);
        for snapShot = 1 : L           
            covII = covMatrix{snapShot}(alphaIndex,alphaIndex);
            covCol = covMatrix{snapShot}(:,alphaIndex);
            muOld = mu{snapShot}(alphaIndex);
            kj = deltaAlpha/(deltaAlpha * covII + 1);
            mu{snapShot} = mu{snapShot} - kj*muOld*covCol;
            covMatrix{snapShot} = covMatrix{snapShot} - kj*(covCol*covCol');
            comm = phi{snapShot}'*(phiSelected{snapShot}*covCol)/sigmaSquare;
            S(:,snapShot) = S(:,snapShot) + kj*comm.^2;
            Q(:,snapShot) = Q(:,snapShot) + kj*muOld*comm;
        end
        % update alpha
        alpha(alphaIndex) = newAlpha(maxDeltaMLIdx);
        % add basis column
    elseif ismember(maxDeltaMLIdx, addIndex)
        if testMode
            disp(['Iteration ', num2str(i), ': Add']);
        end
        for snapShot = 1 : L
            covII = 1/(newAlpha(maxDeltaMLIdx) + S(maxDeltaMLIdx,snapShot));
            muAdd = covII*Q(maxDeltaMLIdx,snapShot);
            phiI = phi{snapShot}(:,maxDeltaMLIdx);
            
            comm1 = covMatrix{snapShot}*(phiSelected{snapShot}'*phiI)/sigmaSquare;
            ei = phiI - phiSelected{snapShot}*comm1;
            mu{snapShot} = [mu{snapShot} - muAdd*comm1; ...
                muAdd];
            comm2 = -covII*comm1;
            covMatrix{snapShot} = [covMatrix{snapShot} + covII*comm1*comm1', comm2; ...
                comm2', covII];
            comm3 = phi{snapShot}'*ei/sigmaSquare;
            S(:,snapShot) = S(:,snapShot) - covII * comm3.^2;
            Q(:,snapShot) = Q(:,snapShot) - muAdd * comm3;
            phiSelected{snapShot} = [phiSelected{snapShot},phiI];
        end
        selectedIndex = [selectedIndex;maxDeltaMLIdx];
        alpha = [alpha;newAlpha(maxDeltaMLIdx)];
        % delete existing basis
    elseif ismember(maxDeltaMLIdx,deleteIndex)
        if testMode
            disp(['Iteration ', num2str(i), ': Delete']);
        end
        for snapShot = 1 : L
            covII = covMatrix{snapShot}(alphaIndex,alphaIndex);
            covCol = covMatrix{snapShot}(:,alphaIndex);
            muOld = mu{snapShot}(alphaIndex);
            covMatrix{snapShot} = covMatrix{snapShot} - 1/covII*(covCol*covCol');
            covMatrix{snapShot}(:,alphaIndex) =[];
            covMatrix{snapShot}(alphaIndex,:) =[];
            mu{snapShot} = mu{snapShot} - muOld/covII*covCol;
            mu{snapShot}(alphaIndex) = [];
            comm = phi{snapShot}'*(phiSelected{snapShot}*covCol)/sigmaSquare;
            S(:,snapShot) = S(:,snapShot) + comm.^2/covII;
            Q(:,snapShot) = Q(:,snapShot) + muOld/covII*comm;            
            phiSelected{snapShot}(:,alphaIndex) = [];
        end
        deleted = [deleted maxDeltaMLIdx];
        selectedIndex(alphaIndex) = [];
        alpha(alphaIndex) = [];
    end
    
end

weights = zeros(phiColNum,L);
recoveredValuesSD = zeros(length(selectedIndex),1);
for i = 1 : L
    weights(selectedIndex,i) = mu{i}; 
    recoveredValuesSD = recoveredValuesSD + sqrt(diag(covMatrix{i}))/L;
end

% to estimate the noise variance 
numeratorSum = 0;
denominatorSum = 0;
for i = 1 : L
    numeratorSum = numeratorSum + sum((y{i} - phiSelected{i}*mu{i}).^2);
    denominatorSum = denominatorSum + (Row(i) - length(selectedIndex) + alpha'*diag(covMatrix{i}));
end
recSigmaSquare = numeratorSum/denominatorSum;
% recSigmaSquare = sigmaSquare;
time = toc; 

end

