function [ time,weights,selectedIndex, recSigmaSquare,recoveredValuesSD] = multiTaskBCS_Hua_alpha0( phi, y, eta )
% multiple tasks

% phi : the measurement matrix (size: N*M)
% phi{i} = [phi(i,1),phi(i,2).....phi(i,M)], where phi(i,n) is N*1 column
% vector
% y : the noised observed measurement vector (size N*1)



% version: October, 2017
% Hua Bai

maxIteration = 10000;

testMode = false;

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
    sigmaSquare = sigmaSquare + var(y{i})*0.1/L;
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

while 1
    [maxML, selectedIndex] = max(ml_sum);
     mlPlot(1) = maxML;
    % equation (39) at multitask compressive sensing
    alpha = L./...
        sum((phiYProd(selectedIndex,:).^2./phiInnerPord(selectedIndex,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex,:));
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
%     for j  =  1 : phiColNum
%       S(j,i) = phiInnerPord(j,i)/sigmaSquare - ...
%           phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * phi{i}(:,j) / sigmaSquare;
%       Q(j,i) = phiBYProd(j) -...
%           phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * y{i} / sigmaSquare;
%     end
end

tic 
for i = 2 : maxIteration
    s = S; q = Q; 
    alphaRow = repmat(alpha,1,L);
    % update the s matrix, q matrix and g matrix
    % s( and q,g) matix, the row stands for the index of slected basis column,
    % the column stands for which snapshot
    s(selectedIndex,:) = alphaRow.*S(selectedIndex,:)./(alphaRow - S(selectedIndex,:));
    q(selectedIndex,:) = alphaRow.*Q(selectedIndex,:)./(alphaRow - S(selectedIndex,:));
    
    
    % search the next candidate basis according to the marginal likelihood
    twoDeltaML = -inf * ones(phiColNum, L);
    theta = sum((q.^2 - s)./(s.^2),2);
    thetaLg0Index = find(theta > 0);
    thetaSm0Index = find(theta < 0);
    newAlpha = L./theta;
    
    % another way to update alpha. 5/15/2018, little change compared with above
%     theta = sum((q.^2 - s),2);
%     thetaLg0Index = find(theta > 0);
%     thetaSm0Index = find(theta < 0);
%     newAlpha = sum(s.^2,2)./theta;
    % calculate the change of marginal likelihood for re-estimation
    [reestimateIndex,~,alphaIndex] = intersect(thetaLg0Index,selectedIndex);    
    if ~isempty(reestimateIndex)      
        newAlphaRow = repmat(newAlpha(reestimateIndex),1,L);
        oldAlphaRow = repmat(alpha(alphaIndex),1,L);
        deltaAlpha = 1./newAlphaRow - 1./oldAlphaRow;
        qReestimate = Q(reestimateIndex,:);
        sReestimate = S(reestimateIndex,:);
        twoDeltaML(reestimateIndex,:) = ...
            (qReestimate.^2)./(sReestimate + 1./deltaAlpha) - ...
            log(1 + sReestimate.*deltaAlpha);
    end
    % calculate the change of marginal likelihood for adding new basis
    addIndex = setdiff(thetaLg0Index,selectedIndex);
    if ~isempty(addIndex)
        qAdd = Q(addIndex,:);
        sAdd = S(addIndex,:);
        twoDeltaML(addIndex,:) = (qAdd.^2 - sAdd)./sAdd + log(sAdd./(qAdd.^2));
    end
    % calculate the change of marginal likelihood for deleting existing basis
    [deleteIndex,~,alphaIndex] = intersect(thetaSm0Index,selectedIndex); 
    if ~isempty(deleteIndex)
        oldAlphaRow = repmat(alpha(alphaIndex),1,L);
        qDel = Q(deleteIndex,:);
        sDel = S(deleteIndex,:);
        twoDeltaML(deleteIndex,:) = ...
            qDel.^2./(sDel - oldAlphaRow) - ...
            log(1 - sDel./oldAlphaRow);
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
