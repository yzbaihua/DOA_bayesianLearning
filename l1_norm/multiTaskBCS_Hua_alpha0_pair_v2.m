function [ time, weights,selectedIndex, selectedIdxAssociated, recSigmaSquare,recoveredValuesSD ]...
    = multiTaskBCS_Hua_alpha0_pair_v2( phi, y, eta )
% multiple tasks

% phi : the measurement matrix (size: N*M)
% phi{i} = [phi(i,1),phi(i,2).....phi(i,M)], where phi(i,n) is N*1 column
% vector
% y : the noised observed measurement vector (size N*1)

% version: October, 2017
% Hua Bai

maxIteration = 20000;

testMode = true;

if iscell(y)
    L = length(y); % number of snapshots
else
    L = 1;
    y = {y};
    phi = {phi};
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
for i = 1 : L
    sigmaSquare = sigmaSquare + var(y{i})*0.1/L;
end

% calculate the phiInnerProd and phiYprod
for i = 1 : L
   phiInnerPord(:,i) = sum(phi{i}.*phi{i})'; 
   phiYProd(:,i) = phi{i}'*y{i};
end

% only for initial likelihood 
ml = phiYProd.^2./phiInnerPord; 
% Oct 2017 version
% ml_sum = sum(ml,2);
% update 11/14/2017 likelihood is considered with pair
ml_sum = sum(ml(1:phiColNum/2,:) + ml(phiColNum/2+1:end,:),2);

while 1
    [maxML, selectedIndex] = max(ml_sum);
    
    % equation (39) at multitask compressive sensing
    % Oct 2017 version
%     alpha = L./...
%         sum((phiYProd(selectedIndex,:).^2./phiInnerPord(selectedIndex,:) - sigmaSquare)./...
%         phiInnerPord(selectedIndex,:));
    % update 11/14/2017 likelihood is considered with pair
    alphaReal = L./...
        sum((phiYProd(selectedIndex,:).^2./phiInnerPord(selectedIndex,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex,:));
    alphaImag = L./...
        sum((phiYProd(selectedIndex + phiColNum/2,:).^2./phiInnerPord(selectedIndex + phiColNum/2,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex + phiColNum/2,:));
    if (alphaReal > 0) && (alphaImag > 0)
        % update, we consider the pair
%         if selectedIndex <= phiColNum/2
%             selectedIdxAssociated = selectedIndex + phiColNum/2;
%         else
%             selectedIdxAssociated = selectedIndex - phiColNum/2;
%         end
        % update 11/14/2017 likelihood is considered with pair
        selectedIndex = [selectedIndex; selectedIndex + phiColNum/2];
        alpha = diag([alphaReal,alphaImag]);
        break;
    else
        ml_sum(selectedIndex) = 0;
    end
end
mlPlot(1) = maxML;

% calculate S,Q and G matrix
% only the selected phi are included 
for i = 1 : L
    phiSelected{i} = phi{i}(:,selectedIndex);
    covMatrix{i} = inv(phiSelected{i}' * phiSelected{i} / sigmaSquare + alpha);
    mu{i} = covMatrix{i} * phiSelected{i}' * y{i} / sigmaSquare;
    phiBPhiProd = phi{i}' * phiSelected{i} / sigmaSquare;
    phiBYProd = phi{i}' * y{i} / sigmaSquare;
    % Oct 2017 version
%     S(:,i) = phiInnerPord(:,i)/sigmaSquare - covMatrix{i}*(phiBPhiProd).^2;
%     Q(:,i) = phiBYProd - covMatrix{i}*phiBPhiProd*phiSelected{i}'*y{i}/sigmaSquare;
    % update 11/14/2017 
    for j  =  1 : phiColNum
      S(j,i) = phiInnerPord(j,i)/sigmaSquare - ...
          phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * phi{i}(:,j) / sigmaSquare;
      Q(j,i) = phiBYProd(j) -...
          phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * y{i} / sigmaSquare;
    end
end
% update 11/14/2017 
alpha = diag(alpha);
selectedIdxAssociated = [];
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
    
    % update with pair 
    % record the status of each basis (it should be added, deleted or reestimated)
    % 1 : add or reestimate
    % 0 : delete and the rest
%     statusIdx = zeros(phiColNum,1);
%     statusIdx(addIndex) = 1;
%     statusIdx(reestimateIndex) = 1;
%     pairStatus = statusIdx(1:phiColNum/2) + ...
%                  statusIdx(phiColNum/2 + 1 : end);
%     sameStatusIdx = find(mod(pairStatus,2) == 0);          
%     diffStatusIdx = setdiff((1:phiColNum/2)',sameStatusIdx);
%     twoDeltaMLPair = zeros(phiColNum/2,1);
    % if the pair are both going to be added or reestimated, then the ML is
    % the sum of two individual ML

%     twoDeltaMLPairSum = zeros(phiColNum/2,1);
%     twoDeltaMLPair(sameStatusIdx,:) = twoDeltaML(sameStatusIdx,:)...
%                                    + twoDeltaML(sameStatusIdx + phiColNum/2,:);
%     twoDeltaMLPairSum(sameStatusIdx) = sum(twoDeltaMLPair(sameStatusIdx,:),2);                           
%     
%     twoDeltaMLPairSum(diffStatusIdx) =...
%         max([sum(twoDeltaML(diffStatusIdx,:),2),sum(twoDeltaML(diffStatusIdx + phiColNum/2,:),2),...
%         sum(twoDeltaML(diffStatusIdx,:)+twoDeltaML(diffStatusIdx + phiColNum/2,:),2)],[],2);  
    twoDeltaMLPairSum = sum(twoDeltaML(1:phiColNum/2,:) + twoDeltaML(phiColNum/2+1:end,:),2); 
    % the max ML of pair    
    [maxML(i),maxDeltaMLIdx] = max(twoDeltaMLPairSum);  
    mlPlot(i) = mlPlot(i-1) + maxML(i);
    % terminate ?
    if maxML(i) == -Inf
        break;
    end
    if (i > 2) && (abs(maxML(i) - maxML(i-1)) < (abs(maxML(i) - maxML(1)) * eta))
        break;
    end
    
    % primaryPairIdx are from the likelihood
    % assiciatedPairIdx is from the pair (its delta ML is -INF)
    primaryPairIdx = []; associatedPairIdx = [];
%     if (ismember(maxDeltaMLIdx, sameStatusIdx))
        primaryPairIdx = maxDeltaMLIdx;
        primaryPairIdx = [primaryPairIdx, maxDeltaMLIdx + phiColNum/2];
%     else
%         if (twoDeltaML(maxDeltaMLIdx) > twoDeltaML(maxDeltaMLIdx + phiColNum/2))
%             primaryPairIdx = maxDeltaMLIdx;
%             associatedPairIdx = maxDeltaMLIdx + phiColNum/2;
%         else
%             associatedPairIdx = maxDeltaMLIdx;
%             primaryPairIdx = maxDeltaMLIdx + phiColNum/2;
%         end
%     end 
    
    % the primary element of the pair
    % update the alpha, mu and covariance and (S, Q)
    for primaryBasis = 1 : length(primaryPairIdx)
        maxMLIndex = primaryPairIdx(primaryBasis);
        alphaIndex = find(selectedIndex == maxMLIndex);
        % remove it from the associated index vectoy if it should be added
        assoIndex = find(selectedIdxAssociated == maxMLIndex);
        if ismember(maxMLIndex, reestimateIndex)
            if testMode
                disp(['Iteration ', num2str(i), ': Reestimate']);
            end
            for snapShot = 1 : L
                deltaAlpha = newAlpha(maxMLIndex) - alpha(alphaIndex);
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
            alpha(alphaIndex) = newAlpha(maxMLIndex);
        % add basis column
        elseif ismember(maxMLIndex, addIndex)
            if testMode
                disp(['Iteration ', num2str(i), ': Add']);
            end
            for snapShot = 1 : L
                covII = inv(newAlpha(maxMLIndex) + s(maxMLIndex,snapShot));
                muAdd = Q(maxMLIndex,snapShot)/(newAlpha(maxMLIndex) + s(maxMLIndex,snapShot));
                phiI = phi{snapShot}(:,maxMLIndex);
                
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
            % update 10/10/2017 
            % if the basis should be added now, remove the corresponding
            % indice from the selectedIdxAssociated array
            if ~isempty(assoIndex)
              selectedIdxAssociated(assoIndex) = [];
            end
            selectedIndex = [selectedIndex;maxMLIndex];
            alpha = [alpha;newAlpha(maxMLIndex)];
            % delete existing basis
        elseif ismember(maxMLIndex,deleteIndex)
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
    
    % the associated element of the pair
    for assoBasis = 1 : length(associatedPairIdx)
        maxMLIndex = associatedPairIdx(assoBasis);
        alphaIndex = find(maxMLIndex == selectedIndex);
        assoIndex = find(maxMLIndex == selectedIdxAssociated);
        if ismember(maxMLIndex,reestimateIndex)
           % equations (32) - (36)
           if testMode
                disp(['Asso. Iteration : ', num2str(i), ': Reestimate']);
           end
           for snapShot = 1 : L
                deltaAlpha = newAlpha(maxMLIndex) - alpha(alphaIndex);
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
            alpha(alphaIndex) = newAlpha(maxMLIndex);
           % update 10/10/2017 
           % mark the corresponding in the associated slected index
           % because the primary element should be deleted
           selectedIdxAssociated = [selectedIdxAssociated,primaryPairIdx];
        elseif ismember(maxMLIndex,addIndex)
           % equations (27) - (31)
           if testMode
               disp(['Asso. Iteration ', num2str(i), ': Add']);
           end
           for snapShot = 1 : L
                covII = 1/(newAlpha(maxMLIndex) + S(maxMLIndex,snapShot));
                muAdd = covII*Q(maxMLIndex,snapShot);
                phiI = phi{snapShot}(:,maxMLIndex);
                
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
           selectedIndex = [selectedIndex;maxMLIndex];
           % update 10/10/2017 
           % mark the primary element in the associated slected index
           selectedIdxAssociated = [selectedIdxAssociated,primaryPairIdx];
           alpha = [alpha;newAlpha(maxMLIndex)];
        elseif ismember(maxMLIndex,deleteIndex)
           % equations (37) - (41)
           if testMode
                disp(['Asso. Iteration ', num2str(i), ': Delete']);
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
           % update 10/10/2017 
           selectedIndex(alphaIndex) = [];
           % mark the current indice in the associated slected index
           selectedIdxAssociated = [selectedIdxAssociated,maxMLIndex];
           alpha(alphaIndex) = [];
        else
            % update 10/10/2017
            % mark the current indice in the associated slected index
            if isempty(assoIndex)
              selectedIdxAssociated = [selectedIdxAssociated,maxMLIndex];
            end
        end
    end   
end

% to estimate the noise variance 
numeratorSum = 0;
denominatorSum = 0;
for i = 1 : L
    numeratorSum = numeratorSum + sum((y{i} - phiSelected{i}*mu{i}).^2);
    denominatorSum = denominatorSum + (Row(i) - length(selectedIndex) + alpha'*diag(covMatrix{i}));
end
recSigmaSquare = numeratorSum/denominatorSum;

time = toc;

weights = zeros(phiColNum,L);
recoveredValuesSD = zeros(length(selectedIndex),1);
for i = 1 : L
    weights(selectedIndex,i) = mu{i}; 
    recoveredValuesSD = recoveredValuesSD + sqrt(diag(covMatrix{i}))/L;
end


% plot part used for slides
% plot(mlPlot,'o-','LineWidth',1.5)
% xlabel('Iteration','interpreter','latex','fontsize',15)
% ylabel('Marginal likelihood','fontsize',15)
% set(gca,'fontsize',15)
% figure()
% plotAlpha = ones(size(phi{1},2),1);
% plotAlpha(selectedIndex) = alpha;
% plot(1:size(phi{1},2),plotAlpha,'color','red','LineWidth',2);
% %xlabel('i','interpreter','latex','fontsize',15)
% ylabel('$1/ \alpha$','interpreter','latex','fontsize',15)
% set(gca,'fontsize',15)
end
