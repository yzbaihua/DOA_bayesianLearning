function [time,weights,used,sigmaSquare,errbars,basis,selected,alpha,lambda] = ...
    MultiTaskFastLaplacePair(phi,y,sigmaSquare,eta,lambda_init, adaptive,optimal,scale, verbose)

%Check inputs
if nargin < 2
    error('Not enough inputs');
end
if nargin < 3
    sigmaSquare =  std(y)^2/1e2;
end
if nargin < 4
    eta = 1e-8;
end
if nargin < 5
    lambda_init = [];
end
if nargin < 6
    adaptive = 0;
end
if nargin < 7
    optimal = 1;
end
if nargin < 8
    scale = 0.1;
end
if nargin < 9
    verbose = 0;
end

if iscell(y)
    L = length(y); % number of snapshots
else
    L = 1;
end

for i = 1 : L
    [Row(i), Col(i)] = size(phi{i});
end
if sum(abs(diff(Col))) ~= 0 
    error('The basis matrix should have same columns.');
else
    phiColNum = Col(1);
end

% calculate the phiInnerProd and phiYprod
for i = 1 : L
   phiInnerPord(:,i) = sum(phi{i}.*phi{i})'; 
   phiYProd(:,i) = phi{i}'*y{i};
end

% only for initial likelihood 
ml = phiYProd.^2./phiInnerPord; 
ml_sum = sum(ml(1:phiColNum/2,:) + ml(phiColNum/2+1:end,:),2);

while 1
    [maxML, selectedIndex] = max(ml_sum);
    mlPlot(1) = maxML;
    alphaReal = L./...
        sum((phiYProd(selectedIndex,:).^2./phiInnerPord(selectedIndex,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex,:));
    alphaImag = L./...
        sum((phiYProd(selectedIndex + phiColNum/2,:).^2./phiInnerPord(selectedIndex + phiColNum/2,:) - sigmaSquare)./...
        phiInnerPord(selectedIndex + phiColNum/2,:));
    if (alphaReal > 0) && (alphaImag > 0)
        selectedIndex = [selectedIndex; selectedIndex + phiColNum/2];
        alpha = diag([alphaReal,alphaImag]);
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
%     left = phi{i}'*phiSelected{i}/sigmaSquare;
    phiBYProd = phi{i}' * y{i} / sigmaSquare;    
    for j  =  1 : phiColNum
        S(j,i) = phiInnerPord(j,i)/sigmaSquare - ...
            phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * phi{i}(:,j) / sigmaSquare;
        Q(j,i) = phiBYProd(j) -...
            phi{i}(:,j)' * phiSelected{i} / sigmaSquare * covMatrix{i} * phiSelected{i}' * y{i} / sigmaSquare;
    end
    
%     S(:,i) = phiInnerPord(:,i)/sigmaSquare - covMatrix{i}*(left).^2;
%     Q(:,i) = phiBYProd - covMatrix{i}*phiBYProd(selectedIndex)*left;
end


alpha = diag(alpha);
selected = [selectedIndex; selectedIndex + phiColNum/2];
deleted = [];
max_it = 10000;

tic
for count = 1:max_it

    s = S; q = Q;
    alphaRow = repmat(alpha,1,L);
    s(selectedIndex,:) = alphaRow.*S(selectedIndex,:)./(alphaRow-S(selectedIndex,:));
    q(selectedIndex,:) = alphaRow.*Q(selectedIndex,:)./(alphaRow-S(selectedIndex,:));

    if isempty(lambda_init)
        lambda = 2*( length(selectedIndex) - 1 ) / sum(1./alpha);
    else
        lambda = lambda_init;
    end
   
        % old part 
%     A = L*lambda + sum(s,2) - sum(q.^2,2);
%     B = 2*lambda.*sum(s,2) + sum(s.^2,2);
%     C = lambda.*sum(s.^2,2);
%     theta = sum(q.^2-s,2);
%     ig0 = find(theta > L*lambda);
%     is0 = find(theta < L*lambda);
    % end of old
    % new part
    A = lambda*sum(1./(s.^2),2) + sum(1./s,2) - sum((q.^2)./(s.^2),2);
    B = L + 2*lambda*sum(1./s,2);
    C = L*lambda;
    theta = sum((q.^2)./(s.^2),2) - sum(1./s,2);
    ig0 = find(theta > lambda*sum(1./(s.^2),2));
    is0 = find(theta < lambda*sum(1./(s.^2),2));
    % end of new part
    
    discriminant = B.^2 - 4.*A.*C;
    
    newAlpha = (-B - sqrt(discriminant) ) ./ (2*A);

    % choose the next alpha that maximizes marginal likelihood
    ml = -inf*ones(phiColNum,L);
   
    % indices for reestimation
    [ire,foo,which] = intersect(ig0,selectedIndex);
    if ~isempty(ire)
        newAlphaRow = repmat(newAlpha(ire),1,L);
        oldAlphaRow = repmat(alpha(which),1,L);
        ml(ire,:) = q(ire,:).^2./ (newAlphaRow + s(ire,:)) + log(newAlphaRow ./ (newAlphaRow + s(ire,:))) ...
                   -lambda./ newAlphaRow ...
                   -q(ire,:).^2./ (oldAlphaRow + s(ire,:)) - log(oldAlphaRow ./ (oldAlphaRow + s(ire,:)))...
                   +lambda./ oldAlphaRow;
    end
    
    % indices for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        addAlphaRow = repmat(newAlpha(iad),1,L);
        ml(iad,:) = log(addAlphaRow ./ (addAlphaRow + s(iad,:)) ) + ...
            q(iad,:).^2 ./ (addAlphaRow + s(iad,:)) - lambda./addAlphaRow;
        which = intersect(deleted,iad);
        ml(which,:) = -inf;
        
    end
    % indices for deleting
    [ide,foo,which] = intersect(is0,selectedIndex);
    if ~isempty(ide)
         oldAlphaRow = repmat(alpha(which),1,L);
         if length(selectedIndex) == 1
             ml(ide,:) = -inf;
         else
             ml(ide,:) = -q(ide,:).^2 ./ (oldAlphaRow + s(ide,:)) ...
                 - log( oldAlphaRow ./(oldAlphaRow + s(ide,:))) + lambda./oldAlphaRow;
         end

    end
    twoDeltaMLPairSum = sum(ml(1:phiColNum/2,:) + ml(phiColNum/2+1:end,:),2); 
    [ML(count),maxDeltaMLIdx] = max(twoDeltaMLPairSum);
     mlPlot(count + 1) = mlPlot(count) + ML(count);
     % check convergence
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta
        break;
    end
    primaryPairIdx = [maxDeltaMLIdx,maxDeltaMLIdx + phiColNum/2];
    % update alphas
    % Choose the basis which results in the largest increase in the
    % likelihood
    for primaryBasis = 1 : length(primaryPairIdx)
        maxMLIndex = primaryPairIdx(primaryBasis);
        alphaIndex = find(selectedIndex==maxMLIndex);
        if ismember(maxMLIndex, ire) % reestimate a basis
            deltaAlpha = newAlpha(maxMLIndex) - alpha(alphaIndex);
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
            alpha(alphaIndex) = newAlpha(maxMLIndex);      
            if verbose,  fprintf(2,'Reestimate %d..\n',maxMLIndex);end

        elseif ismember(maxMLIndex, iad)

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
            alpha = [alpha;newAlpha(maxMLIndex)];
            if verbose,  fprintf(2,'Add %d.. \n',maxMLIndex); end
        elseif ismember(maxMLIndex,ide)
           deleted = [deleted maxMLIndex];
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
            if verbose,  fprintf(2,'Delete %d.. \n',maxDeltaMLIdx); end
        end    
        selected = [selected; maxDeltaMLIdx];
    end
end
time = toc;
weights	= zeros(phiColNum,L);
used = selectedIndex;
% re-estimated sigma2
numeratorSum = 0;
denominatorSum = 0;
for i = 1 : L
    weights(selectedIndex,i) = mu{i};
    numeratorSum = numeratorSum + sum((y{i} - phiSelected{i}*mu{i}).^2);
    denominatorSum = denominatorSum + (Row(i) - length(selectedIndex) + alpha'*diag(covMatrix{i}));
end
sigmaSquare = numeratorSum/denominatorSum; 

meanCov = zeros(length(selectedIndex),length(selectedIndex));
errbars = zeros(length(selectedIndex),length(selectedIndex));
for i = 1 : L
    meanCov = meanCov + covMatrix{i};
    errbars = errbars + sqrt(diag(covMatrix{i}))/L;
end
%   errbars = sqrt(diag(meanCov));

basis = [];
% generate a basis for adaptive CS?
if adaptive
    if optimal
        [V,D] = eig(Sig);
        [foo,idx] = max(diag(D));
        basis = V(:,idx)';
    else
        temp = phi'*phi/sigmaSquare;
        Sig_inv = temp + scale*mean(diag(temp))*eye(length(used));
        [V,D] = eig(Sig_inv);
        [foo,idx] = min(diag(D));
        basis = V(:,idx)';
    end
end
%fprintf(1,'Algorithm converged, # iterations : %d \n',count);



