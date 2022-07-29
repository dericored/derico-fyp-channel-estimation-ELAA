function[hhat, hhat_n, hhat_f, posn, posf, n, f]=Hybrid_SAMP(y,W,Uf,Un,s,SNR,gamma)

% 1. Perform SAMP for near-field components in polar domain
% 2. Perform SAMP for far-field components in angle domain
% 3. Hybrid channel reconstruction

An = W*Un;
Af = W*Uf;
r = y;  % Initialization of residual

[M,N] = size(An);
posn = []; % Final set
In = s; % Active set size
n = 1; % Stage
epsilon_n = norm(y)^2/(SNR+1)*(1+gamma*SNR); % Threshold
% Near-field estimation
while (norm(r)^2>epsilon_n && In<=N)
    % Candidate list
    [~, idx] = sort(abs(An'*r), 'descend');       
    Cn = union(posn, idx(1:In));
    
    % Finalist
    [~, idx] = sort(abs(pinv(An(:,Cn))*y), 'descend');
    posn_new = Cn(idx(1:In));
    hhat_n = zeros(N,1);
    hhat_n(posn_new) = pinv(An(:,posn_new))*y;
    r_new = y-An*hhat_n;
    if (norm(r_new)>=norm(r))
        % Switch to new stage
        n = n+1;
        In = n*s;    
    else          
        % Update residual and active set    
        r = r_new;
        posn = posn_new;    
    end
end

posf=[];
[M,N]=size(Af);
If = s; % Active set size
f = 1; % Stage
epsilon_f = norm(y)^2/(SNR+1); % Threshold
if ~isempty(posn)
    y = y - An*hhat_n;
end
% Far-field estimation
while (norm(r)^2>epsilon_f && If<=N)
    % Candidate list
    [~, idx] = sort(abs(Af'*r), 'descend');       
    Cf = union(posf, idx(1:If));
    
    % Finalist
    [~, idx] = sort(abs(pinv(Af(:,Cf))*y), 'descend');
    posf_new = Cf(idx(1:If));
    hhat_f = zeros(N,1);
    hhat_f(posf_new) = pinv(Af(:,posf_new))*y;
    r_new = y - Af*hhat_f;
    if (norm(r_new)>=norm(r))
        % Switch to new stage
        f = f+1;
        If = f*s;    
    else          
        % Update residual and active set    
        r = r_new;
        posf = posf_new;    
    end
end

% Hybrid channel reconstruction
if ~isempty(posn)
    if ~isempty(posf)
        hhat = Uf*hhat_f+Un*hhat_n;
    else
        hhat = Un*hhat_n;
    end
else
    hhat = Uf*hhat_f;
end
