clear all
close all
clc
parpool(2)
tic
% % NumOfNeurons = 400;
% g = 2.0;
mu = 0.0;
g_target=1.0;

N = [2000];
N_range=1:length(N);
tolerance=10^(-5);
num_steps=50000;

tic
eta_tab=[-0.9:0.1:0.9];



geff=1.5;
g=geff./(1+eta_tab);


rlzTot = 2;
rhoMat = zeros(length(N),length(eta_tab));
ww=0;
dim_mat=zeros(length(N),length(eta_tab),rlzTot);
var_per_mat=zeros(length(N),length(eta_tab),rlzTot);
std_tab=zeros(length(N),length(eta_tab),rlzTot);
std_time=zeros(length(N),length(eta_tab),rlzTot,num_steps);
PR_mat=zeros(length(N),length(eta_tab),rlzTot);
leading_eigenvalue=zeros(length(N),length(eta_tab),rlzTot);
LyapExp_tab=zeros(length(N),length(eta_tab),rlzTot);
LyapExp_FP=zeros(length(N),length(eta_tab),rlzTot);

TC=zeros(length(N),rlzTot,length(eta_tab));
frac_unstable=zeros(length(N),rlzTot,length(eta_tab));
frac_outside_circle=zeros(length(N),rlzTot,length(eta_tab));
No_fp=zeros(length(N), length(eta_tab), rlzTot);
Time_Of_Convergence=zeros(length(N), length(eta_tab), rlzTot);
Path_Length = zeros(length(N), length(eta_tab), rlzTot);
max_eig=zeros(length(N),rlzTot,length(eta_tab));

AngGeom_median      = zeros(length(N), length(eta_tab), rlzTot);  % median(sin^2 θ)
S_tang_over_total   = zeros(length(N), length(eta_tab), rlzTot);  % ∫||v_perp|| / ∫||v||
WindingProxy        = zeros(length(N), length(eta_tab), rlzTot);  % ∫ ω(t) dt proxy

RadRate_median     = zeros(length(N), length(eta_tab), rlzTot);   % median gamma_H
H_tang_over_total  = zeros(length(N), length(eta_tab), rlzTot);   % ∫||v_{H,⊥}|| / ∫||v||
ThetaTrueIntegral  = zeros(length(N), length(eta_tab), rlzTot);   % ∫ (||v_perp|| / ||x||) dt
Response= zeros(length(N),length(eta_tab),rlzTot,num_steps);
Response_sq= zeros(length(N),length(eta_tab),rlzTot,num_steps);
eigJ_th_ex= zeros(length(N),length(eta_tab),rlzTot,4);
lam_J=zeros(1,num_steps);
eigJ_pairs= zeros(length(N),length(eta_tab),rlzTot,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for NumIdx =N_range
    NumOfNeurons = N(NumIdx);

% disp(NumIdx)
 sprintf('N=%d',NumOfNeurons)


eta_ind=0;
for eta = eta_tab
%     disp(k_ind)
eta_ind=eta_ind+1;
     sprintf('eta_ind=%d',eta_ind)

   sigma = sqrt(g(eta_ind)^2/ NumOfNeurons);
% WMat =sigma.*randn(NumOfNeurons,NumOfNeurons) + mu;

   parfor rlzvar = 1:rlzTot
%      for rlzvar = 1:rlzTot

           disp(rlzvar)

 Mat = randn(NumOfNeurons);
 Msymm=(Mat+Mat.')/sqrt(2);
 Masymm = (Mat-Mat.')/sqrt(2);
 MtotPair = @(eta) (sqrt((1+eta)/2)*(Msymm)+ sqrt((1-eta)/2) * (Masymm));
 WMat=sigma*MtotPair(eta);
  

  Mat_power=WMat;
%   Mat_power=squeeze(Mat_power);
trace_Mat_power(1,rlzvar)=trace(mpower(Mat_power,2));

rho_tmp=trace(mpower(Mat_power,2))/NumOfNeurons;
rho=rho_tmp;
  EigenTable=eig(Mat_power);
  g_factor=max(real(EigenTable));
  Mat_power=(geff / g_factor)*Mat_power;
  EigenTable=eig(Mat_power);





A=Mat_power;
leading_eigenvalue(NumIdx,eta_ind,rlzvar)=max(real(eig(A)));
% p0 = -0.5+1*rand(NumOfNeurons,1);
% t_span = [0,1000];
% [t,P] = ode45(@(t,p) myode(t,p,A),t_span,p0);

dt=0.1;
x=zeros(num_steps,NumOfNeurons);
y=zeros(num_steps,NumOfNeurons);


x(1,:)=randn(1,NumOfNeurons);
y(1,:)=x(1,:)+0.01*rand;
epsilon=norm(y(1,:)-x(1,:));

tr=500;

lypexp=0;
firstZeroTime = NaN;   % SENTINEL: will remain NaN if convergence time not detected

for i = 1:num_steps

%     Phi_prime = diag(1 - tanh(x(i,:)).^2);
%     J = A * Phi_prime;
%     eigs_J = eig(J);
      
%     lam_J(1,i)=max(real(eigs_J));

    x(i+1,:) = (x(i,:)'+(A*tanh(x(i,:)')-x(i,:)')*dt)';

    y(i+1,:) = (y(i,:)'+(A*tanh(y(i,:)')-y(i,:)')*dt)';

    delta= y(i+1,:)-x(i+1,:);

if i>tr
lypexp=lypexp+1/((num_steps-tr)*dt)*log(norm(delta)/epsilon);
end
deltaprime=epsilon*(delta/norm(delta));


y(i+1,:)=x(i+1,:)+deltaprime;

end
Response(NumIdx, eta_ind, rlzvar,:)=mean((sech(x(1:num_steps,:)).^2),2);
% Response(NumIdx, eta_ind, rlzvar,:)=mean((1-tanh(x(1:num_steps,:)).^2),2);
Response_sq(NumIdx, eta_ind, rlzvar,:)=mean((sech(x(1:num_steps,:)).^2).^2,2);
% Response_sq(NumIdx, eta_ind, rlzvar,:)=mean((1-tanh(x(1:num_steps,:)).^2).^2,2);

% figure
% plot(lam_J)
rms_slope = sqrt(mean((1 - tanh(x(end,:)).^2).^2));

Phi_prime = diag(1 - tanh(x(end,:)).^2);
J = A * Phi_prime;
eigs_J = eig(J);
max_real_eig = max(real(eigs_J));
max_imag_eig = max(imag(eigs_J));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================================================================
% === ADDED: Theoretical "Dragonfly" Edge & Split Prediction ===
% ================================================================

% 1. Get the vector of slopes s from the final state
s_dist = 1 - tanh(x(end,:)).^2;

% 2. Calculate the Statistical Moments of the slope
mu_1 = mean(s_dist);      % <s^1> (Mean)
mu_2 = mean(s_dist.^2);   % <s^2> (Used for RMS)
mu_3 = mean(s_dist.^3);   % <s^3> (Skew correction)

% 3. Calculate the Shape Factor F
% This correction factor > 1 accounts for the heavy tails of saturation
F_shape = (mu_1 * mu_3) / (mu_2^2);

% 4. Compute the Theoretical Real Edge 'a' (The Waist)
% Formula: a = g * sqrt(<s^2>) * (1 + tau * F)
% Note: Assuming 'eta' is your variable for tau, and 'g' is an array
tau_val = eta; 
g_val = g(eta_ind);

theory_waist = g_val * sqrt(mu_2) * (1 + tau_val * F_shape);

% 5. Predict Bubble Split
% If the waist 'a' <= 0, the real edge has collapsed (Gap Formation)
if theory_waist <= 0
    theory_split_status = 1; % Split (Dragonfly / Complex Bubbles)
    theory_waist_clamped = 0; % Edge cannot be negative, it vanishes
else
    theory_split_status = 0; % Connected Spectrum
    theory_waist_clamped = theory_waist;
end

% Optional: Store for plotting later
% Theory_Edge_Store(NumIdx, eta_ind, rlzvar) = theory_waist_clamped;
% Split_Status_Store(NumIdx, eta_ind, rlzvar) = theory_split_status;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, idxRight] = max(real(eigs_J));
imag_rightmost = imag(eigs_J(idxRight));

% ---- (2) Real part of the eigenvalue with the upper-most imag part ----
[~, idxTop] = max(imag(eigs_J));
real_topmost_imag = real(eigs_J(idxTop));
eigJ_pairs(NumIdx, eta_ind, rlzvar,:)=[max_real_eig imag_rightmost max_imag_eig real_topmost_imag theory_waist];


term_inside = mean((1 - tanh(x(end,:)).^2).^2,2) + 2*eta*mean((1 - tanh(x(end,:)).^2)).^2;
if term_inside < 0
   % Purely imaginary spectrum
   theory_edge = 0; 
else
   theory_edge = g(eta_ind) * (1 + eta) * rms_slope;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Get the empirical slopes s
s_values = 1 - tanh(x(end,:)).^2;
% Filter out numerical zeros to speed up/stabilize
s_values = s_values(s_values > 1e-6); 
N_eff = length(s_values);

% 2. Define the equations as a function of (a, y)
% We effectively do Monte Carlo integration using the actual data points
% Eq 1: y - mean( s ./ (a - tau*g^2*y*s) ) = 0
% Eq 2: 1 - g^2 * mean( s.^2 ./ (a - tau*g^2*y*s).^2 ) = 0

tau_val = eta; % Assuming your variable is eta
g_val = g(eta_ind);
K = tau_val * g_val^2;

% We solve for intersection. 
% Since this is sensitive, we scan 'a' and solve for 'y' iteratively.

% Define range to search for 'a' (start from naive guess)
a_guess = g_val * (1+tau_val) * sqrt(mean(s_values.^2)); 
a_scan = linspace(0, a_guess*2, 100); 
residuals = zeros(size(a_scan));

for i = 1:length(a_scan)
    a = a_scan(i);
    
    % Solve Eq 1 for y using fzero
    % y must be positive real for the edge
    fun_y = @(y) y - mean(s_values ./ (a - K*y*s_values));
    try
        y_sol = fzero(fun_y, 0.5); % Start guess at 0.5
    catch
        y_sol = NaN;
    end
    
    if ~isnan(y_sol)
        % Check Eq 2 (Stability)
        LHS = 1;
        RHS = g_val^2 * mean(s_values.^2 ./ (a - K*y_sol*s_values).^2);
        residuals(i) = LHS - RHS;
    else
        residuals(i) = NaN;
    end
end

% 3. Find where Eq 2 is satisfied (residual crosses 0)
% The rightmost solution is the edge.
[~, idx] = min(abs(residuals));
exact_edge = a_scan(idx);

eigJ_th_ex(NumIdx, eta_ind, rlzvar,:)=[max_real_eig theory_edge exact_edge max_imag_eig];

%%%%path length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf=length(x(1:end,1));

nr=randi(NumOfNeurons,1);
rt=randi(tf,1);
dx = diff(x(:,nr));

cond=abs(dx) < tolerance;

std_tab(NumIdx, eta_ind, rlzvar) =std(x,0,"all");

std_time(NumIdx, eta_ind, rlzvar,:)=std(x(1:num_steps,:),0,2);
if all(cond==0)

    No_fp(NumIdx, eta_ind, rlzvar)=1;
    firstZeroIndex=num_steps;
    firstZeroTime=firstZeroIndex;
%     std_tab(NumIdx, eta_ind, rlzvar) =std(x,0,"all");


else

zeroDerivativeIndices = find(cond);

firstZeroIndex = zeroDerivativeIndices(1);

for o=firstZeroIndex:tf-1

%     for a single x
xi     = x(o, nr);                          % reference value
tail   = x((o+1):end, nr);                  % all future values
diffs  = abs(tail - xi) < tolerance;       % logical vector


if diffs  %all(diffs, 'all')

    firstZeroIndex = o;

    firstZeroTime = firstZeroIndex;
    remainingDerivatives = dx(firstZeroIndex:end);
  
    Time_Of_Convergence(NumIdx, eta_ind, rlzvar)=firstZeroTime;
    break
elseif  o==tf-1
    No_fp(NumIdx, eta_ind, rlzvar)=1;
end

end

conv_idx = Time_Of_Convergence(NumIdx, eta_ind, rlzvar);


% ---------- Path-length of the transient ------------------------------
if conv_idx>1 && conv_idx < 0.9*num_steps
    % Velocity matrix up to convergence (rows = time, cols = neurons)
    v_seg = diff(x(1:firstZeroTime, :), 1, 1) / dt;  % [ (T_conv-1) Ã— N ]
    
    % Euclidean speed at each step (vector of length T_conv-1)
    speed = vecnorm(v_seg, 2, 2);                    % 2-norm across neurons
    
    % Path length per neuron (divide by sqrt(N) for N-invariance)
    Path_Length(NumIdx, eta_ind, rlzvar) = sum(speed) * dt / sqrt(NumOfNeurons);
%     std_tab(NumIdx, eta_ind, rlzvar) =std(x(conv_idx,:));
else
    % If convergence was not detected, mark as NaN (or 0 if you prefer)
    Path_Length(NumIdx, eta_ind, rlzvar) = NaN;
end

end
% ====== Geometry checks (1) angle/tangential, (2) rotational) on the same transient ======
% if ~isnan(firstZeroTime) && firstZeroTime > 2 && firstZeroTime < 0.9*num_steps
    % Build the same segment used for Path_Length
    firstZeroTime=num_steps;
    Xseg = x(1:firstZeroTime, :);               % T×N
    Vseg = diff(Xseg, 1, 1) / dt;               % (T-1)×N

    nr    = vecnorm(Xseg(1:end-1, :), 2, 2);    % ||x|| per time step
    nr(nr==0) = eps;
    vnorm = vecnorm(Vseg, 2, 2);                % ||dx/dt|| per time step
    vnorm(vnorm==0) = eps;

    % Decompose velocity: radial vs tangential components
    vr    = sum(Xseg(1:end-1, :) .* Vseg, 2) ./ nr;                    % radial speed
    vperp = Vseg - Xseg(1:end-1, :) .* (vr ./ nr);                     % (T-1)×N tangential

    % Angle metric: cosθ between x and dx/dt → sin^2θ = 1 - cos^2θ
    cos_th = sum(Xseg(1:end-1, :) .* Vseg, 2) ./ (nr .* vnorm);
    cos_th = max(min(cos_th, 1), -1);                                  % clamp
    ang_geom = 1 - cos_th.^2;

    % Arc-length style integrals, matching your Path_Length convention
    S_total_seg = sum(vnorm) * dt;                                     % ∫||dx/dt|| dt
    S_tang_seg  = sum(vecnorm(vperp, 2, 2)) * dt;                      % ∫||v_perp|| dt

    AngGeom_median(NumIdx, eta_ind, rlzvar)    = median(ang_geom);
    S_tang_over_total(NumIdx, eta_ind, rlzvar) = S_tang_seg / max(S_total_seg, eps);

    % Rotational proxy via antisymmetric part of A (your Mat_power)
    Aconn = A;                                                          % = Mat_power
    Askew = 0.5 * (Aconn - Aconn.');                                    % antisymmetric
    Ax    = (Askew * Xseg(1:end-1, :).').';                             % (T-1)×N
    omega = vecnorm(Ax, 2, 2) ./ nr;                                    % dimensionless
    t_omega = (0:length(omega)-1).' * dt;
    WindingProxy(NumIdx, eta_ind, rlzvar) = trapz(t_omega, omega);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% if  abs(x(num_steps,1)-x(0.9*num_steps,1))<10^(-4) && abs(x(num_steps,1)-x(0.95*num_steps,1))<10^(-4)

if all(abs(diff(x(0.9*num_steps:num_steps,1)))<10^(-4))

    LyapExp_FP(NumIdx,eta_ind, rlzvar)=-1;

end


   LyapExp_tab(NumIdx,eta_ind, rlzvar)=lypexp;


Cov = cov(x);

Cov_eig=eig(Cov);
Cov_eig_tilde=Cov_eig./sum(Cov_eig);
PR=sum(Cov_eig_tilde)^2/sum(Cov_eig_tilde.^2);

timeSteps = 1:numel(x(:, 1));

if abs(diff(x([round(0.95*numel(timeSteps)),numel(timeSteps)], 1)))<10^(-3)
PR=0;
end

var_per=Cov_eig./sum(Cov_eig);
var_exp=cumsum(var_per);


PR_mat(NumIdx,eta_ind,rlzvar)=PR;
var_per_mat(NumIdx,eta_ind,rlzvar)=numel(var_per(var_per>0.001));
std_mat(NumIdx,eta_ind,rlzvar)=mean(std(x'));
TC(NumIdx,rlzvar,eta_ind)=(1/NumOfNeurons)*(sum(log(abs(eig(A)-1))));
mag_eig=abs(eig(A));
retoteig=real(eig(A));

frac_unstable(NumIdx,rlzvar,eta_ind)=sum(retoteig>1)/NumOfNeurons;
frac_outside_circle(NumIdx,rlzvar,eta_ind)=sum((mag_eig>1))/NumOfNeurons;
max_eig(NumIdx,rlzvar,eta_ind)=max(mag_eig);

    end

   sum_trace_Mat_power=sum(trace_Mat_power(1,:),2);
   rhoMat(NumIdx,eta_ind) =  sum_trace_Mat_power./(g(eta_ind)^2*NumOfNeurons*rlzTot); 

end


end
% 
% % save('leading_eigenvalue_N2000_geff15_1.mat','leading_eigenvalue');
% % save('PR_mat_N2000_geff15_1.mat','PR_mat', '-v7.3');
% % save('rhoMat_N2000_geff15_1.mat','rhoMat', '-v7.3');
% % save('LyapExp_tab_N2000_geff15_1.mat','LyapExp_tab', '-v7.3');
% % save('LyapExp_FP_N2000_geff15_1.mat','LyapExp_FP', '-v7.3');
% % save('TC_rho_N2000_geff15_1.mat','TC', '-v7.3');
% % save('frac_unstable_rho_N2000_geff15_1.mat','frac_unstable', '-v7.3');
% % save('frac_outside_circle_N2000_geff15_1.mat','frac_outside_circle', '-v7.3');
% % save('No_fp_N2000_geff15_1.mat','No_fp', '-v7.3');
% % save('Path_Length_N2000_geff15_1.mat','Path_Length', '-v7.3');
% % save('max_eig_N2000_geff15_1.mat','max_eig', '-v7.3');
% % save('AngGeom_median_N2000_geff15_1.mat','AngGeom_median', '-v7.3');
% % save('S_tang_over_total_N2000_geff15_1.mat','S_tang_over_total', '-v7.3');
% % save('WindingProxy_N2000_geff15_1.mat','WindingProxy', '-v7.3');
% % save('std_tab_N2000_geff15_1.mat','std_tab', '-v7.3');
% % save('std_time_N2000_geff15_1.mat','std_time', '-v7.3');
% % save('Response_N2000_geff15_1.mat','Response', '-v7.3');
% % save('Response_sq_N2000_geff15_1.mat','Response_sq', '-v7.3');
% % save('eigJ_th_ex_N2000_geff15_1.mat','eigJ_th_ex', '-v7.3');
% % save('eigJ_pairs_N2000_geff15_1.mat','eigJ_pairs', '-v7.3');

% % 
toc
