% BVAR_FULL.m 
% This code replicates the results from the 1st empirical illustration 
% in Koop and Korobilis (2009).
% Modified by Zizhe, Xia for the FRB disagreement project
%
% You can chose 6 different priors. For some priors, analytical results are
% available, so Monte Carlo Integration is used. For other priors, you need
% to use the Gibbs sampler. For Gibbs sampler models I take a number of
% 'burn-in draws', so that I keep only the draws which have converged.
%
% The specification of the prior hyperparmeters are in the file
% prior_hyper.m. See there for details.
%
% The convention used here is that ALPHA is the K x M matrix of VAR coefficients,
% alpha is the KM x 1 column vector of vectorized VAR coefficients, i.e.
% alpha = vec(ALPHA), and SIGMA is the M x M VAR covariance matrix.
%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models using posterior simulation. Dependent on your choice of forecasting,
% the VAR model is:
%
% In this code we provide direct (as opposed to iterated) forecasts
% Direct h-step ahead foreacsts:
%     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
%
% so that in this case there are also p lags of Y (from 0 to p-1).
%
% In any of the two cases, the model is written as:
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A).
%--------------------------------------------------------------------------
% NOTES: The code sacrifices efficiency for clarity. It follows the
%        theoretical equations in the monograph and the manual.
%
% AUTHORS: Gary Koop and Dimitris Korobilis
% CONTACT: dikorombilis@yahoo.gr
%--------------------------------------------------------------------------

clear all;
clc;
randn('seed',2); %#ok<RAND>
rand('seed',2); %#ok<RAND>

%------------------------------LOAD DATA-----------------------------------
read_data;
%%% load new data
load Yraw.dat;
load W_T.dat;
load W_S.dat;
W_main = W_T;
% if adding any other exogenous var than the disagreement measure
% the FFR should always be placed as the last dependent variables
% if adding any other exogenous interaction var (e.g. uncertainty) 
% it should be included in the W_main matrix, and the main var to 
% be investigated (disagreement) should be placed first.

%Yraw = Yraw(29:189,:);
% or Simulate data from a simple VAR Data Generating process
%[Yraw] = bvardgp();

% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.

%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
p = 8;               % Number of lags on dependent variables
forecasting = 0;     % 1: Compute h-step ahead predictions, 0: no prediction

repfor = 50;         % Number of times to obtain a draw from the predictive 
                     % density, for each generated draw of the parameters                     
h = 1;               % Number of forecast periods
impulses = 1;        % 1: compute impulse responses, 0: no impulse responses
ihor = 12;           % Horizon to compute impulse responses
full_interaction = 0;% 1: full interaction, 0: only disagreement and effective FFR
other_exo = 0;       % 1: use other exogenous variables, 0: no
exo_included = 0;    % 1: include the standalone exo var in the impulse response, 0: no
w_upper_bar = 100;   % set the exogenous w_upper_bar to be 90th percentile
w_lower_bar = 20;    % set the exogenous w_lower_bar to be 10th percentile
jump = 1;            % 1: sample with jumps to avoid correlation (for final results)

% Set prior for BVAR model:
prior = 4;  % prior = 1 --> Diffuse ('Jeffreys') (M-C Integration)
            % prior = 2 --> Minnesota            (M-C Integration)
            % prior = 3 --> Normal-Wishart       (M-C Integration)           
            % prior = 4 --> Independent Normal-Wishart      (Gibbs sampler)
            % prior = 5 --> SSVS in mean-Wishart            (Gibbs sampler)
            % prior = 6 --> SSVS in mean-SSVS in covariance (Gibbs sampler)
            % irrelevent methods will be deleted here

% Gibbs-related preliminaries
if jump
    njump = 4;            % Number of draws to skip to avoid correlation
    nsave = 4000;
    nburn = 4000;
else
    njump = 1;
    nsave = 8000;         % Final number of draws to save
    nburn = 4000;         % Draws to discard (burn-in)
end

% For models using analytical results, there are no convergence issues (You
% are not adviced to change the next 3 lines)
if prior ==1 || prior == 2 || prior == 3
    nburn = 0*nburn;
end
if jump
    ntot = nburn + nsave * njump;
else
    ntot = nsave + nburn;  % Total number of draws
end
it_print = 2000;       % Print on the screen every "it_print"-th iteration
% get other exogenous variables
if other_exo
    try
        load W_others.dat;
        [Traw,n_exo] = size(W_others);
    catch
        other_exo = false;
        n_exo = 0;
    end
end
%load W_others.dat
%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw, M] = size(Yraw);
[Traw, w] = size(W_main);

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0 % Check for wrong (incompatible) user input
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end
    % Now create VAR specification according to forecast method
   
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
   
else
   Y1 = Yraw;
   Y2 = Yraw;
end
        
% Generate lagged Y matrix. This will be part of the X matrix
% The Ylag matrix will include the dependent variables and the interaction
% term between y and w. Only the interaction between y1(interest rate) and
% w is allowed for the simple interaction. Note that w itself is exogenous, 
% and other exogenous variables are added later.
Ylag = mlag2_interact(Y2,p,W_main,full_interaction); 
%   Y is [T x M]. 
%   Ylag with simple interaction is [T x (M*p+w+p*w)]
%   Ylag with full interaction is [T x (M*p+w+M*p*w)]
%   where w is the number of interaction variable in W_main

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies).
% Note that in this example I do not include exogenous variables (other macro
% variables, dummies, or trends). You can load a file with exogenous
% variables, call them, say W, and then extend variable X1 in line 133, as:
%            X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:) W(p+1:Traw,:)];
% original X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
% and line 135 as:
%            X1 = [Ylag(p+1:Traw,:)  W(p+1:Traw,:)];
% original X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
if other_exo
    if constant
        X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:) W_others(p+1:Traw,:)];
    else
        X1 = [Ylag(p+1:Traw,:)  W_others(p+1:Traw,:)];  %#ok<UNRCH>
    end
else
    if constant
        X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
    else
        X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
    end
end
% Get size of final matrix X
[Traw3, K] = size(X1);

% Create the block diagonal matrix Z
Z1 = kron(eye(M),X1);

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if forecasting==1
    Y_pred = zeros(nsave*repfor,M); % Matrix to save prediction draws
    PL = zeros(nsave,1);            % Matrix to save Predictive Likelihood
    
    % Direct forecasts, we only need to keep the last observation for evaluation
        Y = Y1(1:end-1,:);          
        X = X1(1:end-1,:);
        Z = kron(eye(M),X);
        T = T - 1;
   
else % if no prediction is present, keep all observations
    Y = Y1;
    X = X1;
    Z = Z1;
end

%========= IMPULSE RESPONSES SET-UP:
% Create matrices to store forecasts`
if impulses == 1
    imp_rgdp1 = zeros(nsave,M,ihor); % impulse responses to a shock in rGDP
    imp_infl1 = zeros(nsave,M,ihor);  % impulse responses to a shock in inflation (CPI)
    imp_une1 = zeros(nsave,M,ihor);  % impulse responses to a shock in unemployment
    imp_int1 = zeros(nsave,M,ihor); % impulse responses to a shock in the interest rate (effective FFR)
    % same matrices, for upper and lower w_bar results
    imp_rgdp2 = zeros(nsave,M,ihor);
    imp_infl2 = zeros(nsave,M,ihor);
    imp_une2 = zeros(nsave,M,ihor); 
    imp_int2 = zeros(nsave,M,ihor); 
    
    bigj = zeros(M,M*p);
    bigj(1:M,1:M) = eye(M);
end

%-----------------------------PRELIMINARIES--------------------------------
% First get ML estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of parameters, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);   % Sum of squared errors
SIGMA_OLS = SSE./(T-K+1);

% Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;     % This is the single draw from the posterior of alpha
ALPHA = A_OLS;     % This is the single draw from the posterior of ALPHA
SSE_Gibbs = SSE;   % This is the SSE based on each draw of ALPHA
SIGMA = SIGMA_OLS; % This is the single draw from the posterior of SIGMA
IXY =  kron(eye(M),(X'*Y));


% Storage space for posterior draws
alpha_draws = zeros(nsave,K*M);   % save draws of alpha
ALPHA_draws = zeros(nsave,K,M);   % save draws of alpha
SIGMA_draws = zeros(nsave,M,M);   % save draws of ALPHA

%-----------------Prior hyperparameters for bvar model
% load file which sets hyperparameters for chosen prior
prior_hyper;
%-------------------- Prior specification ends here
    
%========================== Start Sampling ================================
%==========================================================================
tic;
disp('Number of iterations');
for irep = 1:ntot  %Start the Gibbs "loop"
    if mod(irep,it_print) == 0 % print iterations
        disp(irep);
        toc;
    end
    
    %--------- Draw ALPHA and SIGMA with Independent Normal-Wishart Prior
    if prior == 4
        VARIANCE = kron(inv(SIGMA),eye(T));
        V_post = inv(V_prior + Z'*VARIANCE*Z);
        a_post = V_post*(V_prior*a_prior + Z'*VARIANCE*Y(:));
        alpha = a_post + chol(V_post)'*randn(n,1); % Draw of alpha
        
        ALPHA = reshape(alpha,K,M); % Draw of ALPHA
        
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
        v_post = T + v_prior;
        S_post = S_prior + (Y - X*ALPHA)'*(Y - X*ALPHA);
        SIGMA = inv(wish(inv(S_post),v_post));% Draw SIGMA
    
    end
    % =============Estimation ends here
    
    
    
    % ****************************|Predictions, Responses, etc|***************************
    if irep > nburn && mod(irep,njump) == 0
        %=========FORECASTING:
        if forecasting==1
          
                Y_temp = zeros(repfor,M);
                % compute 'repfor' predictions for each draw of ALPHA and SIGMA
                for ii = 1:repfor
                    X_fore = [1 Y(T,:) X(T,2:M*(p-1)+1)]; % specify the value of the exogenous disagreement here
                    % Forecast of T+1 conditional on data at time T
                    Y_temp(ii,:) = X_fore*ALPHA + randn(1,M)*chol(SIGMA);
                end
                % Matrix of predictions
                Y_pred(((irep-nburn)-1)*repfor+1:(irep-nburn)*repfor,:) = Y_temp;
                % Predictive likelihood
                PL(irep-nburn,:) = mvnpdf(Y1(T+1,:),X(T,:)*ALPHA,SIGMA);
                if PL(irep-nburn,:) == 0
                    PL(irep-nburn,:) = 1;
                end
          
        end % end forecasting
        %=========Forecasting ends here
        
        %=========IMPULSE RESPONSES:
        if impulses==1
%             % Previous IRF using std dev shocks
%             % ------------Identification code I:
%             % modified to capture the impulse response of the interest rate
%             % (effective FFR) during different levels of disagreement
%             % measure. Use both the st dev matrix decomposition and the
%             % unit variance shock??? (or only the former)
%             Bv = zeros(M,M,p);
%             for i_1=1:p
%                 Bv(:,:,i_1) = ALPHA(1+((i_1-1)*M + 1):i_1*M+1,:);
%             end

            % Impulse response for a 1% increase of effective FFR
            %%%%%%%%%% Beware of b_0 - the coefficent of w_t - without any
            %%%%%%%%%% interaction, pure exogenous variable
            if full_interaction
                By = zeros(M,M,p); % coef matrix for dependent var
                Bx = zeros(M,M,p); % coef matrix for interaction
                if other_exo % size=(w+n_exo,M) coef matrix for exo var
                    Bw = [ALPHA(1+M*p+1:1+M*p+w,:);ALPHA(1+M*p+w+M*p*w+1:1+M*p+w+M*p*w+n_exo,:)];
                else
                    Bw = ALPHA(1+M*p+1:1+M*p+w,:);
                end
                for i=1:p
                    By(:,:,i) = ALPHA(1+((i-1)*M+1):i*M+1,:);
                    Bx(:,:,i) = ALPHA(1+M*p+1+(i-1)*M+1:1+M*p+1+i*M);
                end
            else
                By = zeros(M,M,p);
                Bx = zeros(M,1,p);
                if other_exo % coef matrix for exo var
                    Bw = [ALPHA(1+M*p+1:1+M*p+w,:);ALPHA(1+M*p+w+p*w+1:1+M*p+w+p*w+n_exo,:)];
                else
                    Bw = ALPHA(1+M*p+1:1+M*p+w,:);
                end
                for i=1:p
                    By(:,:,i) = ALPHA(1+((i-1)*M+1):i*M+1,:);
                    Bx(:,:,i) = ALPHA(1+M*p+1+i:1+M*p+1+i,:);
                end
            end
            
            % st dev matrix for structural VAR
            chol_shock = chol(SIGMA)';
            d = diag(diag(chol_shock));
            chol_shock = inv(d)*chol_shock;
            [responses_lower]=impulse_interact_cholesky(By,Bx,Bw,chol_shock,ihor,w_lower_bar,full_interaction,exo_included);
            [responses_upper]=impulse_interact_cholesky(By,Bx,Bw,chol_shock,ihor,w_upper_bar,full_interaction,exo_included);
            
%             simple_shock = eye(4);
%             [responses]=impulse_interact_cholesky(By,Bx,simple_shock,ihor,w_bar,full_interaction);
            
            % Restrict to policy shocks
            % changed to accommodate more shocks
%             responses1 = squeeze(responses_upper(:,1,:));
%             responses2 = squeeze(responses_upper(:,2,:));
%             responses3 = squeeze(responses_upper(:,3,:));
%             responses4 = squeeze(responses_upper(:,4,:));
            
            imp_rgdp1((irep-nburn)/njump,:,:) = squeeze(responses_upper(:,1,:));
            imp_infl1((irep-nburn)/njump,:,:) = squeeze(responses_upper(:,2,:));
            imp_une1((irep-nburn)/njump,:,:) = squeeze(responses_upper(:,3,:));
            imp_int1((irep-nburn)/njump,:,:) = squeeze(responses_upper(:,4,:));
            
            imp_rgdp2((irep-nburn)/njump,:,:) = squeeze(responses_lower(:,1,:));
            imp_infl2((irep-nburn)/njump,:,:) = squeeze(responses_lower(:,2,:));
            imp_une2((irep-nburn)/njump,:,:) = squeeze(responses_lower(:,3,:));
            imp_int2((irep-nburn)/njump,:,:) = squeeze(responses_lower(:,4,:));

        end
               
        %----- Save draws of the parameters
        alpha_draws((irep-nburn)/njump,:) = alpha;
        ALPHA_draws((irep-nburn)/njump,:,:) = ALPHA;
        SIGMA_draws((irep-nburn)/njump,:,:) = SIGMA;

    end % end saving results
       
end %end the main Gibbs for loop
%====================== End Sampling Posteriors ===========================

%Posterior mean of parameters:
ALPHA_mean = squeeze(mean(ALPHA_draws,1)); %posterior mean of ALPHA
SIGMA_mean = squeeze(mean(SIGMA_draws,1)); %posterior mean of SIGMA

%Posterior standard deviations of parameters:
ALPHA_std = squeeze(std(ALPHA_draws,1)); %posterior std of ALPHA
SIGMA_std = squeeze(std(SIGMA_draws,1)); %posterior std of SIGMA

%or you can use 'ALPHA_COV = cov(alpha_draws,1);' to get the full
%covariance matrix of the posterior of alpha (of dimensions [KM x KM] )



% mean prediction and log predictive likelihood
if forecasting == 1
    Y_pred_mean = mean(Y_pred,1); % mean prediction
    Y_pred_std = std(Y_pred,1);   % std prediction
    log_PL = mean((log(PL)),1);

    %This are the true values of Y at T+h:
 
        true_value = Y1(T+1,:);
    
    %(subsequently you can easily also get MSFE and MAFE)
 
    %======PLOT posterior predictive
    %% names already changed
    figure
    bars = 3000;
    subplot(3,1,1)
    hist(Y_pred(:,1),bars);
    title('Inflation')
    text(Y_pred_mean(:,1),max(hist(Y_pred(:,1),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,1)) ', std = ' num2str(Y_pred_std(:,1))],...       
        'HorizontalAlignment','left')
    subplot(3,1,2)
    hist(Y_pred(:,2),bars);
    title('Unemployment')
    text(Y_pred_mean(:,2),max(hist(Y_pred(:,2),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,2)) ', std = ' num2str(Y_pred_std(:,2))],...
        'HorizontalAlignment','left')
    subplot(3,1,3)
    hist(Y_pred(:,3),bars);
    title('Interest Rate')
    text(Y_pred_mean(:,3),max(hist(Y_pred(:,3),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,3)) ', std = ' num2str(Y_pred_std(:,3))],...
        'HorizontalAlignment','left')    
end

% You can also get other quantities, like impulse responses
if impulses==1
    
    % Set quantiles from the posterior density of the impulse responses
    qus = [.05, .16, .5, .84, .95];
    alphas = [.32 .05];
    [imp_resp_rgdp1, imp_resp_rgdp2, pbound_resp_rgdp, imp_resp_diff_qus_rgdp] = posterior_analysis(imp_rgdp1, imp_rgdp2, qus, alphas);
    [imp_resp_infl1, imp_resp_infl2, pbound_resp_infl, imp_resp_diff_qus_infl] = posterior_analysis(imp_infl1, imp_infl2, qus, alphas);
    [imp_resp_une1, imp_resp_une2, pbound_resp_une, imp_resp_diff_qus_une] = posterior_analysis(imp_une1, imp_une2, qus, alphas);
    [imp_resp_int1, imp_resp_int2, pbound_resp_int, imp_resp_diff_qus_int] = posterior_analysis(imp_int1, imp_int2, qus, alphas);

%     figure 
%     plot(squeeze(pbound_resp_int(:,1,:))')
%     plot(squeeze(pbound_resp_int(:,2,:))')
%     plot(squeeze(pbound_resp_int(:,3,:))')

%     plot(squeeze(pbound_resp_int(:,1,:))')
%     imp_resp_rgdp1 = squeeze(quantile(imp_rgdp1,qus));
%     imp_resp_infl1 = squeeze(quantile(imp_infl1,qus));
%     imp_resp_une1 = squeeze(quantile(imp_une1,qus));
%     imp_resp_int1 = squeeze(quantile(imp_int1,qus));
%     
%     imp_resp_rgdp2 = squeeze(quantile(imp_rgdp2,qus));       
%     imp_resp_infl2 = squeeze(quantile(imp_infl2,qus));
%     imp_resp_une2 = squeeze(quantile(imp_une2,qus));
%     imp_resp_int2 = squeeze(quantile(imp_int2,qus));
    
    % Get mean and variance from the posterior and calculate the test
    % statistics using two sample mean test

    %======= Plot impulse responses to FFR
    
    
    % High level disagreement response
    figure
    color_order =  [0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];
    set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder',':|--|-|--|:')%|:|--|-|--|:
    subplot(2,2,1)
    plot(squeeze(imp_resp_rgdp1(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of rGDP, Shock to Interest Rate; High Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,2)
    plot(squeeze(imp_resp_infl1(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Interest Rate; High Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,3)
    plot(squeeze(imp_resp_une1(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Interest Rate; High Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,4)
    plot(squeeze(imp_resp_int1(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Effective FFR, Shock to Effective FFR; High Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    
    % Low level disagreement response
    figure
    color_order =  [0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];
    set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder',':|--|-|--|:')%|:|--|-|--|:
    subplot(2,2,1)
    plot(squeeze(imp_resp_rgdp2(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of rGDP, Shock to Interest Rate; Low Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,2)
    plot(squeeze(imp_resp_infl2(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Interest Rate; Low Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,3)
    plot(squeeze(imp_resp_une2(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Interest Rate; Low Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,4)
    plot(squeeze(imp_resp_int2(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Effective FFR, Shock to Effective FFR; Low Disagreement Level')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    
    % Difference Probability Band
    figure
    color_order =  [0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];
    set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder',':|--|-|--|:')%|:|--|-|--|:
    subplot(2,2,1)
    plot(squeeze(imp_resp_diff_qus_rgdp(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of rGDP, Shock to Interest Rate; Difference Bands')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,2)
    plot(squeeze(imp_resp_diff_qus_infl(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Interest Rate; Difference Bands')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,3)
    plot(squeeze(imp_resp_diff_qus_une(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Interest Rate; Difference Bands')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(2,2,4)
    plot(squeeze(imp_resp_diff_qus_int(:,4,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Effective FFR, Shock to Effective FFR; Difference Bands')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    
%     %======= Plot full impulse responses
%     figure
%     set(0,'DefaultAxesColorOrder',[0 0 0],...
%         'DefaultAxesLineStyleOrder','--|:|-|:|--')
%     subplot(4,4,1)
%     plot(squeeze(imp_resp_rgdp(:,1,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of rGDP, Shock to rGDP')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,2)
%     plot(squeeze(imp_resp_infl(:,1,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Inflation, Shock to rGDP')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,3)
%     plot(squeeze(imp_resp_une(:,1,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Unemployment Rate, Shock to rGDP')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,4)
%     plot(squeeze(imp_resp_int(:,1,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Effective FFR, Shock to rGDP')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     
%     subplot(4,4,5)
%     plot(squeeze(imp_resp_rgdp(:,2,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of rGDP, Shock to Inflation')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,6)
%     plot(squeeze(imp_resp_infl(:,2,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Inflation, Shock to Inflation')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,7)
%     plot(squeeze(imp_resp_une(:,2,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Unemployment Rate, Shock to Inflation')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,8)
%     plot(squeeze(imp_resp_int(:,2,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Effective FFR, Shock to Inflation')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     
%     subplot(4,4,9)
%     plot(squeeze(imp_resp_rgdp(:,3,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of rGDP, Shock to Unemployment')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,10)
%     plot(squeeze(imp_resp_infl(:,3,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Inflation, Shock to Unemployment')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,11)
%     plot(squeeze(imp_resp_une(:,3,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Unemployment, Shock to Unemployment')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,12)
%     plot(squeeze(imp_resp_int(:,3,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Effective FFR, Shock to Unemployment')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     
%     subplot(4,4,13)
%     plot(squeeze(imp_resp_rgdp(:,4,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of rGDP, Shock to Effective FFR')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,14)
%     plot(squeeze(imp_resp_infl(:,4,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Inflation, Shock to Effective FFR')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,15)
%     plot(squeeze(imp_resp_une(:,4,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Unemployment, Shock to Effective FFR')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
%     subplot(4,4,16)
%     plot(squeeze(imp_resp_int(:,4,:))')
%     hold;
%     plot(zeros(1,ihor),'-')
%     title('Response of Effective FFR, Shock to Effective FFR')
%     xlim([1 ihor])
%     set(gca,'XTick',0:4:ihor)
    

    
end



% Clear screen and print elapsed time
clc;
toc;

% Print some directions to the user
disp('Please find the means and variances of the VAR parameters in the vectors')
disp('ALPHA_mean and ALPHA_std for the VAR regression coefficients, and ')
disp('SIGMA_mean and SIGMA_std for the VAR covariance matrix. The predictive')
disp('mean and standard deviation are in Y_pred_mean and Y_pred_std, respectively.')
disp('The log Predictive Likelihood is given by variable log_PL. The true value')
disp('of y(t+h) is given in the variable true_value. For example the mean squared')
disp('forecast error can be obtained using the command')
disp('                MSFE = (Y_pred_mean - true_value).^2')
disp('If you are using the SSVS prior, you can get the averages of the restriction')
disp('indices $\gamma$ and $\omega$. These are in the variables gammas_mat and omegas_mat') 


if prior == 1; name = 'diffuse';
elseif prior == 2; name = 'minnesota';
elseif prior == 3; name = 'normal_wishart' ;
elseif prior == 4; name = 'indep_normal_wishart';
elseif prior == 5; name = 'SSVS_wishart';
elseif prior == 6; name = 'SSVS_full';
end    
save(sprintf('%s_%s.mat','results',name),'-mat');