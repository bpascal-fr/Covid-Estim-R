% Compute the global infectiousness time series independently for each
% country monitored and during the entire time period:
%
%    Zphi(t) = sum_s Phi(s) * Z(t-s) for s = 1, ..., tau_phi,
%
% where Phi is the serial interval function, modeling the distribution of
% the random delay between primary and secondary infection discretized at
% the week scale.
%
% Border effects at t = 1 are handled by cropping and normalizing the
% serial interval function as described in
% - Du, J., Pascal, B., & Abry, P. (2023, August). Compared performance of
% Covid19 reproduction number estimators based on realistic synthetic data.
% GRETSI’23 XXIXème Colloque Francophone De Traitement Du Signal Et Des Images.
%
% Default serial interval function is chosen as the serial interval
% function of Covid19 estimated in
% - D. Cereda, M. Tirani, F. Rovida, V. Demicheli, M. Ajelli, P. Poletti,
% F. Trentini, G. Guzzetta, V. Marziano, A. Barone et al. (2020). The early
% phase of the COVID-19 outbreak in Lombardy, Italy.
% Preprint arXiv:2003.09320
% - F. Riccardo, M. Ajelli, X. D. Andrianou, A. Bella, M. Del Manso,
% M. Fabiani, S. Bellino, S. Boros, A. M. Urdiales, V. Marziano et al. (2020).
% Epidemiological characteristics of COVID-19 cases and estimates of the
% reproductive numbers 1 month into the epidemic, Italy, 28 January
% to 31 March 2020. Euro Surveillance
% discretized at a week scale, and not at a day scale as in usual studies.
%
%
% Implementation B. Pascal and P. Abry
% May 2024


function [Zphi, Z, Phi] = Phi_Week_normal(Z,Phi)

    % Inputs:  - Z: weekly new infection counts of size C x W
    %          - Phi: weekly discretized serial interval function (optional)
    %            (default Gamma distribution of mean 6.6 days and standard
    %            deviation 3.5 days cropped at 4 weeks and normalized
    %            modeling the serial interval function of Covid19)
    %
    % Outputs: - Zphi: global infectiousness (not defined on first week) hence of size C x W-1
    %          - Z: new infection counts with first week cropped, hence of size C x W-1
    %          - Phi: weekly discretized serial interval function
    %


    % Handle possibly univariate Z of size T x 1
    [d1,d2] = size(Z);
    if d2 == 1
        Z = reshape(Z,1,d1);
    end

    if nargin <=1

        % Serial interval function modeled by a Gamma distribution with
        % - mean: 6.6 days
        % - standard deviation 3.5 days
        % truncated at 25 days and normalized.

        shape    = 1/0.28;
        scale    = 1.87;
        tau_day  = 25;
        Phi_Day  = gampdf(0:tau_day,shape,scale);
        Phi_Day  = Phi_Day/sum(Phi_Day); % normalize the weights applied to past tau_phi infection counts

        % tau - 1 should be a multiple of seven
        if mod(tau_day,7)
            add_day  = 7 - mod(tau_day,7);
            tau_day  = tau_day + add_day;
            Phi_Day  = [Phi_Day, zeros(1,add_day)];
        end
        
        % Aggregation over seven days to get a weekly discretized interval function
        Phi          = zeros(1,tau_day/7);
        for w = 1:tau_day/7
            tmp_Phi  = Phi_Day(7*(w-1)+2:7*w+1);
            Phi(w+1) = sum(tmp_Phi);
        end
        tau_week     = length(Phi)-1;

    else

        if ~(Phi(1) == 0)
            tau_week = length(opts.Phi);
            Phi      = reshape(opts.Phi,1,tau_week);
            Phi      = [0, Phi];
        else
            tau_week = length(opts.Phi)-1;
            Phi      = reshape(opts.Phi,1,tau_week+1);
        end

    end



    % compute the convolution of Z by Phi, territory per territory if multivariate Z
    Zphi             = convn(Z,Phi);
    Zphi             = Zphi(:,1:size(Z,2));  save('Zphi2','Zphi')

    % recompute the first tau_phi elements with normalized Phi
    Phi              = reshape(Phi,length(Phi),1);
    Zphi(1:tau_week) = 0;
    for t = 2:tau_week

        Phi_crop    = Phi(1:t)./sum(Phi(1:t));
        fZ          = fliplr(Z(:,1:t));
        Zphi(:,t)   = fZ*Phi_crop;

    end

    % Crop the first day at which Zphi = 0
    Z               = Z(:,2:end);
    Zphi            = Zphi(:,2:end);

    % Reshape serial interval function for output
    Phi              = reshape(Phi,1,length(Phi));

    % Return Z, Zphi of size T-1 x 1 if input Z is univariate of size T x 1
    if d2 == 1
        Z    = reshape(Z,d1-1,d2);
        Zphi = reshape(Zphi,d1-1,d2);
    end
end