function [in,inci,pd] = intensity(P,alpha)

%INTENSITY Calculate intensity (density) of points on the stream network
%
% Syntax
%
%     in = intensity(P)
%     [in,inci,pd] = intensity(P,alpha);
%
% Description
%
%     intensity measures the expected number of points per unit length, 
%     which is referred to as the intensity of a point process.
% 
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%     alpha  significance level (default = 0.05)
%
% Output arguments
%
%     in     intensity (scalar)
%     inci   confidence bounds of intensity (two element vector)
%     pd     instance of prob.GammaDistribution. Contains the probability
%            distribution of the intensity (see Furbish 2021)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     [int,inci,pd] = intensity(P);
%
%     % Plot intensity estimate
%     ci = icdf(pd,[.0001 .9999]);
%     inti = linspace(ci(1),ci(2),500);
%     p  = pdf(pd,inti);
%     plot(inti,p)
%     xlabel('\lambda [m^{-1}]')
%     ylabel('f_{\lambda} [m]')
%     hold on
%     xline(inci,'k--')
%     xline(int,'k')
%     legend('PDF','lower CI 95%','upper CI 95%','mean') 
%
% See also: PPS, PPS/npoints 
%
% Reference: Furbish (2021): The “uncertainty principle” of a Poisson point 
%            process.
% https://my.vanderbilt.edu/davidjonfurbish/files/2013/06/Uncertainty-Principle.pdf
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 4. March, 2026

arguments
    P PPS
    alpha (1,1) {mustBePositive} = 0.05
end

n  = npoints(P);
s  = tlength(P);
in = n/s;   
if nargout > 1
    % If error bounds required
    if nargin == 1
        alpha = 0.05;
    else
        validateattributes(alpha,{'numeric'},{'>',0,'<',1},...
            'intensity','alpha',2);
    end
    
    pd = makedist('Gamma','a',n,'b',1/s);
    % int = mean(pd);
    inci  = icdf(pd,[alpha/2 1-alpha/2]);

%     % Alternative way to derive confidence intervals of the intensity
%     % Fit loglinear model with intercept only
%     mm = fitloglinear(P,getnal(P.S)+1,'intercept',false);
%     [p,pci]   = predict(mm,1,'alpha',alpha);
%     d   = l./numel(P.S.x);
%     in = p./d; %.cellsize;   
%     inci = pci./d;
end