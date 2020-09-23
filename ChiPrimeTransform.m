function c = ChiPrimeTransform(S,A,varargin)

%ChiPrimeTransform Coordinate transformation using the integral approach
%
% Syntax
%
%     c = ChiPrimeTransform(S,A)
%     c = ChiPrimeTransform(S,A,pn,pv,...)
%
% Description
%
%     ChiPrimeTransform transforms the horizontal spatial coordinates of a river
%     longitudinal profile using an integration in upstream direction of
%     drainage area (chi, see Perron and Royden, 2013).
%
% Input arguments
%
%     S     STREAMobj
%     A     upslope area as returned by the function flowacc
%     
% Parameter name/value pairs
%
%     'a0'     reference area (default=1e6)
%     'mn'     mn-ratio (default=0.45)
%     'plot'   0 : no plot (default)
%              1 : chimap
%     'UoverK' GRIDobj or node attribute list of spatial pattern of uplift (U) divided by erodibility (K) (default=[])
%
% Output argument
%
%     c     node attribute list (nal) of chi values
%
% Example: chi computation considering differential uplift, uniform K.
% 
%
%     DEM = GRIDobj('Diff_EXP_17hr.tif');
%     DEM.Z(DEM.Z<-9998)=NaN;
%     FD = FLOWobj(DEM,'preprocess','carve');
%     ST = STREAMobj(FD,'minarea',5000);
%     A = flowacc(FD);
%     Ufast = 0.021; % Fast uplift rate is 0.021 m/hr.
%     Umid = 0.016; % Medium uplift rate is 0.016 m/hr.
%     Uslow = 0.008; % Slow uplift rate is 0.008 m/hr.
%     U_K = DEM;
%     U_K.Z = NaN(DEM.size); 
%     columns = U_K.size(2); 
%     Third = floor(columns/3); % represents 1/3 of the width of the box.
%     U_K.Z(:,1:Third) = Ufast;
%     for i=Third+1:2*Third
%       U_K.Z(:,i) = ((i-Third)/(Third))*(Umid-Ufast)+Ufast;
%     end
%     for i=2*Third+1:columns
%       U_K.Z(:,i) = ((i-2*Third)/(Third))*(Uslow-Umid)+Umid;
%     end
%     chi = ChiPrimeTransform(ST,A,'mn',0.45,'UoverK', U_K);
%     plotc(ST,chi)    
%
% See also: chiplot
%
% References:
%     Willett, Sean D., et al. (2014): Dynamic reorganization of river basins. Science 343.6175.?
%     Perron, J. & Royden, L. (2013): An integral approach to bedrock river 
%     profile analysis. Earth Surface Processes and Landforms, 38, 570-576.
%     [DOI: 10.1002/esp.3302]

% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. December, 2015
% Modified by: Kobi Havusha & Liran Goren 
% Date: 9. September, 2020

% Parse Inputs
p = inputParser;         
p.FunctionName = 'ChiPrimeTransform';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj') || isnal(S,x));

addParamValue(p,'mn',0.45,@(x) isscalar(x) || isempty(x));
addParamValue(p,'a0',1e6,@(x) isscalar(x) && isnumeric(x));
addParamValue(p,'plot',false);
addParamValue(p,'correctcellsize',true,@(x) isscalar(x));
addParamValue(p,'UoverK',[],@(x) isempty(x) || isnal(S,x) || isa(x,'GRIDobj'));


parse(p,S,A,varargin{:});

% get node attribute list with elevation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = getnal(S,A);
elseif isnal(S,A)
    a = A;
else
    error('Imcompatible format of second input argument')
end

if ~isempty(p.Results.UoverK)
    calcwithk = true;
    if isnal(S,p.Results.UoverK)
        UoverK = p.Results.UoverK;
    else 
        UoverK = getnal(S,p.Results.UoverK);
    end
else
    calcwithk = false;
end
        
       
if p.Results.correctcellsize
    a = a.*S.cellsize^2;
end

if ~calcwithk
    a = ((p.Results.a0) ./a).^p.Results.mn;
else
    % This transformation is only possible if we assume that n in the
    % mn-ratio is one.
    a = UoverK.*((p.Results.a0)./a).^p.Results.mn;
end
c = cumtrapz(S,a);

if p.Results.plot
    plotc(S,c)
end

