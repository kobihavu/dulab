function C = ChiAtNearestStream(FD,S,DEM,chi)

% ChiAtNearestStream computes Chi values at the nearest stream.
%
% Syntax
%
%     C = ChiAtNearestStream(FD,S,DEM,chi)
%
% Description
%
%     ChiAtNearestStream calculates the chi value of each cell in
%     a digital elevation model DEM at the nearest stream cell in S along
%     the flow path in FD.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)
%     chi   node attribute list (nal) of chi values
%
% Output arguments
%
%     C    Chi value at nearest streams (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     chi = chitransform(S,A,'mn',0.45);
%     C = ChiAtNearestStream(FD,S,DEM,chi);
%     imageschs(DEM,C)
%     hc = colorbar;
%     hc.Label.String = '\chi [m] at nearest stream';
%     
% See also: FLOWobj, FLOWobj/flowdistance, FLOWobj/mapfromnalGRIDobj, 
%           STREAMobj, FLOWobj/vertdistance2stream
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016
% Modified by: Kobi Havusha (kobihavu@post.bgu.ac.il)
% Date: 9. September, 2020
% 

narginchk(4,4)

validatealignment(S,DEM);
validatealignment(FD,DEM);
C = DEM; % C is an elevation grid.
C.Z = -inf(DEM.size); % C elevation values are -inf.
C.Z(S.IXgrid) = chi; % C elevation values at the stream network are chi.
ix = FD.ix; % givers
ixc = FD.ixc; % receivers
for r = numel(ix):-1:1 % Goes through all the givers (bottom to top).
    C.Z(ix(r)) = max(C.Z(ix(r)),C.Z(ixc(r))); % C values at the givers index is the maximum value between the giver and receiver (bottom to top).
end