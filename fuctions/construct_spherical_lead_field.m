function [LFM_sphere] = construct_spherical_lead_field(EEG,r0,r1,r2,r3,sig1,sig2,sig3,DipN)
% This functions creates a spherical lead-field that can be used in SOUND
% to correct the data. 
% 
% Input variable:
%
% EEG = The studied EEGLAB dataset. The function expects that EEG contains
% the chanlocs on a spherical surface.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

% r0 = distance of the current dipoles from the origin in [meters]
% (default, if input r0=[], r0=76*1e-3, as in the original article)
if nargin < 2 || isempty(r0)
   r0 = 76*1e-3;
end

% r1 = Radius of the inner skull surface in [meters] (default, if input
% r1=[], r1=81*1e-3, as in the original article)
if nargin < 3 || isempty(r1) 
   r1 = 81*1e-3;
end

% r2 = Radius of the outer skull surface in [meters] (default, if input
% r2=[], r2=85*1e-3, as in the original article)
if nargin < 4 || isempty(r2) 
   r2 = 85*1e-3;
end

% r3 = Radius of the scalp surface in [meters] (default, if input
% r3=[], r3=88*1e-3, as in the original article)
if nargin < 5 || isempty(r3) 
   r3 = 88*1e-3;
end

% sig1 = Conductivity of th brain in [1/0mega*1/m] (default, if input
% sig1=[], sig1=0.33, as in the original article)
if nargin < 6 || isempty(sig1) 
   sig1 = 0.33;
end

% sig2 = Conductivity of the skull in [1/0mega*1/m] (default, if input
% sig2=[], sig2=0.33/50, as in the original article)
if nargin < 7 || isempty(sig2) 
   sig2 = 0.33/50;
end

% sig3 = Conductivity of th scalp in [1/0mega*1/m] (default, if input
% sig3=[], sig3=0.33, as in the original article)
if nargin < 8 || isempty(sig3) 
   sig3 = 0.33;
end

% DipN = The number of the brain noise source dipoles (default, in input 
% DipN=[], DipN=5000, as in the original article)
if nargin < 9 || isempty(DipN)
   DipN = 5000;
end


% Reading the electrode locations from the chanlocs field in the given EEG
% file:

k = 1;
for i=1:length(EEG.chanlocs)
%     if strcmp(EEG.chanlocs(i).type,'EEG')
    elec_coords(k,1) = EEG.chanlocs(i).X;
    elec_coords(k,2) = EEG.chanlocs(i).Y;
    elec_coords(k,3) = EEG.chanlocs(i).Z;
    k = k+1;
%     end
end

% Projecting the coordinates on a unit sphere:
elec_coords = elec_coords./repmat(sqrt(sum(elec_coords.^2,2)),[1,3]);

%%%% Computing the spherical lead field using function leadfield1

% THE SPHERICAL HEAD MODEL
rad=[r1,r2,r3];
sig=[sig1,sig2,sig3];

 % the radius of the cortex hemisphere on which the dipoles lie
R=r3*elec_coords; % electrode positions,

P=randn(3,DipN);
P=r0*P./(ones(3,1)*sqrt(sum(P.^2))); 
Pns=[P(1,:);P(2,:);(P(3,:))]; % dipole positions on the upper hemisphere
                               % with radius r0,
%Moments of the unit dipoles                               
Q = Pns;                              
Qns = Q./repmat(sqrt(sum(Q.^2,1)),[3,1]);

%Computing the leadfield accoring to [1]
LFM_sphere = leadfield1(R',Pns,Qns,rad,sig,30);

%[1] Mutanen, Tuomas P., et al. "Recovering TMS-evoked EEG responses masked by muscle artifacts." Neuroimage 139 (2016): 157-166.

end
