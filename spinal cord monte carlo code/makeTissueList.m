function tissue = makeTissueList(nm)

%   Returns the tissue optical properties at the wavelength nm:


%% Load spectral library
load spectralLIB.mat
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muamel,nm);
MU(:,5) = interp1(nmLIB,muafat,nm);
LOADED = 1;

%% Create tissueList

j=1;
tissue(j).name  = 'air';
tissue(j).mua   = 0.0001; % Negligible absorption yet still tracks, 
tissue(j).mus   = 1.0;    % but take steps in air
tissue(j).g     = 1.0;    % and don't scatter.

j=2;
tissue(j).name  = 'CSF';
%MU(3)
tissue(j).mua   = MU(3); %water absorbance coefficient
tissue(j).mus   = 0.12; %0.12
tissue(j).g =     0.16; %0.16

j=3;
tissue(j).name  = 'Dura Mater';
B       = 0.002;
S       = 0.68;
W       = 0.65;
F       = 0;
M       = 0;
musp500 = 43.6;
fray    = 0.41;
bmie    = 0.562;
gg      = 0.8; %should be 0.8

musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M F]';
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j = 4;
tissue(j).name = 'White Matter';
B = 0.027; 
S = 0.68;
W = 0.70;
M = 0;
F = 0.19;
musp500 = 27.4;
fray    = 0.315;
bmie    = 1.087;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M F]';
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=5;
tissue(j).name  = 'Gray Matter';
B = 0.052;
S = 0.68;
W = 0.85;
M = 0;
F = 0.06;
musp500 = 27.4;
fray    = 0.315;
bmie    = 1.087;
gg      = 0.87;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M F]';
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

%We need to add one of these for our hydrogel
%Lets assume it would look similar to how
%Air and CSF is defined
%WHAT WE SHOULD ACCOUNT FOR
%Our gel is made of a few components
%Photoinitializer (RU SPS)
%Gelma components
%visually I expect very low scattering from this
%As when I look at the thing it doesn't look blurry
%Idea: Model as water with the only thing that causes absorbance being
%RU SPS
%Ru is an atom, SPS is a tiny molecule
%This leads to primarily rayleigh scattering

j=6;
tissue(j).name = 'Hydrogel';

musp500 = 43.6; %???
gg      = 0.0; %almost emtirely rayleigh
fray    = 1; %all rayleigh

musp = musp500*(fray*(nm/500).^-4); %scattering

exCoeff = 14600; %M^-1 cm^-1, at 450 nm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8201352/
concRu = 1.0/1000; %1 mM%typical conc of RU. Only Ru absorbs in this scenario
concSPS = 5/1000; % 5 mM%SPS conc

%https://omlc.org/classroom/ece532/class3/bilirubin.html
%This most definitely only works for wavelength = 450 nm
tissue(j).mua   = exCoeff*concRu*log(10); % blame chemists for the ln(10)
tissue(j).mus   = 0.12; %this is the stuff for CSF
tissue(j).g =     0.16; %We may want to make this better at some point



disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

