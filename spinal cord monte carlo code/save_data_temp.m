i = 2;

%load('continuous_data.mat')

continuous_data(i).wavelength = nm;
%continuous_data(i).model = T;
%continuous_data(i).deposition = F.*muav(T);
continuous_data(i).fluence = F;
%continuous_data(i).Fzx = Fzx;
%continuous_data(i).Fzy = Fzy; 

save('continuous_data.mat',"continuous_data");