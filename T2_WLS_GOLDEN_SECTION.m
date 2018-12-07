%----------------------------------------------------------------------------------------------------------------------------------
%           MATLAB script implementing the T2 mapping algorithm described in the article:
%----------------------------------------------------------------------------------------------------------------------------------
%
%           FAST AND ACCURATE COMPENSATION OF SIGNAL OFFSET FOR T2 MAPPING
%
%----------------------------------------------------------------------------------------------------------------------------------
%               to be published in "Magnetic Resonance Materials in Physics, Biology and Medicine"
%----------------------------------------------------------------------------------------------------------------------------------
%   Copyright: 	Jan Michálek, Dr.sc.techn.ETH
%               Centre for Biomedical Image Analysis
%               Faculty of Informatics, Masaryk University
%               Botanicka 68a
%               Brno, 602 00
%               Czech Republic
%               jan.michalek@fi.muni.cz
%----------------------------------------------------------------------------------------------------------------------------------

clear variables;close all hidden
warning off
format shortG
slice_directory = input('Enter the directory name of your CPMG slice data :','s');
cd(slice_directory)
first_filename = input('Enter the filename of the first echo:','s');

warning off images:initSize:adjustingMag
warning off MATLAB:MKDIR:DirectoryExists
fullscreen = get(0,'ScreenSize');
screen_left=fullscreen(1);
screen_bottom=fullscreen(2);
screen_width=fullscreen(3);
screen_height=fullscreen(4);
DCM_info = dicominfo(first_filename);
H=DCM_info.Height;% number of rows
W=DCM_info.Width;% number of columns
first_echo_time=20;
last_echo_time = 320;

display_range=[0 400];

T_E=10;
nr_of_echoes=(last_echo_time-first_echo_time)/T_E+1;
t=(first_echo_time:T_E:last_echo_time)';
T_array=repmat(t,[1 H W]);
T_array=permute(T_array,[2 3 1]);
t_ord=1:nr_of_echoes;

Echo_array=zeros(H,W,nr_of_echoes);
layer=1;
for echo_time=t'
    filename=first_filename;
    if echo_time<100
        filename(2:3)=num2str(echo_time);
    else
        filename(1:3)=num2str(echo_time);
    end
    Echo_array(:,:,layer)  = dicomread(filename);
    layer=layer+1;
end

Nonzero_mask_array=(Echo_array>0);

nr_of_iterations=20;

A_i=zeros(H,W,nr_of_iterations);
B_i=zeros(H,W,nr_of_iterations);
C_i=zeros(H,W,nr_of_iterations);
G_i=zeros(H,W,nr_of_iterations);

golden_section=(sqrt(5)-1)/2;
C_min=zeros(H,W);
Min_Echo=min(Echo_array,[],3);
Max_Echo=max(Echo_array,[],3);

layer=10;
C_max=Echo_array(:,:,layer);
C=C_min;
tic

for iteration=1:nr_of_iterations
    C_array=repmat(C,[1 1 size(t)]);
    D_array=Echo_array-C_array;% offset-corrected echo
    F_array=log(double(max(D_array,realmin('double'))));% prevent log from being -Inf, since 0*Inf=NaN
    
    D_array=max(D_array,0);% for weighting, negative offset-corrected echos are replaced with zeros
    Nonzero_mask_array=(D_array>0);
    N=sum(Nonzero_mask_array,3);
    TT_D_array=T_array.*D_array;
    FT_D_array=F_array.*D_array;
    
    D_square=sum(D_array.*D_array,3);
    TT_DD_1=sum(TT_D_array.*D_array,3);
    FT_DD_1=sum(FT_D_array.*D_array,3);
    TT_DD_F=sum(TT_D_array.*FT_D_array,3);
    TT_DD_T=sum(TT_D_array.*TT_D_array,3);
    NUM=TT_DD_1.*FT_DD_1-D_square.*TT_DD_F;
    DEN=D_square.*TT_DD_T-TT_DD_1.*TT_DD_1;
    DEN_plus=DEN.*double(DEN>0)+double(DEN<=0).*realmin('single');%replace denominator with -> 0 wherever denominator==0;
    B=NUM./DEN_plus;%prevents B from being NaN
    B=B.*(DEN~=0)+realmax('single').*(DEN==0);% make T2_GS=0 whenever B infinite
    B=max(B,0);%only non-negative decay rates admissible by the physics
    DEN_B0=D_square;
    DEN_B0_plus=max(DEN_B0,realmin('single'));% realmin will prevent division by zero
    NUM_B0=(FT_DD_1+B.*TT_DD_1);
    B0=NUM_B0./DEN_B0_plus;
    B0=min(B0,realmax('single'));%replace Inf by realmax('single') to prevent 0*Inf=NaN
    B0=max(B0,-realmax('single'));
    B0=B0.*(DEN~=0)-realmax('single').*(DEN==0);% DEN==0 pushes A-> zero
    A=exp(B0);
    A=min(A,realmax('single'));
    A=A.*(DEN~=0)+0.*(DEN==0);
    
    A_array=repmat(A,[1 1 size(t)]);
    B_array=repmat(B,[1 1 size(t)]);
    
    Fit_array=A_array.*exp(-B_array.*T_array)+C_array;
    R_array=Echo_array-Fit_array;% error between the echos and the current fit
    G=sum(R_array.*R_array,3);% cost function values for all pixels
    G=min(G,realmax('single'));
    
    if(iteration==1)%initialize inner GS point
        G_inner=G;
        x_left=C_min;
        x_right=C_max;
        x_inner=x_right-(x_right-x_left)*golden_section;%will lie to the left of x_probe
        A_i(:,:,iteration)=A;
        B_i(:,:,iteration)=B;
        C_i(:,:,iteration)=C;
        C=x_inner;
        A_LL=A;
        B_LL=B;
        C_LL=C;
        T2_LL=1./B;
        T2_LL=min(T2_LL,realmax('single'));
        T2_LL=T2_LL.*(NUM~=0)+realmax('single').*(NUM==0);
        T2_LL=T2_LL.*(DEN~=0)+0.*(DEN==0);
        % imwrite_tiff_single_LZW(single(T2_LL),'T2_LL.tif')
        % imwrite_tiff_single_LZW(single(A_LL),'A_LL.tif')
        % imwrite_tiff_single_LZW(single(B_LL),'B_LL.tif')
        % imwrite_tiff_single_LZW(single(C_LL),'C_LL.tif')
    elseif(iteration==2)% initialize probe GS point
        G_inner=G;
        x_probe=x_left+(x_right-x_left)*golden_section;%will lie to the right of x_inner
        A_i(:,:,iteration)=A;
        B_i(:,:,iteration)=B;
        C_i(:,:,iteration)=C;
        C=x_probe;
        G_max=G;
    else % next GS iteration
        G_probe=G;
        G_probe_LT_G_inner=G_probe<G_inner;
        G_probe_GE_G_inner=~G_probe_LT_G_inner;
        A_i(:,:,iteration)=G_probe_LT_G_inner.*A+G_probe_GE_G_inner.* A_i(:,:,iteration-1);
        B_i(:,:,iteration)=G_probe_LT_G_inner.*B+G_probe_GE_G_inner.* B_i(:,:,iteration-1);
        C_i(:,:,iteration)=G_probe_LT_G_inner.*C+G_probe_GE_G_inner.* C_i(:,:,iteration-1);
        
        x_inner_new=G_probe_LT_G_inner.*x_probe+G_probe_GE_G_inner.*x_inner;
        G_inner=G_probe_LT_G_inner.*G_probe+G_probe_GE_G_inner.*G_inner;
        
        x_probe=G_probe_LT_G_inner.*(x_inner+(x_probe-x_inner)/golden_section)...
            +G_probe_GE_G_inner.*(x_probe+(x_inner-x_probe)/golden_section);
        x_inner=x_inner_new;
        C=x_probe;
    end
    G_i(:,:,iteration)=G_inner;
end %for iteration=1:nr_of_iterations

ElapsedTime=toc;
disp(['ElapsedTime=' num2str(ElapsedTime)...
    ' nr_of_iterations='  num2str(nr_of_iterations)]);

G_last_LT_G_first=G_i(:,:,nr_of_iterations)<G_i(:,:,1);
G_last_GE_G_first=~G_last_LT_G_first;
A_GS=G_last_LT_G_first.*A+G_last_GE_G_first.* A_i(:,:,1);
B_GS=G_last_LT_G_first.*B+G_last_GE_G_first.* B_i(:,:,1);
C_GS=G_last_LT_G_first.*C+G_last_GE_G_first.* C_i(:,:,1);
B_GS=max(B_GS,realmin('single'));
T2_GS=1./B_GS;
T2_GS=min(T2_GS,realmax('single'));

string_nr_of_echoes=[num2str(nr_of_echoes) ' echoes'];

% imwrite_tiff_single_LZW(single(T2_GS),'T2_GS.tif')

Mrange=Max_Echo-Min_Echo;

Extended_Echo_array=zeros(H,W,nr_of_echoes+1);
Extended_Echo_array(:,:,1:nr_of_echoes)=Echo_array(:,:,1:nr_of_echoes);
Extended_Echo_array(:,:,nr_of_echoes+1)=Echo_array(:,:,nr_of_echoes);

Shifted_Echo_array=zeros(H,W,nr_of_echoes+1);
Shifted_Echo_array(:,:,2:nr_of_echoes+1)=Echo_array(:,:,1:nr_of_echoes);
Shifted_Echo_array(:,:,1)=Echo_array(:,:,1);
Difference_Echo_array=Extended_Echo_array-Shifted_Echo_array;
abs_Difference_Echo_array=abs(Difference_Echo_array);
TV_Echo=sum(abs_Difference_Echo_array,3);% variability=100;

Mrange=max(Mrange,realmin('double'));
TV_Mrange_ratio=TV_Echo./Mrange;
variability=2;
% variability=1.4;
Noise_mask=double(TV_Mrange_ratio<variability);% Noise_mask == 1 means acknowledged pixels
% imwrite_tiff_single_LZW(single(Noise_mask),'Noise_mask.tif')

T2_GS_clipped=min((T2_GS),display_range(2));
T2_GS_denoised=T2_GS_clipped.*Noise_mask;
image_nr=0;
image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','Noise_mask',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
imshow(Noise_mask,'InitialMagnification','fit','Border','tight')

image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','T2_GS_clipped',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
imshow(T2_GS_clipped,display_range,'InitialMagnification','fit','Border','tight')
% imwrite(uint8(T2_GS_clipped),'T2_GS_clipped_imwrite.tif')
% export_fig T2_GS_clipped.tif -transparent -native

image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','T2_GS_denoised',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
imshow(T2_GS_denoised,display_range,'InitialMagnification','fit','Border','tight')
% export_fig T2_GS_denoised.tif -transparent -native

image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','T2_GS_colormap' ,...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
imshow(T2_GS_denoised,display_range,'InitialMagnification','fit')
colormap(gca,hot)
colorbar
box on
% export_fig T2_GScolor.tif -transparent -native

image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','T2_LL_colormap',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
imshow(T2_LL.*Noise_mask,display_range,'InitialMagnification','fit')
colormap(gca,hot)
colorbar
box on
% export_fig T2_LLcolor.tif -native

image_nr=image_nr+1;
figure(image_nr);
set(image_nr,'Name','T2_GS_red_denoised' ,...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
Noise_mask_red = cat(3,~Noise_mask, 0*Noise_mask, 0*Noise_mask);
T2_GS_gray=mat2gray(T2_GS_denoised);
T2_GS_rgb=cat(3,T2_GS_gray,T2_GS_gray,T2_GS_gray);
T2_GS_red_denoised=T2_GS_rgb.*~Noise_mask_red+Noise_mask_red;
imshow(T2_GS_red_denoised,'InitialMagnification','fit','Border','tight')
% export_fig T2_GS_red_denoised.tif  -transparent -native
commandwindow

