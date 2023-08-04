clear all
close all

% imaging parameters
PSF_pixel_size = 53.5;
output_pixel_size = 53.5;
N = 256;
line_width = 3;
file = 'C:\Users\ew535\Documents\Example_Image_Bactera_Liposome\488_561_50uM_fL6-fairSIM-1.tif';


x = [90 170]; x = x/2;
y = [126 126]; y = y/2;


bacteria_data = double(imread(file,1));
bacteria_data = make_square(bacteria_data,N/2);
bacteria_data(bacteria_data<0)=0;
bacteria_data = bacteria_data./max(bacteria_data(:));
f = figure;
f.Position = [200 200 1400 700];
subplot(2,4,1)
imagesc(bacteria_data); axis off; axis square; hold on;
plot([x(1),x(2)],[y(1),y(2)],'Color','r','LineWidth',2);
plot([x(1),x(2)],[y(1)+line_width,y(2)+line_width],'Color','r','LineWidth',2);
drawnow
mask = imbinarize(bacteria_data);
mask = bwareafilt(mask,[200 600]);
stats_data = regionprops('table',mask,'Centroid','MajorAxisLength','MinorAxisLength');

membrane_data = double(imread(file,2));
membrane_data = make_square(membrane_data,N/2);
membrane_data(membrane_data<0)=0;
membrane_data = membrane_data./max(membrane_data(:));

subplot(2,4,2)
imagesc(membrane_data); axis off; axis square; hold on;
plot([x(1),x(2)],[y(1),y(2)],'Color','r','LineWidth',2);
plot([x(1),x(2)],[y(1)+line_width,y(2)+line_width],'Color','r','LineWidth',2);
drawnow;

% bacteria parameters
heights = stats_data.MajorAxisLength*output_pixel_size;
width = stats_data.MinorAxisLength*output_pixel_size;
m_thickness = PSF_pixel_size;

% convert to pixels

heights = (heights-width)/PSF_pixel_size;
m_thickness = m_thickness/PSF_pixel_size;

radius = width/(2*PSF_pixel_size);

% begin calculations

disp('Calculating shapes...')

% bacteria membrane
bacteria_pill = getPill(N-1,radius, radius, radius, heights);
membrane_pill = getPill(N-1,radius+m_thickness, radius+m_thickness, radius+m_thickness, (heights+(m_thickness/2)));

membrane_pill = membrane_pill-bacteria_pill;
membrane_pill(membrane_pill<0)=0;

saveResult(membrane_pill,'membranee',1)
saveResult(bacteria_pill,'bacteria',1)

%% calculate imaging parameters
% The following code generates model SIM images for sample under 647nm illumination
% Change the wavelength used in calculating the PSF when different fluorophores are used

PSFem = zeros(N,N,N,'single'); % pre-allocate for speed
for page = 1:N
    textwaitbar(page, N, 'Loading detection PSF')
    PSFem(:,:,page) = single(imread('647 PSF.tif',page));
end
disp('Calculating detection OTF...')

% Just for fun let's find the best fft algorithm
fftw('swisdom',[]);
fftw('planner','patient');
OTFem = fftshift(fftn(PSFem));

% Get SIM OTF from WF OTF
[xsize,ysize,pages]=size(OTFem);
[Y,X]=meshgrid(1:ysize,1:xsize);

xc=floor(xsize/2+1);% the x-coordinate of the center
yc=floor(ysize/2+1);% the y-coordinate of the center
yr=Y-yc;
xr=X-xc;
R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)
n_filt = 1 - exp(-0.5*R.^1.2); % notch filter used by FairSIM for stripe suppression and optical sectioning
% Ensure the width of the notch filter matches the "OTF supression" parameter used in FairSIM

figure(); imagesc(n_filt);
OTFem = OTFem.*n_filt;
OTF_sim= zeros(size(OTFem));

% generate the SIM OTF from the kx-ky SIM pattern parameters reported in FairSIM
for i = 1:6
    shift_x = round(30.8 * sin(0.5+(2*pi/6)*(i-1))); 
    shift_y = round(30.8 * cos(0.5+(2*pi/6)*(i-1)));
    OTF_sim = OTF_sim + circshift(OTFem,[shift_x shift_y 0]);
end

% figure(); imagesc(log(abs(OTF_sim(:,:,128))));
OTF_sim = fftshift(OTF_sim);
clear PSFem OTFem n_filt R yr xr xc yc xsize ysize X Y shift_x shift_y% clean the workspace

% Calculate membrane profile
disp('Calculating sample spectra...')
mem_temp = (fftn(membrane_pill));
mem_temp = mem_temp.*OTF_sim;
disp('Calculating image response...')
mem_temp =  real(fftshift(fftn(mem_temp)));
mem_temp(mem_temp<0)=0;
mem_temp = mem_temp./max(mem_temp(:));
saveResult(mem_temp,'membrane data',1)
mem_temp = reshape(mem_temp, 2, N/2, 2, N/2, 2, N/2);
mem_temp = squeeze(sum(mem_temp, [1,3,5])) / 8;

membrane_pill = reshape(membrane_pill, 2, N/2, 2, N/2, 2, N/2);
membrane_pill = squeeze(sum(membrane_pill, [1,3,5])) / 8;

% saveResult(mem_temp,'membrane data',1)
% saveResult(mem_temp(:,:,pages/2),'membrane section',1);

% Calculate bacteria profile
disp('Calculating sample spectra...')
bacteria_temp = (fftn(bacteria_pill));
bacteria_temp = bacteria_temp.*OTF_sim;
disp('Calculating image response...')
bacteria_temp =  real(fftshift(fftn(bacteria_temp)));
bacteria_temp(bacteria_temp<0)=0;
bacteria_temp = bacteria_temp./max(bacteria_temp(:));
saveResult(bacteria_temp,'bacteria data',1)
bacteria_temp = reshape(bacteria_temp, 2, N/2, 2, N/2, 2, N/2);
bacteria_temp = squeeze(sum(bacteria_temp, [1,3,5])) / 8;

bacteria_pill = reshape(bacteria_pill, 2, N/2, 2, N/2, 2, N/2);
bacteria_pill = squeeze(sum(bacteria_pill, [1,3,5])) / 8;

%saveResult(bacteria_temp,'bacteria data',1)
%saveResult(bacteria_temp(:,:,pages/2),'bacteria section',1);

% mask bacteria

mask = imbinarize(bacteria_temp(:,:,N/4));
mask = bwareafilt(mask,[200 600]);
stats = regionprops('table',mask,'Centroid','MajorAxisLength','MinorAxisLength');


% show bacteria 

subplot(2,4,3)
imagesc(bacteria_temp(:,:,N/4)); axis square; hold on
plot([x(1),x(2)],[y(1)+line_width,y(2)+line_width],'Color','r','LineWidth',2);
plot([x(1),x(2)],[y(1),y(2)],'Color','r','LineWidth',2);

% show membrane 
subplot(2,4,4);
imagesc(mem_temp(:,:,N/4)); axis square; hold on
plot([x(1),x(2)],[y(1)+line_width,y(2)+line_width],'Color','r','LineWidth',2);
plot([x(1),x(2)],[y(1),y(2)],'Color','r','LineWidth',2);

%% Unused but the following code automatically generates line profiles across bacteria

bacteria_data_profile = improfile(bacteria_data,x,y);
for i = 1:line_width
    bacteria_data_profile = bacteria_data_profile + improfile(bacteria_data,x+i,y);
end
bacteria_data_profile = bacteria_data_profile./line_width;
band_x = linspace(0,length(bacteria_data_profile),length(bacteria_data_profile));
ft = fittype( 'a*(1/(1+exp(-b*(x-c))))*(1-(1/(1+exp(-b*(x-d)))))',...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.957166948242946 1.51 10 40];
% Fit model to data.
sigmoid_data = fit( band_x.', bacteria_data_profile, ft, opts );
bacteria_data_width = abs(sigmoid_data.d - sigmoid_data.c);

subplot(2,4,5);
plot(sigmoid_data,band_x,bacteria_data_profile); hold on



% plot membrane_data profile
mem_data_profile = improfile(membrane_data,x,y);
for i = 1:line_width
    mem_data_profile = mem_data_profile + improfile(membrane_data,x+i,y);
end
mem_data_profile = mem_data_profile./line_width;

%%  plot membrane data profile

% plot membrane profile
subplot(2,4,6);
mem_profile = improfile(mem_data_profile,x,y);
mid_point = round(sigmoid_data.c+abs(sigmoid_data.d-sigmoid_data.c)/2);
half_width = [1:mid_point];
plot(mem_data_profile(1:mid_point),'b*','DisplayName','left data');hold on
plot((half_width+mid_point-1),mem_data_profile(mid_point:2*mid_point-1),'b*','DisplayName','right data');

% compute left side
f_mem_left = fit(half_width.',mem_data_profile(1:mid_point),'smoothingspline');
xFit = linspace(1,mid_point,400);
yFit = f_mem_left(xFit);
plot(xFit,yFit,'DisplayName','left fit')

left_peak_max = xFit(yFit==max(yFit));
left_membrane_sep = sigmoid_data.c - left_peak_max;

% compute right side
half_width = [1:mid_point];
f_mem_right = fit(half_width.',mem_data_profile(mid_point:2*mid_point-1),'smoothingspline');
xFit = linspace(1,mid_point,400);
yFit = f_mem_right(xFit);
plot(xFit+mid_point-1,yFit)

right_peak_max = xFit(yFit==max(yFit))+mid_point-1;
right_membrane_sep = right_peak_max - sigmoid_data.d;

mean_membrane_sep = (right_membrane_sep+left_membrane_sep)/2;

%%  plot bacteria profile
bacteria_profile = improfile(bacteria_temp(:,:,N/4),x,y);
band_x = linspace(0,length(bacteria_profile),length(bacteria_profile));
ft = fittype( 'a*(1/(1+exp(-b*(x-c))))*(1-(1/(1+exp(-b*(x-d)))))',...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.957166948242946 1.51 10 40];
% Fit model to data.
sigmoid = fit( band_x.', bacteria_profile, ft, opts );
subplot(2,4,7);
plot(sigmoid,band_x,bacteria_profile); hold on
plot(improfile(bacteria_pill(:,:,N/4),x,y))
bacteria_width = abs(sigmoid.d - sigmoid.c);

% plot membrane profile
subplot(2,4,8);
mem_profile = improfile(mem_temp(:,:,N/4),x,y);
mid_point = round(sigmoid.c+abs(sigmoid.d-sigmoid.c)/2);
half_width = [1:mid_point];
plot(mem_profile(1:mid_point),'b*','DisplayName','left data');hold on
plot((half_width+mid_point-1),mem_profile(mid_point:end),'b*','DisplayName','right data');

% compute left side
f_mem_left = fit(half_width.',mem_profile(1:mid_point),'smoothingspline');
xFit = linspace(1,mid_point,400);
yFit = f_mem_left(xFit);
plot(xFit,yFit,'DisplayName','left fit')

% compute right side
f_mem_right = fit(half_width.',mem_profile(mid_point:end),'smoothingspline');
xFit = linspace(1,mid_point,400);
yFit = f_mem_right(xFit);
plot(xFit+mid_point-1,yFit)

%% show results table
figure()
bacteria_number = {'Bacteria 1';'average'};
width = [bacteria_width;bacteria_width];
x_shift = [71;69];
Weight = [176;163];
T = table(width,x_shift,Weight,'RowNames',bacteria_number);
% Get the table in string form.
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

%% Helper functions

function output = make_square(image,target_size)
    if length(size(image)) == 2 
        [X, Y] = size(image);
        if X > target_size || Y > target_size
            disp('Incompatible image size for padding')
            return
        end
        x_pad = round((target_size-X)/2);
        y_pad = round((target_size-Y)/2);
        temp = padarray(image,[x_pad,y_pad],min(image(:)),'both');
        [X, ~] = size(temp);

        if X ~= target_size 
           PSF = fspecial('gaussian',20,5);
           temp = edgetaper(temp,PSF);
           temp = temp(1:target_size,1:end);
           [x, y] = size(temp);
           [x_co, y_co] = meshgrid(1:y,1:x);
           H=fftshift(fft2(temp)); %// Compute 2D Fourier Transform
           x0=0.5; %// Define shifts
           y0=0;
           %// Define shift in frequency domain
           xF = (x_co-(x/2))/x; yF = (y_co-(x/2))/y;    
           %// Perform the shift
           H=H.*exp(-1i*2*pi.*(xF*x0+yF*y0)/200);   
           %// Find the inverse Fourier Transform
           temp=real(ifft2(ifftshift(H)));
        end
        [~, Y] = size(temp);
        if Y ~= target_size 
           PSF = fspecial('gaussian',20,5);
           temp = edgetaper(temp,PSF);
           temp = temp(1:end,1:target_size);
           [x, y] = size(temp);
           [x_co, y_co] = meshgrid(1:y,1:x);
           H=fftshift(fft2(temp)); %// Compute 2D Fourier Transform
           x0=0.5; %// Define shifts
           y0=0;
           %// Define shift in frequency domain
           xF = (x_co-(x/2))/x; yF = (y_co-(x/2))/y; 
           %// Perform the shift
           H=H.*exp(-1i*2*pi.*(xF*x0+yF*y0)/200); 
           %// Find the inverse Fourier Transform
           temp=real(ifft2(ifftshift(H)));
        end

        output = temp;

    elseif length(size(image)) == 3
       [X, Y, pages] = size(image);
       if X > target_size || Y > target_size
            disp('Incompatible image size for padding')
            return
       end
       x_pad = round((target_size-X)/2);
       y_pad = round((target_size-Y)/2);
       output = zeros([target_size,target_size,pages]);
       for p = 1:pages
           temp = padarray(image(:,:,pages),[x_pad,y_pad],min(image(:)),'both');
           [X, Y] = size(temp);
           if X > target_size || Y > target_size
                disp('Incompatible image size for padding')
                return
           end
           
           [X, ~] = size(temp);
           if X ~= target_size 
               PSF = fspecial('gaussian',20,5);
               temp = edgetaper(temp,PSF);
               temp = temp(1:target_size,1:end);
               [x, y] = size(temp);
               [x_co, y_co] = meshgrid(1:y,1:x);
               H=fftshift(fft2(temp)); %// Compute 2D Fourier Transform
               x0=0.5; %// Define shifts
               y0=0;
               %// Define shift in frequency domain
               xF = (x_co-(x/2))/x; yF = (y_co-(x/2))/y;    
               %// Perform the shift
               H=H.*exp(-1i*2*pi.*(xF*x0+yF*y0)/200);   
               %// Find the inverse Fourier Transform
               image=real(ifft2(ifftshift(H)));
            end

            [~, Y] = size(temp);
            if Y ~= target_size 
               PSF = fspecial('gaussian',20,5);
               temp = edgetaper(temp,PSF);
               temp = temp(1:end,1:target_size);
               [x, y] = size(temp);
               [x_co, y_co] = meshgrid(1:y,1:x);
               H=fftshift(fft2(temp)); %// Compute 2D Fourier Transform
               x0=0.5; %// Define shifts
               y0=0;
               %// Define shift in frequency domain
               xF = (x_co-(x/2))/x; yF = (y_co-(x/2))/y; 
               %// Perform the shift
               H=H.*exp(-1i*2*pi.*(xF*x0+yF*y0)/200); 
               %// Find the inverse Fourier Transform
               image=real(ifft2(ifftshift(H)));
            end

        output(:,:,p) = temp;
       end
    else
        disp('Incompatible image size for padding')
    end
end
function v = getPill(N,xr,yr,zr,length)
    
    x = -N/2:1:N/2;
    y = -N/2:1:N/2;
    z = -N/2:1:N/2;
    [X,Y,Z] = meshgrid(x,y,z);

    v = zeros(size(X));

    elipse = ((X.^2)/(xr^2))+(((Y+length/2).^2)/(yr^2))+((Z.^2)/(zr^2));
    v(elipse<=1)=1;

    elipse = ((X.^2)/(xr^2))+(((Y-length/2).^2)/(yr^2))+((Z.^2)/(zr^2));
    v(elipse<=1)=1;

    [~,rho,z] = cart2pol(X,Y,Z);

    compoundCondInd = (rho<xr) & (z>(-length/2)) & (z<length/2);
    compoundCondInd = permute(compoundCondInd, [3 2 1]);
    v(compoundCondInd) =1;

end

function output = FT_upsample(image)
    PSF = fspecial('gaussian',20,5);
    image = edgetaper(image,PSF);
    image_FT = fftshift(fft2(image));
    image_FT = padarray(image_FT,2*size(image_FT),0,'both');
    image = real((ifft2(fftshift(image_FT))));
    image(image<0)=0;
    output = image./max(image(:));
end

function [] = saveResult(x,name,progress)
    x = x - min(x(:));
    x = 65535 * x/max(x(:));
    x = uint16(x);
    path = strcat(name,'.tif');
    if isfile(path)
     % File exists.
    delete(path);
    end
    
    msg = strcat('Saving image: ',path);
    frames = size(x,3);
    if frames >50
        writeFast(x,path,progress)
    else
        for p =1:frames
            if progress == 1
            textwaitbar(p, frames, msg);
            end
            imwrite(x(:,:,p),path,'writemode','append');
        end
    end
end
function [] = writeFast(x,path,progress)
    msg = strcat('Saving image: ',path);
    fTIF = Fast_Tiff_Write(path,0.125,0);
	for page =1:size(x,3)
        if progress == 1
        textwaitbar(page, size(x,3), msg);
        end
        fTIF.WriteIMG(x(:,:,page));
    end
    fTIF.close;
end
function textwaitbar(i, n, msg)
% A command line version of waitbar.
% Usage:
%   textwaitbar(i, n, msg)
% Input:
%   i   :   i-th iteration.
%   n   :   total iterations.
%   msg :   text message to print.
%
% Date      : 05/23/2019
% Author    : Xiaoxuan He   <hexxx937@umn.edu>
% Institute : University of Minnesota
%
% Previous percentage number.
persistent i_prev_prct;
% Current percentage number.
i_prct = floor(i ./ n * 100);
% Print message when counting starts.
if isempty(i_prev_prct) || i_prct < i_prev_prct
    i_prev_prct = 0;
    S_prev = getPrctStr(i_prev_prct);
    
    fprintf('%s: %s',msg, S_prev);
end
% Print updated percentage.
if i_prct ~= i_prev_prct
    S_prev = getPrctStr(i_prev_prct);
    fprintf(getBackspaceStr(numel(S_prev)));
    
    S = getPrctStr(i_prct);
    fprintf('%s', S);
    
    i_prev_prct = i_prct;
end
% Clear percentage variable.
if i_prct == 100
    fprintf(' Done.\n');
    clear i_prev_prct;
end
end
function S = getPrctStr(prct)
S = sprintf('%d%%  %s',prct,getDotStr(prct));
if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end
function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '.';
S = ['[',S,']'];
end
function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end



