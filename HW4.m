clear all;
close all;
clc;

%% Video 1: Study Hall
video = [];
v = VideoReader('IMG_3132.MOV');
while hasFrame(v)
    frame = readFrame(v);
    frame = rgb2gray(frame);
    frame = reshape (frame, [], 1);
    video = [video, frame];
end
video = reshape(video, [1920,1080,80]);
video = imresize(video,.25);
video = double(video);
video = reshape(video, [480*270,80]);
v1 = video(:,1:end-1);
v2 = video(:,2:end);
[U, Sigma, V] = svd(v1, 'econ');
%%
r=1;
Sr = Sigma(1:r, 1:r);
Ur = U(:, 1:r);
Vr = V(:, 1:r);
Stilde = Ur'*v2*Vr*diag(1./diag(Sr));
[eV, D] = eig(Stilde);
mu = diag(D);
omega = log(mu);
Phi = v2*Vr/Sr*eV;

y0 = Phi\video(:,1);
v_modes = zeros(r,length(v1(1,:)));
for i = 1:length(v1(1,:))
    v_modes(:,i) = (y0.*exp(omega*i));
end
v_dmd = Phi*v_modes;
v_dmd = abs(v_dmd);
v_sparse = v1 - v_dmd;

residual_matrix = v_sparse.*(v_sparse < 0);
v_dmd = residual_matrix + abs(v_dmd);
v_sparse = v_sparse - residual_matrix;
uvid = reshape(v_dmd, [480, 270, 79]);

figure(1)

for i = 1:12
    subplot(3,4,i)
    vidtype = floor((i-1)/4);
    timeframe = mod(i,4);
    if timeframe == 0
        timeframe = 4;
    end
    if vidtype == 0
        temp = video(:,timeframe*15);
    elseif vidtype == 1
        temp = v_dmd(:,timeframe*15); 
    elseif vidtype == 2
        temp = v_sparse(:,timeframe*15); 
    end
    temp = reshape(temp, 480, 270);
    imagesc(temp);
    colormap(gray);
    axis off;
end

%% Video 2: Cars Running Downhill

video = [];
v = VideoReader('IMG_3154.MOV');
while hasFrame(v)
    frame = readFrame(v);
    frame = rgb2gray(frame);
    frame = reshape (frame, [], 1);
    video = [video, frame];
end
video = reshape(video, [1920,1080,147]);
video = imresize(video,.25);
video = double(video);
video = reshape(video, [480*270,147]);
v1 = video(:,1:end-1);
v2 = video(:,2:end);
[U, Sigma, V] = svd(v1, 'econ');
%%
r=3;
Sr = Sigma(1:r, 1:r);
Ur = U(:, 1:r);
Vr = V(:, 1:r);
Stilde = Ur'*v2*Vr*diag(1./diag(Sr));
[eV, D] = eig(Stilde);
mu = diag(D);
omega = log(mu);
Phi = v2*Vr/Sr*eV;

y0 = Phi\video(:,1);
v_modes = zeros(r,length(v1(1,:)));
for i = 1:length(v1(1,:))
    v_modes(:,i) = (y0.*exp(omega*i));
end
v_dmd = Phi*v_modes;
v_dmd = abs(v_dmd);
v_sparse = v1 - v_dmd;

residual_matrix = v_sparse.*(v_sparse < 0);
v_dmd = residual_matrix + abs(v_dmd);
v_sparse = v_sparse - residual_matrix;
uvid = reshape(v_dmd, [480, 270, 146]);

figure(2)

for i = 1:12
    subplot(3,4,i)
    vidtype = floor((i-1)/4);
    timeframe = mod(i,4);
    if timeframe == 0
        timeframe = 4;
    end
    if vidtype == 0
        temp = video(:,timeframe*15);
    elseif vidtype == 1
        temp = v_dmd(:,timeframe*15); 
    elseif vidtype == 2
        temp = v_sparse(:,timeframe*15); 
    end
    temp = reshape(temp, 480, 270);
    imagesc(temp);
    colormap(gray);
    axis off;
end

%% Video 3: Eating an Egg

video = [];
v = VideoReader('52359343966__6F3FF01D-F7C5-473A-A049-FD070AC65BD6.MOV');
while hasFrame(v)
    frame = readFrame(v);
    frame = rgb2gray(frame);
    frame = reshape (frame, [], 1);
    video = [video, frame];
end
video = reshape(video, [1280,720,95]);
video = imresize(video,.25);
video = double(video);
video = reshape(video, [320*180,95]);
v1 = video(:,1:end-1);
v2 = video(:,2:end);
[U, Sigma, V] = svd(v1, 'econ');
%%
r=1;
Sr = Sigma(1:r, 1:r);
Ur = U(:, 1:r);
Vr = V(:, 1:r);
Stilde = Ur'*v2*Vr*diag(1./diag(Sr));
[eV, D] = eig(Stilde);
mu = diag(D);
omega = log(mu);
Phi = v2*Vr/Sr*eV;

y0 = Phi\video(:,1);
v_modes = zeros(r,length(v1(1,:)));
for i = 1:length(v1(1,:))
    v_modes(:,i) = (y0.*exp(omega*i));
end
v_dmd = Phi*v_modes;
v_dmd = abs(v_dmd);
v_sparse = v1 - v_dmd;

residual_matrix = v_sparse.*(v_sparse < 0);
v_dmd = residual_matrix + abs(v_dmd);
v_sparse = v_sparse - residual_matrix;
uvid = reshape(v_dmd, [320, 180, 94]);

figure(3);

for i = 1:12
    subplot(3,4,i)
    vidtype = floor((i-1)/4);
    timeframe = mod(i,4);
    if timeframe == 0
        timeframe = 4;
    end
    if vidtype == 0
        temp = video(:,timeframe*15);
    elseif vidtype == 1
        temp = v_dmd(:,timeframe*15); 
    elseif vidtype == 2
        temp = v_sparse(:,timeframe*15); 
    end
    temp = reshape(temp, 320, 180);
    imagesc(temp);
    colormap(gray);
    axis off;
end




