%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                            Grupo 14                                 %%%
%%%                     Joao Henriques - 81633                          %%%
%%%                     Beatriz Toscan0 - 83434                         %%%
%%%                     Denis Voicu - 83443                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all;
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             PATH                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiagoPath = 'C:\Users\tiago\OneDrive\Área de Trabalho\CVI-G14\';
denisPath = 'C:\Users\Denis\Documents\CVI\CVI-G14\';
biaPath = 'file:///Users/beatriztoscano/CVI-G14\';

%%%%% change path  here %%%%%
myPath = strcat(tiagoPath, 'proj2\3DMOT2015\train\PETS09-S2L1\img1\');
bckgPath = strcat(tiagoPath, 'proj2\BackGround.jpg');
txtPath = strcat(tiagoPath, 'proj2\3DMOT2015\train\PETS09-S2L1\gt\gt.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se = strel('disk', 3);
thr = 60;
minArea = 100;

nameSize = 6; %6 digits -> name of file
str = ['%s%.' num2str(nameSize) 'd.%s'];

step = 1;
nFrame = 795; %number of files

file = importdata(txtPath,',');

%contadores das frames maiores que uj threshold
cont0 = 0;
cont25 = 0;
cont50 = 0;
cont75 = 0;
cont1 = 0;

try
    %try to read the backgorund imagem if calculated before
    bckg = imread(bckgPath);
catch ME
    %if background image doesnt exist then calculate the image 
    for i = 1 : 1 : nFrame
        str1 = sprintf(str, myPath, i, 'jpg');
        img = imread(str1);
        fotoSum(:, :, :, i) = img;  
    end
    bckg = median(fotoSum, 4);
    %saves image
    imwrite(bckg, bckgPath);
end

pos_x = [];
pos_y = [];
pVec = [];
rVec = [];
frameVec = [];

%[rows, columns] = size(bckg);
heatMap_aux = zeros(580, 770);
heatMap_max = 0;
       
for i=1 : 1 : nFrame
    bbAlg = [];
    h =figure(1);
    imgfr = imread(sprintf(str, myPath, i, 'jpg'));
    subplot(2,2,1);
    imshow(imgfr);
    title_aux = sprintf('%d Frame', i);
    title(title_aux);
    hold on;
    imgdif = ...
        (abs(double(bckg(:,:,1))-double(imgfr(:,:,1)))>thr) | ... %R
        (abs(double(bckg(:,:,2))-double(imgfr(:,:,2)))>thr) | ... %G
        (abs(double(bckg(:,:,3))-double(imgfr(:,:,3)))>thr);      %B
    
    bw = imclose(imgdif, se);
    
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb, 'Area', 'Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);

    if regnum
        for j = 1:regnum
            [lin col] = find(lb == inds(j));
            upLPoint = min([lin col]);
            dWindow = max([lin col]) - upLPoint + 1;
            %Bw_aux(lin, col) = 1;
            %bounding box
            r = rectangle('Position', [fliplr(upLPoint) fliplr(dWindow)], 'EdgeColor',  [1 1 0], 'linewidth', 2);
            yellowLine = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color',[1 1 0]);
            redLine = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color',[1 0 0]);
            legend([yellowLine, redLine], 'Bounding box', 'Ground truth');
            pos_x = cat(1, pos_x, regionProps(inds(j)).Centroid(1));
            pos_y = cat(1, pos_y, regionProps(inds(j)).Centroid(2));
            %add rectangle to vector
            bbAlg = cat(1,bbAlg,r);
        end
        
        subplot(2,2,1);
        tpF = 0;
        fnF = 0;
        fpF = 0;
        contOverlaps = 0;
        iou = 0;
        
        %metrics
        [lin col] = find(file(:,1) == i);
        if i == 1
            for v = 0:1:100
                for a = 1 : length(bbAlg)
                    contfpF = 0;
                    for p = 1 : length(lin)
                        %get BB from file
                        rec = rectangle('Position', [[file(lin(p),3) file(lin(p),4)] [file(lin(p),5) file(lin(p),6)]], 'EdgeColor',  [1 0 0], 'linewidth', 2);
                        % overlap between alg BB and red BB
                        overlap = rectint(rec.Position, bbAlg(a).Position);
                        gtArea = file(lin(p),5) * file(lin(p),6);
                        intr = (overlap / gtArea) * 100;
                        if overlap > 0 && intr <= v
                            fpF = fpF + 1;
                        elseif overlap > 0 && intr >= v
                            difArea = bbAlg(a).Position(3) * bbAlg(a).Position(4) - gtArea;
                            if difArea < 0
                                difArea = difArea * -1;
                            end
                            if intr == 100 && difArea < 1000
                                tpF = tpF + 1;
                            elseif intr < 100 && difArea < 1000
                                tpF = tpF + 1;
                            else
                                fpF = fpF + 1;
                            end    
                        elseif overlap == 0 
                           contfpF = contfpF + 1; 
                        end
                    end
                    vermelhos = length(lin);
                    if contfpF == vermelhos
                       fpF = fpF + 1; 
                    end
                end 
                vermelhos = length(lin);
                fnF = vermelhos - tpF;

                precisionF = tpF /(tpF + fpF);
                recallF = tpF / ( tpF + fnF);

                pVec = cat(1,pVec,precisionF);
                rVec = cat(1,rVec,recallF);
            end
        end
          
    end
    
    drawnow;
    hold off;
    
    figure(2);
    subplot(2,2,1);
    plot(rVec, pVec);
   
end