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
cont5 = 0;
cont10 = 0;
cont15 = 0;
cont20 = 0;
cont25 = 0;
cont30 = 0;
cont35 = 0;
cont40 = 0;
cont45 = 0;
cont50 = 0;
cont55 = 0;
cont60 = 0;
cont65 = 0;
cont70 = 0;
cont75 = 0;
cont80 = 0;
cont85 = 0;
cont90 = 0;
cont95 = 0;
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
            %optimizaçao
%             %Trajectories Map
%             subplot(2, 2, 3);
%             viscircles(regionProps(inds(j)).Centroid, 1, 'Color', [0 1 0]); %[color_aux(j) color_aux(j+1) color_aux(j+2)]
% 
%             %heatmap calculation
%              x = uint8(regionProps(inds(j)).Centroid(1));
%              y = uint8(regionProps(inds(j)).Centroid(2));
%              %preenche as posicoes passadas
%              heatMap_aux(x, y) = heatMap_aux(x, y) + 1;
%              
%              %calcula o maximo
%              if heatMap_aux(x, y) > heatMap_max
%                  heatMap_max = heatMap_aux(x, y);
%              end
%             axis([0 576 0 768]);
%             set(gca, 'YDir', 'reverse');
%             title('Trajectories Map', 'FontSize', 14);
%             subplot(2,2,1);
        end
        
        subplot(2,2,1);
        tp = 0;
        fn = 0;
        fp = 0;
        contOverlaps = 0;
        iou = 0;
        
        %metrics
        [lin col] = find(file(:,1) == i);
        for a = 1 : length(bbAlg)
            contFp = 0;
            for p = 1 : length(lin)
                %get BB from file
                rec = rectangle('Position', [[file(lin(p),3) file(lin(p),4)] [file(lin(p),5) file(lin(p),6)]], 'EdgeColor',  [1 0 0], 'linewidth', 2);
                % overlap between alg BB and red BB
                overlap = rectint(rec.Position, bbAlg(a).Position);
                gtArea = file(lin(p),5) * file(lin(p),6);
                intr = (overlap / gtArea) * 100;
                %for the iou metric
                if overlap > 0
                    contOverlaps = contOverlaps + 1;
                    union = gtArea + bbAlg(a).Position(3) * bbAlg(a).Position(4) - overlap;
                    iou = iou + overlap / union;
                    if contOverlaps == length(lin)
                        iou = iou / contOverlaps;
                        if iou > 0
                            cont0 = cont0 + 1;
                        end
                        if iou > 0.05
                            cont5 = cont5 + 1;
                        end
                        if iou > 0.10
                            cont10 = cont10 + 1;
                        end
                        if iou > 0.15
                            cont15 = cont15 + 1;
                        end
                        if iou > 0.20
                            cont20 = cont20 + 1;
                        end
                        if iou > 0.25
                            cont25 = cont25 + 1;
                        end
                        if iou > 0.30
                            cont30 = cont30 + 1;
                        end
                        if iou > 0.35
                            cont35 = cont35 + 1;
                        end
                        if iou > 0.40
                            cont40 = cont40 + 1;
                        end
                        if iou > 0.45
                            cont45 = cont45 + 1;
                        end
                        if iou > 0.50
                            cont50 = cont50 + 1;
                        end
                        if iou > 0.55
                            cont55 = cont55 + 1;
                        end
                        if iou > 0.60
                            cont60 = cont60 + 1;
                        end
                        if iou > 0.65
                            cont65 = cont65 + 1;
                        end
                        if iou > 0.70
                            cont70 = cont70 + 1;
                        end
                        if iou > 0.75
                            cont75 = cont75 + 1;
                        end
                        if iou > 0.80
                            cont80 = cont80 + 1;
                        end
                        if iou > 0.85
                            cont85 = cont85 + 1;
                        end
                        if iou > 0.90
                            cont90 = cont90 + 1;
                        end
                        if iou > 0.95
                            cont95 = cont95 + 1;
                        end
                        if iou == 1
                            cont1 = cont1 + 1;
                        end
                        figure(3);
                        perc0 = cont0 / i * 100;
                        perc5 = cont5 / i * 100;
                        perc10 = cont10 / i * 100;
                        perc15 = cont15 / i * 100;
                        perc20 = cont20 / i * 100;
                        perc25 = cont25 / i * 100;
                        perc30 = cont30 / i * 100;
                        perc35 = cont35 / i * 100;
                        perc40 = cont40 / i * 100;
                        perc45 = cont45 / i * 100;
                        perc50 = cont50 / i * 100;
                        perc55 = cont55 / i * 100;
                        perc60 = cont60 / i * 100;
                        perc65 = cont65 / i * 100;
                        perc70 = cont70 / i * 100;
                        perc75 = cont75 / i * 100;
                        perc80 = cont80 / i * 100;
                        perc85 = cont85 / i * 100;
                        perc90 = cont90 / i * 100;
                        perc95 = cont95 / i * 100;
                        perc1 = cont1 / i * 100;
                        x = 0:5:100;
                        bar(x, [perc0 perc5 perc10 perc15 perc20 perc25 perc30 perc35 perc40 perc45 perc50 ...
                                perc55 perc60 perc65 perc70 perc75 perc80 perc85 perc90 perc95 perc1], 0.4);
                        axis([-1 100 0 100]);
                        title('IOU');
                        xlabel('Threshold');
                        ylabel('Percentage%');
                        set(0,'CurrentFigure',h);
                    end
                end
                if overlap > 0 && intr <= 60
                    fp = fp + 1;
                elseif overlap > 0 && intr >= 60
                    difArea = bbAlg(a).Position(3) * bbAlg(a).Position(4) - gtArea;
                    if difArea < 0
                        difArea = difArea * -1;
                    end
                    if intr == 100 && difArea < 1000
                        tp = tp + 1;
                    elseif intr < 100 && difArea < 1000
                        tp = tp + 1;
                    else
                        fp = fp + 1;
                    end    
                elseif overlap == 0 
                   contFp = contFp + 1; 
                end
            end
            vermelhos = length(lin);
            if contFp == vermelhos
               fp = fp + 1; 
            end
        end      
    end
    vermelhos = length(lin);
    fn = vermelhos - tp;
    
    precision = tp /(tp + fp);
    recall = tp / ( tp + fn);
    
    pVec = cat(1,pVec,precision);
    rVec = cat(1,rVec,recall);
    frameVec = cat(1,frameVec,i);
    
    drawnow;
    hold off;
    
    figure(2);
    subplot(2,2,1);
    plot(rVec, pVec, 'O');
    title('Recall-Precision', 'FontSize', 14);
    xlabel('Recall');
    ylabel('Precision');
    
    subplot(2,2,3);
    if length(rVec) > 2
        f = fit(rVec,pVec,'smoothingspline');
        plot(f,rVec,pVec);
        title('Recall-Precision Fitted Curve', 'FontSize', 14);
        set(gca,'NextPlot','replacechildren') ;
    end
    
    
    subplot(1,2,2)
    plot(frameVec, pVec);   
    hold on;
    plot(frameVec,rVec,'g');
    title('Recall-Precision by Frame', 'FontSize', 14);
    xlabel('Frame');
    ylabel('Precision and Recall');
    hold off;
    
    
    
    set(0,'CurrentFigure',h);
    %dinamycall draw of pedestrians move 
    [m n] =  size(pos_x); %size(pos_x == pos_y)
    if m > 20
        aux_x = pos_x(end-19:end);
        aux_y = pos_y(end-19:end);

        subplot(2,2,2), axis ij;
        scatter(aux_x, aux_y, 40, linspace(1, 10, 20), 'filled');
        axis([0 768 0 576]);
        set(gca, 'YDir', 'reverse');
        title('Dinamycall draw', 'FontSize', 14);
    end
    
    
    %Trajectories Map
    subplot(2, 2, 3);
    if regnum
        for j = 1:regnum
            h = viscircles(regionProps(inds(j)).Centroid, 1, 'Color', [0 1 0]); %[color_aux(j) color_aux(j+1) color_aux(j+2)]
           
            %heatmap calculation
             y = uint16(regionProps(inds(j)).Centroid(1));
             x = uint16(regionProps(inds(j)).Centroid(2));
             %preenche as posicoes passadas
             heatMap_aux(x, y) = heatMap_aux(x, y) + 10;
             
             %calcula o maximo
             if heatMap_aux(x, y) > heatMap_max
                 heatMap_max = heatMap_aux(x, y);
             end

        end
    end
    axis([0 768 0 576]);
    set(gca, 'YDir', 'reverse');
    title('Trajectories Map', 'FontSize', 14);
    
    %heatmap
    subplot(2, 2, 4);
    scatter_kde(pos_x, pos_y, 'filled', 'MarkerSize', 50);
    axis([0 768 0 576]);
    set(gca, 'YDir', 'reverse');
    title('HeatMap', 'FontSize', 14);
      
end

%%nao funciona
% [grad,im] = colorGradient([0 0 1], [1 0 0], heatMap_max + 1);
% subplot(, 2, 1);
% for j = 1:length(pos_x)
%     value = heatMap_aux(uint16(pos_x(j)),uint16(pos_y(j)));
%     viscircles([pos_x(j), pos_y(j)], 1, 'Color', grad((value + 1), :));
% end
% axis([0 576 0 768]);
% set(gca, 'YDir', 'reverse');
% title('hhf', 'FontSize', 14);
         
% [grad,im] = colorGradient([0 0 1], [1 0 0], heatMap_max + 1);
% subplot(2, 2, 4);
% title('Heat map', 'FontSize', 14);
% [m n] =  size(heatMap_aux);
% for j = 1:m
%     for k = 1:n  
%         value = heatMap_aux(j, k);
%         h = viscircles([j k], 1, 'Color', grad((value + 1), :));
%     end
% end
figure;

%subplot(2, 2, 2);
imshow(bckg), hold on;
scatter_kde(pos_x, pos_y, 'filled', 'MarkerSize', 50);
title('Trajectory', 'FontSize', 14);

heat_map = HeatMap(heatMap_aux, 'Symmetric', 'false');     
heat_map.Colormap = 'jet';