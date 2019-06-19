%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                            Grupo 14                                 %%%
%%%                     João Henriques - 81633                          %%%
%%%                     Beatriz Toscan0 - 83434                         %%%
%%%                     Denis Voicu - 83443                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             PATH                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiagoPath = 'C:\Users\tiago\OneDrive\Área de Trabalho\CVI-G14\ficheiros prof\MATERIAL\database\';
denisPath = 'C:\Users\Denis\Documents\CVI\CVI-G14\ficheiros prof\MATERIAL\database\';
biaPath = 'file:///Users/beatriztoscano/CVI-G14/ficheiros%20prof/MATERIAL/database/';

%%%%% change path  here %%%%%
myPath = tiagoPath;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image = selectImage(myPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Gloabl Variables                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grayImg = rgb2gray(image);
redImage = image(:,:,1); %grab only red component
bw_gray = imbinarize(grayImg);

%Struct element -> erosion
se = strel('disk', 12);
se2 = strel('disk',8);
se3 = strel('disk', 3);
se4 = strel('disk',6);

%%%
threshold1 = floor(graythresh(redImage)*256); %200
threshold = 115;
minArea = 15;
%%%

bw_red = redImage > threshold;
bw_imerode_gray = imerode(grayImg, se);
bw_imerode_red = imerode(bw_red, se);
bw_opening = imdilate(bw_imerode_red, se2);
cleanImg = imclose(bw_gray,se4);

[lb, num] = bwlabel(bw_opening);
regionProps = regionprops(lb ,'Area','FilledImage','Centroid','Perimeter', 'BoundingBox');
inds = find([regionProps.Area]>minArea);
            
%color
white = [1.0 1.0 1.0];
red = [1.0 0 0];

 color_codes_green = [
     0     255/255    0          
     0     225/255    0          
     0     195/255    0
     0     165/255    0          
     0     135/255    0            
     0     105/255    0          
     0     75/255     0        
     0     45/255     0
     0     15/255     0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(true)
    disp('Please, choose a picture:');
    disp('     0 - Exit program;'); %done
    disp('     1 - Number of objects contained in each image;');  %done 
    disp('     2 - Visualization centroid, perimeter and area of each container;'); %done
    disp('     3 - Relative distance of the objects;'); %done
    disp('     4 - Derivative of the objects boundary (with binary image);'); %done
    disp('     5 - Objects boundary (internal and external boundary);'); %done
    disp('     6 - Ordering the object depending on area;'); %done
    disp('     7 - Ordering the object depending on perimeter;'); %done
    disp('     8 - Ordering the object depending on circularity;'); %done
    disp('     9 - Ordering the object depending on sharpness;'); %done
    disp('     10 - Compute the amount of money;'); %done
    disp('     11 - Show objects area;'); %done
    disp('     12 - Show objects perimeter;'); %done
    disp('     13 - Show objects circularity;'); %done
    disp('     14 - Show objects sharpness;'); %done
    disp('     15 - Show original image, and filters;'); %done
    disp('     16 - Histogram red image;'); %done
    disp('     17 - Histogram RGB;'); %done
    disp('     18 - Histogram grayscale;'); %done
    disp('     19 - Order by user (area);'); %done
    disp('     20 - Order by user (perimeter);'); %done
    disp('     21 - Order by user (circularity);'); %done
    disp('     22 - Order by user (sharpness);'); %done
    disp(' ');
    option = input('Your choice: ');
    switch option
        case 0
            break;
        case 1
            %Number of objects
            close all;
            figure;
            thr = floor(graythresh(redImage)*256);
            imshow(bw_opening);

            hold on
            for i=1:length(inds)
                [B,L] = bwboundaries(bw_opening,'noholes');
                hold on
                if i<=length(B)
                    boundary = B{i};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
                end
                %centroid

                plot(regionProps(inds(i)).Centroid(1),regionProps(inds(i)).Centroid(2),'g*')

            end
            str1  = sprintf('%s%d','O nº de objectos é ', length(inds));
            title(str1);
            
        case 2
            %Perimeter, area and centroid
            close all;
            figure;
            thr = floor(graythresh(redImage)*256); %200
            %a imagem binaria neste caso e a area
            imshow(bw_opening);

            hold on
            for i=1:length(inds)
%                 props = regionprops(double(regionProps(inds(i)).FilledImage),...
%                     'Orientation','MajorAxisLength','MinorAxisLength');
%                 %perimetro
%                  ellipse(props.MajorAxisLength/2,props.MinorAxisLength/2,-props.Orientation*pi/180,...
%                   regionProps(inds(i)).Centroid(1),regionProps(inds(i)).Centroid(2),'r');
                [B,L] = bwboundaries(bw_opening,'noholes');
                hold on
                if i<=length(B)
                    boundary = B{i};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
                end
                %centroid
                plot(regionProps(inds(i)).Centroid(1),regionProps(inds(i)).Centroid(2),'g*');
            
            end
            
            figure;
            subplot(1, 3, 1);
            imshow(bw_opening), hold on;
            for i=1:length(inds)
                plot(regionProps(inds(i)).Centroid(1),regionProps(inds(i)).Centroid(2),'g*');
            end
            title('Centroid', 'FontSize', 14);
            
            subplot(1, 3, 2);
            imshow(bw_opening), hold on;
            for i=1:length(inds)
                [B,L] = bwboundaries(bw_opening,'noholes');
                hold on
                if i<=length(B)
                    boundary = B{i};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
                end

            end
            title('Perimeter', 'FontSize', 14);
            
            subplot(1, 3, 3);
            imshow(bw_opening);
            title('Area', 'FontSize', 14);
            
        case 3
            %Relative distance between objects
            close all;
            figure;
            imshow(image);
            title('Select a coin');
            but = 1;

            %left mouse buttom to select coin
            while(but ==1)
                [ci, li, but] = ginput(1);
                if but == 1 %add point
                    imshow(image);
                    title('Select a coin');
                    hold on;
                    
                    %calcucla o objeto
                    object_1 = lb(li,ci);
                    % if no coin selected
                    if(object_1 == 0)
                        continue;
                    end
                    % coordinates of the centrid
                    x1 = regionProps(object_1).Centroid(1);
                    y1 = regionProps(object_1).Centroid(2);
                    %distance between the coin selected and the others
                    aux = [];
                    for i = 1: length(regionProps)
                        x2 = regionProps(i).Centroid(1);
                        y2 = regionProps(i).Centroid(2);
                        dist = sqrt((x1 - x2)^2 + (y1 - y2)^2);
                        if(dist ~= 0)
                            p = plot([x1 x2], [y1 y2], '-', 'LineWidth', 2000/dist,'DisplayName', int2str(dist));
                            aux = [aux p];
                        end   
                    end
                    legend(aux);
                end
                
                %right mouse buttom to exit
                if( but == 2)
                    break;
                end
            end
            close;

        case 4
            % Boundary with Binary Image
            close all;
            figure;
            [B,L,N,A] = bwboundaries(bw_opening);
            imshow(bw_opening); hold on;
            for k=1:length(B)
               if(~sum(A(k,:)))
                 boundary = B{k};
                 %ext Boundary
                 plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth',2);
                 for l=find(A(:,k))'
                   boundary = B{l};
                   %in Boundary
                   plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
                 end 
               end
            end
            title(' Boundary with  Opening on red image');
            legend('Ext Boundary');
                        
            for k=1:length(B)
               boundary = B{k};
               x = boundary(:, 2);
               y = boundary(:, 1);
               dx = gradient(x(:), 1); %derivada usando gradiente
               
               figure, hold on;
               %subplot(length(B)/2, 2, k);
               plot(dx), hold on;
               str1 = 'Boundary ';
               str2 = num2str(k);
               str = strcat(str1, str2);
               legend(str);
               title( ['Graph of the derivative of the object ' num2str(k) ' boundary'] )
            end
               
        case 5
            close all;
            %Boundary with no Binary Image (with inner Boundary)
            figure;
            thr = floor(graythresh(grayImg)*256); %200

            %clear operation
            se = strel('disk',3);
            bw1 = grayImg > thr;

            [B,L,N,A] = bwboundaries(bw1);
            imshow(bw1); hold on;
            for k=1:length(B)
               if(~sum(A(k,:)))
                 boundary = B{k};
                 plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
                 for l=find(A(:,k))'
                   boundary = B{l};
                   plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
                 end 
               end
            end
            title(' Boundary with a gray image');
            %legend('Ext Boundary','In Boundary');
            
        case 6
            %Ordered by ascending area
            close all;
            [obj_sorted, indsAux] = sort([regionProps.Area]);
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending area:', 'FontSize', 14);
            
        case 7
            %Ordered by ascending perimeter
            close all;
            [obj_sorted, indsAux] = sort([regionProps.Perimeter]);
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending perimeter:', 'FontSize', 14);
            
        case 8 
            %Ordered by ascending circularity
            close all;
            objs_area = [regionProps.Area];
            objs_perimeter = [regionProps.Perimeter];
            circularity = [(objs_perimeter .^ 2) ./ (4 * pi * objs_area)];
            [obj_sorted, indsAux] = sort(circularity); 
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending circularity:', 'FontSize', 14);
            
        case 9
            %Orderend by ascending sharpness
            close all;
            sharpness_array = [];
            for i=1:length(inds)
               sharpness = estimate_sharpness(double(regionProps(inds(i)).FilledImage));
               sharpness_array = [sharpness_array, sharpness]; %push sharpness
            end
            [obj_sorted, indsAux] = sort(sharpness_array);
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Orderend by ascending sharpness:', 'FontSize', 14);
            
        case 10
            %Total money
            close all;
            figure;
            imshow(image), hold on;
            total = 0;
            for i=1:length(inds)
                %area
                area = regionProps(inds(i)).Area;
                %circularity
                perimeter = regionProps(inds(i)).Perimeter;
                circularity = (perimeter .^ 2) ./ (4 * pi * area);
                if( circularity >= 0.98 && circularity < 1.03) %if is a con
                    value = coinValue(area);
                    total = total + value;
                    
                    %get centroid               
                    x1 = regionProps(inds(i)).Centroid(1);
                    y1 = regionProps(inds(i)).Centroid(2);
                    
                    moneyStr = num2str(value);
                    str = strcat (moneyStr, '€');
                    %define text properitis
                    t = text(x1 - 50,y1,str);
                    t.Color = [0 0 0];
                    s = t.FontSize;
                    t.FontSize = 30;
                    
                end
            end
            caption = sprintf('Total of money is: %.2f€', total);
            title(caption, 'FontSize', 14);
            
        case 11
            %Area of each object
            close all;
            figure;
            imshow(image), hold on; 
            for i=1:length(inds)
                %get centroid             
                x1 = regionProps(inds(i)).Centroid(1);
                y1 = regionProps(inds(i)).Centroid(2);
                str = num2str(regionProps(inds(i)).Area);
                %define text properitis
                t = text(x1 - 70,y1,str);
                t.Color = [0.0 0.0 0.0];
                s = t.FontSize;
                t.FontSize = 30;
                
                title('Area:', 'FontSize', 10);
            end
            
        case 12
            %Perimeter of each object
            close all;
            figure;
            %imshow(logical(image).*uint8(bw_opening)), hold on;
            imshow(image), hold on;
            for i=1:length(inds)
                %get centroid              
                x1 = regionProps(inds(i)).Centroid(1);
                y1 = regionProps(inds(i)).Centroid(2);
                str = num2str(regionProps(inds(i)).Perimeter);
                %define text properitis
                t = text(x1 - 70,y1,str);
                t.Color = [0 0 0];
                s = t.FontSize;
                t.FontSize = 25;
                
                title('Perimeter:', 'FontSize', 14);
            end
            
        case 13
            %Circularity of each object
            close all;
            figure;
            imshow(image), hold on;
            for i=1:length(inds)
                %get centroid
                x1 = regionProps(inds(i)).Centroid(1);
                y1 = regionProps(inds(i)).Centroid(2);
                
                area = regionProps(inds(i)).Area;
                perimeter = regionProps(inds(i)).Perimeter;
                circularity = (perimeter .^ 2) ./ (4 * pi * area);
                
                str = num2str(circularity);
                %define text properitis
                t = text(x1- 80,y1,str);
                s = t.Color;
                t.Color = [0 0 0];
                s = t.FontSize;
                t.FontSize = 25;
                
                title('Circularity:', 'FontSize', 14);
            end
            
        case 14
            %Sharpness of each object
            close all;
            figure;
            imshow(image), hold on;
            for i=1:length(inds)
                %get centroid
                x1 = regionProps(inds(i)).Centroid(1);
                y1 = regionProps(inds(i)).Centroid(2);
                
                sharpness = estimate_sharpness(double(regionProps(inds(i)).FilledImage));
                str = num2str(sharpness);
                
                %define text properitis
                t = text(x1 - 80,y1,str);
                %s = t.Color;
                t.Color = [0 0 0];
                s = t.FontSize;
                t.FontSize = 25;
                
                title('Sharpness:', 'FontSize', 14);
            end
            
        case 15
            %Initial img and filter
            close all;
            %show original image
            fig1 = figure(1);
            subplot(3,3,1);
            imshow(image);
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            title('Original image', 'FontSize', 14);

            %convert image to grayscale
            subplot(3, 3, 2);
            imshow(grayImg);
            title('Grayscale image', 'FontSize', 14);

            %red image
            subplot(3, 3, 3);
            imshow(redImage);
            title('Grayscale red image', 'FontSize', 14);

            %BW Red image
            subplot(3, 3, 4);
            imshow(bw_red);
            title('BW Red image', 'FontSize', 14);

            %binarize image
            subplot(3, 3, 5);
            imshow(bw_gray);
            title('Binarized image', 'FontSize', 14);

            %Close operation
            subplot(3, 3, 6)
            imshow(cleanImg);
            title('Close operation', 'FontSize', 14);

            %imerode red image
            subplot(3, 3, 7);
            imshow(bw_imerode_red);
            title('imerode red image', 'FontSize', 14);

            %imerode gray image
            subplot(3, 3, 8);
            imshow(bw_imerode_gray);
            title('imerode gray image', 'FontSize', 14);

            %Opening red image
            subplot(3, 3, 9);
            imshow(bw_opening);
            title('Opening red image', 'FontSize', 14);
            
        case 16
            %histogram red img
            close all;
            figure, hold on;
            subplot(1, 2, 1);
            imshow(redImage);
            subplot(1, 2, 2);
            imhist(redImage); hold on
            thr = floor(graythresh(redImage)*256); %200
            plot(thr, 0, 'r.', 'markersize', 35)
            title('Histogram Red image');
            
        case 17
            %Histogram RGB
            close all;            
            R = imhist(image(:,:,1));
            G = imhist(image(:,:,2));
            B = imhist(image(:,:,3));
            figure, hold on;
            plot(R,'r');
            plot(G,'g');
            plot(B,'b');
            title('RGB Histogram');
            legend('Red channel','Green channel','Blue channel');
            
        case 18
            %Histogram greyscale
            close all;
            figure, hold on;
            imhist(grayImg);
            title('Histogram of gray Image');
            
        case 19
            %Order by area after coin selected
            close all;
            figure;
            imshow(image); hold on;
            title('Select one coin:');
            n = 0; but = 1;

            %1 = tecla esquerda do rato
            while(but ==1)
                [ci, li, but] = ginput(1);
                if but == 1 %add point
                    n = n+1;
                    cp(n) = ci;
                    lp(n) = li;
                    plot(ci, li, 'r.', 'MarkerSize', 8); drawnow;
                end
                if( n == 1)
                    break;
                end
            end
            close;

            %calcucla o objeto
            obj1 = lb(uint16(lp(1)),uint16(cp(1)));
            %calcula as coordenadas os centroids
            x1 = regionProps(obj1).Centroid(1);
            y1 = regionProps(obj1).Centroid(2);
            %calcula a area do objeto selecionado
            obj1_area = regionProps(obj1).Area;

            %calculas as outras areas
            %objs_area = [regionProps.Area] %inclui a do obejto selecionado
            %objs_centroid = [regionProps.Centroid] %inclui a do obejto selecionado
            for i=1:length(inds)
                if regionProps(inds(i)).Area ~= obj1_area && x1 ~= regionProps(inds(i)).Centroid(1) && y1 ~= regionProps(inds(i)).Centroid(2)
                    dif = abs(obj1_area - regionProps(inds(i)).Area);
                    regionProps(inds(i)).Dif_areas = dif;
                    regionProps(inds(i)).Dif_areas;
                else
                    regionProps(inds(i)).Dif_areas = 0;
                end
            end
            
            objs_dif_areas = [regionProps.Dif_areas];
            [obj_sorted, indsAux] = sort(objs_dif_areas); 
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending area similarity (from user selection):', 'FontSize', 14), hold on;
            %desenha a nova imagem, marca os centros e uma linha entre eles
            plot(x1, y1, 'r.-', 'MarkerSize', 8);
            legend('Selected object');
            
        case 20
            %Ordered perimeter after coin selected
            close all;
            figure;
            imshow(image); hold on;
            title('Selecet one coin:');
            n = 0; but = 1;

            %1 = tecla esquerda do rato
            while(but ==1)
                [ci, li, but] = ginput(1);
                if but == 1 %add point
                    n = n+1;
                    cp(n) = ci;
                    lp(n) = li;
                    plot(ci, li, 'r.', 'MarkerSize', 8); drawnow;
                end
                if( n == 1)
                    break;
                end
            end
            close;

            %calcucla o objeto
            obj1 = lb(uint16(lp(1)),uint16(cp(1)));
            %calcula as coordenadas os centroids
            x1 = regionProps(obj1).Centroid(1);
            y1 = regionProps(obj1).Centroid(2);
            %calcula a area do objeto selecionado
            obj1_perimeter = regionProps(obj1).Perimeter;

            for i=1:length(inds)
                if regionProps(inds(i)).Perimeter ~= obj1_perimeter && x1 ~= regionProps(inds(i)).Centroid(1) && y1 ~= regionProps(inds(i)).Centroid(2)
                    dif = abs(obj1_perimeter - regionProps(inds(i)).Perimeter);
                    regionProps(inds(i)).Dif_perimeter = dif;
                    regionProps(inds(i)).Dif_perimeter;
                else
                    regionProps(inds(i)).Dif_perimeter = 0;
                end
            end
            
            objs_dif_perimeter = [regionProps.Dif_perimeter];
            [obj_sorted, indsAux] = sort(objs_dif_perimeter); 
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending perimeter similarity (from user selection):', 'FontSize', 14), hold on;
            %desenha a nova imagem, marca os centros e uma linha entre eles
            plot(x1, y1, 'r.-', 'MarkerSize', 8);
            legend('Selected object');
            
        case 21
            %Ordered by circularity after coin selected
            close all;
            figure;
            imshow(image); hold on;
            title('Selecet one coin:');
            n = 0; but = 1;

            %1 = tecla esquerda do rato
            while(but ==1)
                [ci, li, but] = ginput(1);
                if but == 1 %add point
                    n = n+1;
                    cp(n) = ci;
                    lp(n) = li;
                    plot(ci, li, 'r.', 'MarkerSize', 8); drawnow;
                end
                if( n == 1)
                    break;
                end
            end
            close;

            %calcucla o objeto
            obj1 = lb(uint16(lp(1)),uint16(cp(1)));
            %calcula as coordenadas os centroids
            x1 = regionProps(obj1).Centroid(1);
            y1 = regionProps(obj1).Centroid(2);
            %calcula a area do objeto selecionado
            obj1_perimeter = regionProps(obj1).Perimeter;
            obj1_circularity = (obj1_perimeter .^ 2) ./ (4 * pi * regionProps(obj1).Area);
            
            for i=1:length(inds)
                if regionProps(inds(i)).Perimeter ~= obj1_perimeter && x1 ~= regionProps(inds(i)).Centroid(1) && y1 ~= regionProps(inds(i)).Centroid(2)
                    obj_circularity = (regionProps(inds(i)).Perimeter .^ 2) ./ (4 * pi * regionProps(inds(i)).Area);
                    
                    dif = abs(obj1_circularity - obj_circularity);
                    regionProps(inds(i)).Dif_circularity = dif;
                    regionProps(inds(i)).Dif_circularity;
                else
                    regionProps(inds(i)).Dif_circularity = 0;
                end
            end
            
            objs_dif_circularity = [regionProps.Dif_circularity];
            [obj_sorted, indsAux] = sort(objs_dif_circularity); 
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending circularity similarity (from user selection):', 'FontSize', 14), hold on;
            %desenha a nova imagem, marca os centros e uma linha entre eles
            plot(x1, y1, 'r.-', 'MarkerSize', 8);
            legend('Selected object');
            
        case 22
            %Ordered by sharpness after coin selected
            figure;
            imshow(image); hold on;
            title('Selecet one coin:');
            n = 0; but = 1;

            %1 = tecla esquerda do rato
            while(but ==1)
                [ci, li, but] = ginput(1);
                if but == 1 %add point
                    n = n+1;
                    cp(n) = ci;
                    lp(n) = li;
                    plot(ci, li, 'r.', 'MarkerSize', 8); drawnow;
                end
                if( n == 1)
                    break;
                end
            end
            close;

            %calcucla o objeto
            obj1 = lb(uint16(lp(1)),uint16(cp(1)));
            %calcula as coordenadas os centroids
            x1 = regionProps(obj1).Centroid(1);
            y1 = regionProps(obj1).Centroid(2);
            %calcula a area do objeto selecionado
            obj1_area= regionProps(obj1).Area;
            obj1_sharpness = estimate_sharpness(double(regionProps(obj1).FilledImage));
            for i=1:length(inds)
                if regionProps(inds(i)).Area ~= obj1_area && x1 ~= regionProps(inds(i)).Centroid(1) && y1 ~= regionProps(inds(i)).Centroid(2)
                    obj_sharpness = estimate_sharpness(double(regionProps(inds(i)).FilledImage));
                    
                    dif = abs(obj1_sharpness - obj_sharpness);
                    regionProps(inds(i)).Dif_sharpness = dif;
                    regionProps(inds(i)).Dif_sharpness;
                else
                    regionProps(inds(i)).Dif_sharpness = 0;
                end
            end
            
            objs_dif_sharpness = [regionProps.Dif_sharpness];
            [obj_sorted, indsAux] = sort(objs_dif_sharpness); 
            
            figure, hold on;
            aux = zeros(size(lb)); %lb aux
            for i = 1:length(inds)
               %get the row and the columns of the objects
               [row, column] = find(lb == inds(i));
               matrix_aux = [row column];
               for j=1:length(matrix_aux)
                   %fill the matrix with the row and collomns, by area
                   aux(matrix_aux(j, 1), matrix_aux(j, 2)) = find(indsAux == i);
               end
            end
            imshow(label2rgb(aux, color_codes_green, [0 0 0])), hold on;
            title('Ordered by ascending circularity sharpness (from user selection):', 'FontSize', 14), hold on;
            
            %desenha a nova imagem, marca os centros e uma linha entre eles
            plot(x1, y1, 'r.-', 'MarkerSize', 8);
            legend('Selected object');
        otherwise
            disp('Please, insert a valid image.');
    end
end
