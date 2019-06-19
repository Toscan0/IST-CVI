function image = selectImage(path)
    img1 = 'Moedas1.jpg';
    img2 = 'Moedas2.jpg';
    img3 = 'Moedas3.jpg';
    img4 = 'Moedas4.jpg';
    
    while(true)
        disp('Please, choose a picture:');
        disp('     0 - Exit program;');
        disp('     1 - Moedas1.jpg;');
        disp('     2 - Moedas2.jpg;');
        disp('     3 - Moedas3.jpg;');
        disp('     4 - Moedas4.jpg;');
        disp(' ');
        option = input('Your choice: ');
        switch option
            case 0
                break;
            case 1
                s = strcat(path, img1);
                image = imread(s);
                break;
            case 2
                s = strcat(path, img2);
                image = imread(s);
                break;
            case 3
                s = strcat(path, img3);
                image = imread(s);
                break;
            case 4
                s = strcat(path, img4);
                image = imread(s);
                break;
            otherwise
                disp('Please, insert a valid image.');
        end
    end
end