function value = coinValue(area)
    if area >= 8900 && area < 9100
        value = 0.01;
    elseif area >= 12500 && area < 12600
        value = 0.02;
    elseif area >= 15800 && area < 15900
        value = 0.05;
    elseif area >= 13400 && area < 13650
        value = 0.10;
    elseif area >= 17300 && area < 18000
        value = 0.20;
    elseif area >= 20000 && area < 21100
        value = 0.50;
    elseif area >= 19100 && area < 19500
        value = 1.00;
    else
        value = 0.00;
    end
end