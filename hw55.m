f = @(t, y)(y^2 + y)/t;
a = 1;
b = 3;
hmax = .5;
hmin = .02;
TOL = 10^-4;
i = 1;
t(i) = a;
w(i) = -2;





h = hmax;
FLAG = 1;

%Step 2:
while FLAG == 1
    %step 3
    
    k1 = h * f(t(i), w(i));
    k2 = h * f(t(i)+h/4 , w(i) + k1/4);
    k3 = h * f(t(i)+(3*h/8), w(i) + (3*k1/32) + (9*k2/32));
    k4 = h * f(t(i) + (12*h/13), w(i)+(1932*k1/2197) - (7200*k2/2197) + (7296*k3/2197));
    k5 = h * f(t(i)+h, w(i)+(439*k1/216)-8*k2 +(3680*k3/513) -(845*k4/4104));
    k6 = h * f(t(i) + h/2 , w(i) - (8*k1/27) +2*k2 - (3544*k3/2565) + (1859*k4/4104) - 11*k5/40);

    %step 4
    g = (k1/360)-(128*k3/4275)-(2197*k4/75240)+(k5/50)+(2*k6/55);
    R = (1/h) * abs(g);

    %step 5
    if R <= TOL
        %step 6
        t(i+1) = t(i)+h; %approx accepted, step forward
        w(i+1) = w(i) + (25*k1/216)+(1408*k3/2565)+(2197*k4/4104)-k5/5;
        i = i +1;
       
        
        
    end %step 7

    %step 8 
    d = .84*(TOL/R)^.25;

    %step 9 
    if d <= .1
        h = .1*h;

    elseif d >= 4
            h = 4*h;
    else
        h = d*h;
    end

    %step 10
    if h > hmax
        h = hmax;
    end

    %step 11
    if t(i) >= b
        FLAG = 0;
    elseif t(i) + h > b
        h = b - t(i);
    elseif h < hmin
        FLAG = 0;
    end
  %  i = i+1;

end
figure
plot(t, w);
title('Runge-Kutta-Fehlberg 2c');
xlabel('1 < t < 3');
ylabel('y values');
%% 4d
f2 = @(t2, y2)(-y2 + t2*y2^.5);

a2 = 2;
b2 = 4;
hmax2 = .5;
hmin2 = .05;
TOL2 = 10^-6;

%Step 1:
i = 2;
t2(i) = a2;
w2(i) = 2;
h2 = hmax2;
FLAG = 1;

%Step 2:
while FLAG == 1
    %step 3
    
    k1 = h2 * f2(t2(i), w2(i));
    k2 = h2 * f2(t2(i)+h2/4 , w2(i) + k1/4);
    k3 = h2 * f2(t2(i)+(3*h2/8), w2(i) + (3*k1/32) + (9*k2/32));
    k4 = h2 * f2(t2(i) + (12*h2/13), w2(i)+(1932*k1/2197) - (7200*k2/2197) + (7296*k3/2197));
    k5 = h2 * f2(t2(i)+h2, w2(i)+(439*k1/216)-8*k2 +(3680*k3/513) -(845*k4/4104));
    k6 = h2 * f2(t2(i) + h2/2 , w2(i) - (8*k1/27) +2*k2 - (3544*k3/2565) + (1859*k4/4104) - 11*k5/40);

    %step 4
    g = (k1/360)-(128*k3/4275)-(2197*k4/75240)+(k5/50)+(2*k6/55);
    R = (1/h2) * abs(g);

    %step 5
    if R <= TOL2
        %step 6
        t2(i+1) = t2(i)+h2; %approx accepted, step forward
        w2(i+1) = w2(i) + (25*k1/216)+(1408*k3/2565)+(2197*k4/4104)-k5/5;
        i = i +1;
       
        
        
    end %step 7

    %step 8 
    d = .84*(TOL2/R)^.25;

    %step 9 
    if d <= .1
        h2 = .1*h2;

    elseif d >= 4
            h2 = 4*h2;
    else
        h2 = d*h2;
    end

    %step 10
    if h2 > hmax2
        h2 = hmax2;
    end

    %step 11
    if t2(i) >= b2
        FLAG = 0;
    elseif t2(i) + h2 > b2
        h2 = b2 - t2(i);
    elseif h2 < hmin2
        FLAG = 0;
    end
  %  i = i+1;

end

%actual solution
for j = 2:length(t2)
    yexact(j) = (t2(j)-2+(2^.5)*exp(1)*exp(-t2(j)/2))^2;
end

figure
plot(t2, w2,' r+-','DisplayName','RKF');
hold on
title('Runge-Kutta-Fehlberg 4d');
xlabel('2 < t < 4');
ylabel('y values');
plot(t2, yexact,'k*-- ','DisplayName','exact');
legend
hold off
%% 5a, 5b
f3 = @(t, y1)(2*10^-6)*(100000-y1)*y1;
a3 = 1;
b3 = 30;
i = 1;
t3(i) = 1;
w3(i) = 1000;
hmax3 = .5;
hmin3 = .05; 
TOL3 = 10^-6;
h3 = hmax3;
FLAG = 1;

while FLAG == 1
    %step 3
    
    k1 = h3 * f3(t3(i), w3(i));
    k2 = h3 * f3(t3(i)+h3/4 , w3(i) + k1/4);
    k3 = h3 * f3(t3(i)+(3*h3/8), w3(i) + (3*k1/32) + (9*k2/32));
    k4 = h3 * f3(t3(i) + (12*h3/13), w3(i)+(1932*k1/2197) - (7200*k2/2197) + (7296*k3/2197));
    k5 = h3 * f3(t3(i)+h3, w3(i)+(439*k1/216)-8*k2 +(3680*k3/513) -(845*k4/4104));
    k6 = h3 * f3(t3(i) + h3/2 , w3(i) - (8*k1/27) +2*k2 - (3544*k3/2565) + (1859*k4/4104) - 11*k5/40);

    %step 4
    g = (k1/360)-(128*k3/4275)-(2197*k4/75240)+(k5/50)+(2*k6/55);
    R = (1/h3) * abs(g);

    %step 5
    if R <= TOL3
        %step 6
        t3(i+1) = t3(i)+h3; %approx accepted, step forward
        w3(i+1) = w3(i) + (25*k1/216)+(1408*k3/2565)+(2197*k4/4104)-k5/5;
        i = i +1;
       
        
        
    end %step 7

    %step 8 
    d = .84*(TOL3/R)^.25;

    %step 9 
    if d <= .1
        h3 = .1*h3;

    elseif d >= 4
            h3 = 4*h3;
    else
        h3 = d*h3;
    end

    %step 10
    if h3 > hmax3
        h3 = hmax3;
    end

    %step 11
    if t3(i) >= b3
        FLAG = 0;
    elseif t3(i) + h3 > b3
        h3 = b3 - t3(i);
    elseif h3 < hmin3
        FLAG = 0;
    end
  %  i = i+1;

end
m = 100000;
k = 2*10^-6;
c = 9.9*10^-4;
%actual solution, found by hand
for i=1:length(t3)
    yexact3(i) = ((1/m) + c*exp(-m*k*t3(i)))^-1;
end

%estimated value for t = 30
disp('Part 5a: ');
disp(w3(end));

%part b
figure
plot(t3, w3, 'k*--', 'DisplayName','RKF45');
hold on
plot(t3, yexact3, 'r+-', 'DisplayName', 'Actual Solution');
title('5a and 5b');
ylabel('y values');
xlabel('0 < t < 30');
legend
hold off






