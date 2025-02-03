clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% control points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% defined by mouse
% hold on
% flag = 1;
% k = 1;
% while(flag == 1)
%     [x, y, flag] = ginput(1);
%     if(flag == 1)
%        cp(k,1) = x;
%        cp(k,2) = y;
%        plot(cp(k,1), cp(k,2), 'k--o'); hold on;
%        if(k>=2)
%            line([cp(k-1,1) cp(k,1)], [cp(k-1,2) cp(k,2)]);
%        end
%        k = k + 1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%control points sampled by function
a = [5 12 14 16 -5];
x = [];
y = [];
co = 1;
for i = 1:0.1:5
    x(co) = i;
    y(co) = a(5) * i^5 + a(4) * i^4 + a(3) * i^3 + a(2) * i^2 + a(1) * i;
    plot(x(co), y(co), 'k--o'); hold on;
    co = co + 1;
end
ga = (co-1) / 10;
for i = 1:10
    num = floor(ga*i);
    cp(i,1) = x(num);
    cp(i,2) = y(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 bc = Bspline(cp, 3);
 plot(bc(:,1), bc(:,2), 'r-');
 
function bsplineCurve = Bspline(cp, degree)
num_cp = size(cp, 1);
num_knot = num_cp + degree + 1;
kv = linspace(0, 1, num_knot);
num_gp = 1000;
gap = 1 / num_gp;
bsplineCurve = zeros(num_gp,2);
t = 1;
for u=0:gap:1
    sum = 0;
    BasicF = zeros(1, num_cp);
    for i=1:num_cp
        BasicF(1,i) = BasicFunction(i-1, degree, u, kv);
        sum = sum + BasicF(1,i);
    end
            
    BasicF = BasicF / sum;
    bsplineCurve(t,:) = BasicF * cp;
    t = t + 1;
end
end
    
function BF = BasicFunction(i, degree, u, kv)

if(degree == 0)
    if(kv(i+1) <= u && kv(i + 2) > u)
       BF = 1;
    else
       BF = 0;
    end
else
    d1 = kv(i + degree + 1) - kv(i + 1);
    d2 = kv(i + degree + 2) - kv(i + 2);
    if(d1 == 0)
        d1 = 1;
    end
    
    if(d2 == 0)
        d2 = 1;
    end
    
    BF = (u - kv(i+1)) / d1 * BasicFunction(i, degree - 1, u, kv) + (kv(i + degree + 2) - u) / d2 * BasicFunction(i+1, degree - 1, u, kv);
end
end

     