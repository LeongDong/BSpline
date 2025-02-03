clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sampling from known 3D function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count = 1; %the total number of control points
% county = 1; %the number of control points along y-axis
% for i=-3:0.5:3
%     for j = -3:0.5:3
%         cp(count,1) = j;
%         cp(count,2) = i;
%         cp(count,3) = (j*j - 2*j) * exp(-j*j-i*i-i*j);
%         count = count + 1;
%     end
%     county = county + 1;
% end
% county = county -1;
% countx = (count - 1) / county;%the number of control points along x-axis
% sc = 1;
% 
% figure;
% plot3(cp(1,:),cp(2,:), cp(3,:),'bo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%sampling from gray image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgpath = '/media/liangdong/C14D581BDA18EBFA/111.png';
img=imread(imgpath);
img = rgb2gray(img);
img = imresize(img,0.1);
[m,n] = size(img);
count = 1;
county = 1;
for i=1:floor(m/32):m
    for j=1:floor(n/32):n
        cp(count,1) = j;
        cp(count,2) = i;
        cp(count,3) = img(i,j);
        count = count + 1;
    end
    county = county + 1;
end
county = county -1;
countx = (count - 1) / county;%the number of control points along x-axis
sc = 1;

figure;
plot3(cp(:,1),cp(:,2), cp(:,3), '.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
bc = Bspline(cp, 3, countx, county);
% plot3(bc(:,1), bc(:,2), bc(:,3), '.');
[X,Y,Z]=griddata(bc(:,1), bc(:,2), bc(:,3),linspace(min(bc(:,1)),max(bc(:,1)))',linspace(min(bc(:,2)),max(bc(:,2))),'v4');
surf(X,Y,Z);

function bsplineSurface = Bspline(cp, degree, countx, county) %cp: control points;
num_cp = size(cp, 1);
num_knot_y = county + degree + 1; %the number of knots along y-axis
num_knot_x = countx + degree + 1; %the number of knots along x-axis
kvy = linspace(0, 1, num_knot_y); %knots along y-axis
kvx = linspace(0, 1, num_knot_x); %knots along y-axis
num_gp = 50; 
gap = 1 / num_gp; % sampling gap, smaller more precision;otherwise.
bsplineSurface = zeros(num_gp,3);
t = 1;
for uy=0:gap:1
    for ux=0:gap:1
        sum = 0;
        cpt = 1;
        BasicF = zeros(1, num_cp);
        for i=1:county
            for j=1:countx
                BasicF(1,cpt) = BasicFunction(i-1, degree, uy, kvy) * BasicFunction(j-1, degree, ux, kvx);
                sum = sum + BasicF(1,cpt);
                cpt = cpt + 1;
            end
        end
            
        BasicF = BasicF / (sum+0.01); %re-normalized weights
        bsplineSurface(t,:) = BasicF * cp; %calculate the coordinate in (ux,uy), which is distance weighted sum of all control point, weights are from basis functions.
        t = t + 1;
    end
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