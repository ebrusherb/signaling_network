global e eps;
eps = .5;
e = 1;
bvals=linspace(0,4*eps*e,30);
tvals=linspace(0,3*eps,30);
threshold = meshgrid(bvals,tvals);
for i=1:length(bvals)
    for j=1:length(tvals)
        threshold(j,i)=e*(tvals(j)+eps);
    end
end
figure
surf(tvals,bvals,transpose(threshold))
hold on

for i=1:length(bvals)
    for j=1:length(tvals)
        if tvals(j)>eps
            threshold(j,i)=e*(tvals(j)-eps);
        else threshold(j,i)=abs(e*(2*eps*e-bvals(i))/(2*eps*e+bvals(i))*(tvals(j)-eps));
        end
    end
end

surf(tvals,bvals,transpose(threshold))
